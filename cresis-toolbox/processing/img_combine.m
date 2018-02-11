function [newData, newTime] = img_combine(combine)
% [newData, newTime] = img_combine(combine)
%
% img_combine: Blends and combines together echogram datafiles into a
%  combined datafile according to the parameter set contained in the
%  'combine' object.
%
% Surface tracking image combine
%  combine.img_comb(1): time after surface return where
%    combine will happen
%  combine.img_comb(2): minimum time that combine will occur
%  combine.img_comb(3): guard time which specifies how
%    many seconds at the end of img1 will not be used... this is
%    important because the last samples of img1 will have low signal
%    power and blurred because they will only have captured a portion
%    of the chirp energy (typically this will be set to something
%    close to the pulse duration for img1)
%  combine.img_comb(4-6, 7-9, etc.): same fields as above
%    except between images 2 and 3, 3 and 4, etc.
%
% Authors: John Paden, Victor Berger
%
% See also: run_update_img_combine.m, update_img_combine.m

% combine.Data, combine.Time             => already loaded data
% combine.appendData, combine.appendTime => new data to append
% newData, newTime                       => combined result

if ~isempty(combine.img_comb_weights)
  combine.appendData = combine.appendData.*10.^(combine.img_comb_weights(idx)/10);
end

% Interpolate image N onto already loaded data (assumption is that image
% N-1 always comes before image N)
dt = combine.Time(2)-combine.Time(1);
newTime = (combine.Time(1) : dt : combine.appendTime(end)).';
combine.appendData = interp1(combine.appendTime,combine.appendData,newTime,'linear',0);
combine.imb_comb_surf = interp_finite(combine.imb_comb_surf,0);

% Determine guard at end of image 1 that will not be used
blend_bins = combine.img_comb_bins;
guard_bins = 1 + round(combine.img_comb((combine.idx-2)*3+3)/dt);
max_good_time = length(combine.Time)*ones(1,size(combine.Data,2));

% First row of img_bins indicates the start of the blend-region
if ~(combine.imb_comb_surf == 0)
  img_bins = round(interp1(newTime, 1:length(newTime), ...
    max(min(combine.imb_comb_surf * combine.img_comb_mult, ...
    combine.imb_comb_surf + combine.img_comb((combine.idx-2)*3 + 1)), ...
    combine.img_comb((combine.idx-2)*3 + 2)), 'linear', 'extrap'));
else
  img_bins = max_good_time-guard_bins;
end

% Check to make sure requested time is inside window and just
% force the combination bin to occur at the second to last bin
%   img_bins outside the img1 time window will be NaN due to interp1
%   img_bins inside the img1 time window may still be larger than
%     the guard allows
invalid_rlines = find(isnan(img_bins) ...
  | img_bins > max_good_time-guard_bins-blend_bins);
img_bins(invalid_rlines) = max_good_time(invalid_rlines)-guard_bins-blend_bins;

% Second row of img_bins indicates the end of the blend-region
img_bins(2,:) = img_bins(1,:) + 1 + blend_bins;

% Estimate difference
if strcmpi(combine.img_comb_weights_mode,'auto')
  newData = zeros(size(combine.appendData),'single');
  difference = NaN*zeros(1,size(newData,2));
  for rline = 1:size(newData,2)
    trans_bins = img_bins(1,rline)+1:img_bins(2,rline);
    if trans_bins <= size(combine.appendData,1)
      difference(rline) = mean(combine.Data(trans_bins,rline) ./ combine.appendData(trans_bins,rline));
    end
  end
  difference = nanmedian(difference);
  fprintf('  Difference: %.1f\n', lp(difference));
  difference_report(frm) = difference;
else
  difference = 10^(-0/10); % For debugging
end

% Combine images
newData = zeros(size(combine.appendData),'single');
for rline = 1:size(newData,2)
  trans_bins = img_bins(1,rline)+1:img_bins(2,rline);
  weights = 0.5+0.5*cos(pi*linspace(0,1,length(trans_bins)).');
  if trans_bins <= size(combine.appendData,1)
    newData(:,rline) = [combine.Data(1:img_bins(1,rline),rline); ...
      weights.*combine.Data(trans_bins,rline) ...
      + difference*(1-weights).*combine.appendData(trans_bins,rline); ...
      difference*combine.appendData(img_bins(2,rline)+1:end,rline)];
  else
    newData(:,rline) = combine.Data(1:size(newData,1),rline);
  end
end
