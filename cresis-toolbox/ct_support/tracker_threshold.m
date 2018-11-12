function [surface,quality] = tracker_threshold(data,surf)
% surface = tracker_threshold(data,surf)
%
% data = Nt by Nx 2D matrix of nonnegative linear power values
% surf = structure controlling operating of the tracker
%  .noise_rng = [100 -700 -300];
%    Specifies the region relative to the peak power to estimate the noise
%    power from. It is a 1 by 3 integer array.
%    first element: All data points that are zero at the beginning of the
%      record will be ignored in the noise calculation. This specifies a
%      buffer beyond this in sample bins that will also be ignored.
%    second and third elements: Specifies a range of bins relative to the
%      max bin that will be used to estimate the noise (usually these
%      are negative such as [-50 -10] so that the noise estimate uses data
%      before the peak).
%  .threshold_noise_dB: double scalar, use this value for the noise estimate
%     rather than the value estimated from noise_rng
%  .threshold = double scalar, relative threshold above noise estimate
%    in log10 scale. Default is 17.
%  .sidelobe	= double scalar, sidelobe value that specifies the minimum
%    threshold value as max_value - surf.sidelobe. In log10 scale. Default
%    is 13.
%  .max_diff = maximum deviation of the surface from the inital surface (in
%    range bins), set to inf to not have a limit
%  .filter_len = if defined and not empty, it should be a 2 element vector
%    that specifies the size of a 2D boxcar filter to apply to the data
%    before tracking. The first dimension is cross track and the second
%    dimension is along-track. Each dimension must be a positive odd
%    integer (e.g. [1 7], [3 3], [3 7], etc).
%  .min_bin = the minimum row (range bin) data that the tracker will look at
%  .init.method = 'medfilt' or 'snake' or 'dem'
%     .init.search_rng = required for tracker_snake_simple (e.g. [-120:120])
%     .init.medfilt = required for medfilt1 (e.g. 101)
%  .dem = surface dem used by .init.method == 'dem'
%  .detrend = polynomial order of detrend function (polyfit to power data)
%     zero or lower or leaving empty/undefined does no detrending
%  .search_rng: search range around threshold bin to find a max value
%  .threshold_rng: threshold for a given range line uses this range of
%     range values around the current range line to estimate the threshold.
%     Set to inf to include all range lines in estimate. Default is inf.
%
% surface = 1 by size(data,2) vector of range bins where the surface
%   was tracked at
%
% See also: get_heights.m, update_surface_with_tracker.m,
% tracker_snake_simple.m, tracker_max.m, tracker_threshold.m
%
% Author: John Paden

if ~exist('surf','var')
  surf = struct();
end
if ~isfield(surf,'threshold_noise_rng') || isempty(surf.threshold_noise_rng)
  surf.threshold_noise_rng = [0 -inf -1];
end
if ~isfield(surf,'threshold_noise_dB') || isempty(surf.threshold_noise_dB)
  surf.threshold_noise_dB = [];
end
if ~isfield(surf,'threshold') || isempty(surf.threshold)
  surf.threshold = 17;
end
if ~isfield(surf,'threshold_rng') || isempty(surf.threshold_rng)
  surf.threshold_rng = inf;
end


%% Determine the threshold value using range bins specified by noise_rng
median_mdata = zeros(1,size(data,2));
for rline=1:size(data,2)
  % Remove range bins that are zero
  good_mask = data(:,rline) ~= 0;
  % Remove range bins at start and end of record
  good_mask(1:min(length(good_mask),find(good_mask,1)-1+surf.threshold_noise_rng(1))) = 0;
  good_mask(max(1,find(good_mask,1,'last')+1-surf.threshold_noise_rng(1)):end) = 0;
  % Grab the specified range of bins
  good_mask(1:min(length(good_mask),max(0,surf.dem(rline)+surf.threshold_noise_rng(2)))) = 0;
  good_mask(max(1,surf.dem(rline)+surf.threshold_noise_rng(3)):end) = 0;
  if any(good_mask)
    median_mdata(rline) = nanmedian(data(good_mask,rline));
  end
end

if ~isempty(surf.threshold_noise_dB) && all(median_mdata==0)
  warning('Insufficient information to generate threshold since all values are zero.');
  surface = surf.dem;
  return
end
if isempty(surf.threshold_noise_dB)
  if isfinite(surf.threshold_rng)
    THRESHOLD = NaN;
  else
    THRESHOLD = surf.threshold + median(median_mdata(~isnan(median_mdata)));
  end
else
  THRESHOLD = surf.threshold + surf.threshold_noise_dB;
end

%% Perform thresholding
surface = NaN*zeros(1,size(data,2));
quality = 3*ones(1,size(data,2));
for rline=1:size(data,2)
  rbins = max(1,round(surf.dem(rline) - surf.init.max_diff)) : min(size(data,1),round(surf.dem(rline) + surf.init.max_diff));

  % Fast Sidelobe Check
  if isnan(THRESHOLD)
    mask = median_mdata ~= 0 & (1:size(data,2)) - rline <= surf.threshold_rng;
    RLINE_THRESHOLD = surf.threshold + median(median_mdata(mask));
    threshold_rline = RLINE_THRESHOLD;
  else
    threshold_rline = THRESHOLD;
  end
  
  thresh_bin = find(data(rbins,rline) > threshold_rline,1);
  if ~isempty(thresh_bin)
    quality(rline) = 1;
    surface(rline) = rbins(1) - 1 + thresh_bin;
  end
end

if 0
  % Debug plot
  figure(1); clf;
  imagesc(data);
  hold on;
  plot(surf.dem,'b');
  plot(find(quality==1), surface(quality==1),'g.','LineWidth',2);
  plot(find(quality==3), surface(quality==3),'r.','LineWidth',2);
  cc = caxis;
  colorbar;
  if THRESHOLD < cc(2)
    caxis([THRESHOLD THRESHOLD+3]);
  else
    title(sprintf('ERROR: no value in image is larger than THRESHOLD %g.', THRESHOLD));
  end
  colormap(1-gray(256));
  keyboard
end
