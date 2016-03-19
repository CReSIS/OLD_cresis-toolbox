function surface = tracker_threshold(data,surf)
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
%  .noise_override: double scalar, use this value for the noise estimate
%     rather than the value estimated from noise_rng
%  .threshold = double scalar, relative threshold above noise estimate
%    in log10 scale
%  .sidelobe	= double scalar, sidelobe value that specifies the minimum
%    threshold value as max_value - surf.sidelobe. In log10 scale.
%  .max_diff = maximum deviation of the surface from the inital surface (in
%    range bins), set to inf to not have a limit
%  .filter_len = an along track boxcar filter is applied, must be a
%    positive odd integer (e.g. 1, 3, 5, etc)
%  .min_bin = the minimum row (range bin) data that the tracker will look at
%  .init.method = 'medfilt' or 'snake' or 'dem'
%     .init.search_rng = required for tracker_snake_simple (e.g. [-120:120])
%     .init.medfilt = required for medfilt1 (e.g. 101)
%  .dem = surface dem used by .init.method == 'dem'
%  .detrend = polynomial order of detrend function (polyfit to power data)
%     zero or lower or leaving empty/undefined does no detrending
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
if ~isfield(surf,'noise_rng')
  surf.noise_rng = [0 -inf -1];
end
if ~isfield(surf,'noise_override')
  surf.noise_override = [];
end
if ~isfield(surf,'search_rng')
  surf.search_rng = 0;
end
if ~isfield(surf,'threshold')
  surf.threshold = 17;
end
if ~isfield(surf,'sidelobe')
  surf.sidelobe = 13;
end
if ~isfield(surf,'max_diff')
  surf.max_diff = inf;
end
if ~isfield(surf,'filter_len')
  surf.filter_len = 1;
end
if mod(surf.filter_len,2) == 0
  error('Surface filter length must be odd');
end

if ~isfield(surf,'min_bin') || isempty(surf.min_bin)
  surf.min_bin = 1;
end

if ~isfield(surf,'max_bin') || isempty(surf.max_bin) || ~isfinite(surf.max_bin)
  surf.max_bin = size(data,1);
end

if ~isfield(surf,'detrend') || isempty(surf.detrend)
  surf.detrend = 0;
end

if ~isfield(surf,'dem')
  % Constrain maximum search with this expected surface
  surf.dem = [];
end

%% Clip data according to min and max bins that surface is allowed to be in
data = data(surf.min_bin:surf.max_bin,:);

%% Filter data and convert to log domain
surf.normalize_each_bin = false;
if ~surf.normalize_each_bin
  surf_data = fir_dec(data,ones(1,surf.filter_len(1))/surf.filter_len(1),1);
  if length(surf.filter_len) == 2
    surf_data = fir_dec(surf_data.',ones(1,surf.filter_len(2))/surf.filter_len(2),1).';
  end
  surf_data = lp(surf_data,1);
else
  surf_data = fir_dec(data,ones(1,surf.filter_len)/surf.filter_len,1);
  surf_data = surf_data - repmat(mean(surf_data,2),[1 size(surf_data,2)]);
  surf_data = lp(surf_data);
end

%% Detrend the data
if surf.detrend > 0
  poly_x = (-size(surf_data,1)/2+(1:size(surf_data,1))).';
  mean_power = lp(mean(10.^(surf_data/10),2));
  good_mask = isfinite(mean_power);
  p = polyfit(poly_x(good_mask),mean_power(good_mask),surf.detrend);
  surf_data = surf_data - repmat(polyval(p,poly_x),[1 size(surf_data,2)]);
end

%% Determine the threshold value using range bins specified by noise_rng
median_mdata = zeros(1,size(surf_data,2));
new_surface_max = zeros(1,size(surf_data,2));
new_surface_max_val = zeros(1,size(surf_data,2));
for rline=1:size(surf_data,2)
  [new_surface_max_val(rline),new_surface_max(rline)] = max(surf_data(:,rline));
  % Remove range bins that are zero
  good_mask = surf_data(:,rline) ~= 0;
  % Remove range bins at start and end of record
  good_mask(1:min(length(good_mask),find(good_mask,1)-1+surf.noise_rng(1))) = 0;
  good_mask(max(1,find(good_mask,1,'last')+1-surf.noise_rng(1)):end) = 0;
  % Grab the specified range of bins
  good_mask(1:min(length(good_mask),max(0,new_surface_max(rline)+surf.noise_rng(2)))) = 0;
  good_mask(max(1,new_surface_max(rline)+surf.noise_rng(3)):end) = 0;
  if any(good_mask)
    median_mdata(rline) = median(surf_data(good_mask,rline));
  end
end

%% Create Initial Surface
if isfield(surf,'init')
  if strcmp(surf.init.method,'snake')
    new_surface_max = tracker_snake_simple(surf_data,surf.init);
  elseif strcmp(surf.init.method,'medfilt')
    new_surface_max = medfilt1(new_surface_max,surf.init.medfilt);
  elseif strcmp(surf.init.method,'dem') & length(surf.dem) == size(surf_data,2)
    new_surface_max = surf.dem - surf.min_bin + 1;
  else
    error('Unsupported surface init method or surf.dem is wrong length');
  end
end

if all(median_mdata==0)
  warning('Insufficient information to generate threshold since all values are zero.');
  surface = new_surface_max;
  return
end
if isempty(surf.noise_override)
  THRESHOLD = surf.threshold + median(median_mdata(median_mdata ~= 0));
else
  THRESHOLD = surf.threshold + surf.noise_override;
end

%% Perform thresholding
surface = NaN*zeros(1,size(surf_data,2));
for rline=1:size(surf_data,2)
  rbins = max(1,round(new_surface_max(rline) - surf.max_diff)) : min(size(surf_data,1),round(new_surface_max(rline) + surf.max_diff));
  threshold_rline = max(max(surf_data(:,rline))-surf.sidelobe,THRESHOLD);
  thresh_bin = find(surf_data(rbins,rline) > threshold_rline,1);
  if ~isempty(thresh_bin)
    surface(rline) = rbins(1) - 1 + thresh_bin;
  end
end

%% Do not permit thresholded surface from exceeding max_diff from initial surface
surface(abs(surface - new_surface_max) > surf.max_diff) = NaN;
surface = interp_finite(surface,0);
if all(surface == 0)
  surface = new_surface_max;
end

if (length(surf.search_rng) > 1 || surf.search_rng ~= 0)
  if 1
    % Find the next peak after the threshold
    for rline = 1:size(surf_data,2)
      search_bins = round(surface(rline)) + surf.search_rng;
      search_bins = search_bins(find(search_bins >= 1 & search_bins <= size(surf_data,1)));
      %offset = find(diff(surf_data(search_bins,rline)) < 0,1);
      [~,offset] = max(surf_data(search_bins,rline));
      if isempty(offset)
        surface(rline) = search_bins(end);
      else
        surface(rline) = search_bins(offset);
      end
    end
  end
end
surface = surface + surf.min_bin - 1;

if 0
  % Debug plot
  surf_plot = surface - surf.min_bin;
  figure(1); clf;
  imagesc(surf_data);
  hold on;
  plot(new_surface_max,'b');
  plot(surf_plot,'k','LineWidth',2);
  ylims = [min(min(surf_plot),min(new_surface_max)) max(max(surf_plot),max(new_surface_max))];
  ylims = ylims + diff(ylims)*[-0.1 0.1];
  ylim(ylims)
  cc = caxis;
  caxis([THRESHOLD-20 cc(2)]);
  colormap(1-gray(256));
  keyboard
end

end
