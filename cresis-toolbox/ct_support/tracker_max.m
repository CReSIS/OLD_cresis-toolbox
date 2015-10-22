function surface = tracker_max(data,surf)
% surface = tracker_max(data,surf)
%
% data = 2D matrix of nonnegative linear power values
% surf = structure controlling operating of the tracker
%  .min_bin = positive integer, the minimum range bin to search for the
%    max in;
%  .search_rng = range of values to search for exceeding a threshold
%    value from the max (set to 0 if you just want to return the max)
%    Usually of the form [-120:120] or [-5:5] to search +/-120 and +/-5
%    range bins from the max.
%  .threshold = multiplying factor in log10 scale to the max value.
%    E.g. a value of -C finds the first value in the search range which
%    is above max_value-C.
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
if ~isfield(surf,'min_bin')
  surf.min_bin = 1;
end
if ~isfield(surf,'max_bin') | isempty(surf.max_bin)
  surf.max_bin = size(data,1);
end
if ~isfield(surf,'search_rng')
  surf.search_rng = 0;
end
if ~isfield(surf,'threshold')
  surf.threshold = 5;
end

[max_values surface]= max(data(surf.min_bin:surf.max_bin,:));
surface = surface + surf.min_bin - 1;

if length(surf.search_rng) > 1 || surf.search_rng ~= 0
  % Following code does a simple threshold search
  for rline = 1:size(data,2)
    search_bins = surface(rline) + surf.search_rng;
    search_bins = search_bins(search_bins >= surf.min_bin & search_bins <= size(data,1));
    threshold_idx = find(data(search_bins,rline) > max(data(search_bins,rline)) .* 10.^(-surf.threshold/10),1);
    surface(rline) = search_bins(threshold_idx);
  end
end

if 0
  % Debug plot
  figure(1); clf;
  imagesc(lp(data));
  hold on;
  plot(surface,'b');
  ylims = [min(surface) max(surface)];
  ylims = ylims + diff(ylims)*[-0.1 0.1];
  ylim(ylims)
  colormap(1-gray(256));
end

end
