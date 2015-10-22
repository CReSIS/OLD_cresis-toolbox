function mean_val = mean_without_outliers(vals, std_threshold)
% mean_val = mean_without_outliers(vals, std_threshold)
%
% vals is an array of values, operates on first non-singleton dimension
%   NaN values are ignored
%
% Uses median filter and standard deviation based threshold to cull bad
% values and take the mean of the good values only.
%
% Author: John Paden

if ~exist('std_threshold','var')
  std_threshold = 2;
end

median_vals = median(vals);
std_vals = std(vals);

dim = find(size(vals) > 1,1);
if isempty(dim)
  mean_val = vals;
end

repmat_size = ones(size(size(vals)));
repmat_size(dim) = size(vals,dim);

std_vals(std_vals == 0) = inf;
good_mask = abs(vals - repmat(median_vals,repmat_size)) < repmat(std_vals,repmat_size) * std_threshold;

vals(~good_mask) = NaN;

% 2011 Matlab does not support nanmean
%mean_val = nanmean(vals,dim);
vals_nan = isnan(vals);
vals(vals_nan) = 0;
mean_val = sum(vals,dim) ./ sum(~vals_nan,dim);

end
