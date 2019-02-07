function vals = interp_finite(vals,default_val,interp_method,interp_vals_fh)
% vals = interp_finite(vals,default_val,interp_method,interp_vals_fh)
%
% vals: a vector of numbers
% default_val: if all vals are ~isfinite, then the whole vector will
%   be set to this value
% interp_method: interpolation method passed to interp1 (default is
%   'linear')
% interp_vals_fh: function handle that returns true for values that will
% not be interpolated. Default is @isfinite. Other common examples are:
%   @isfinite       % interpolate NaN, inf, and -inf
%   @(x) ~isnan(x)  % interpolate only NaN
%   @(x) x~=inf     % interpolate only inf
% 
%
% vals = all isfinite elements remain unchanged, all ~isfinite values are
%   interpolated from the isfinite values using this scheme:
%   1. values on the end are interpolated using nearest neighbor
%   2. values in the middle will be linearly interpolated
%   3. if no isfinite values exist, the whole vector is set to default_val
%
% Author: John Paden

if ~exist('interp_method','var') || isempty(interp_method)
  interp_method = 'linear';
end

if ~exist('interp_vals_fh','var') || isempty(interp_vals_fh)
  interp_vals_fh = @isfinite;
end

good_mask = interp_vals_fh(vals);

%% For bad values at the beginning and end, use nearest neighbor interpolation
first_good = find(good_mask,1);
if isempty(first_good)
  % The whole vector is ~isfinite, so set to default_val and return
  vals(:) = default_val;
  return;
end
vals(1:first_good) = vals(first_good);
good_mask(1:first_good) = 1;

last_good = find(good_mask,1,'last');
vals(last_good:end) = vals(last_good);
good_mask(last_good:end) = 1;

%% Use linear interpolation for everything in between
vals(~good_mask) = interp1(find(good_mask),vals(good_mask),find(~good_mask),interp_method);

end
