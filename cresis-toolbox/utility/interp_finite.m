function vals = interp_finite(vals,default_val,interp_method)
% vals = interp_finite(vals,default_val,interp_method)
%
% vals: a vector of numbers
% default_val: if all vals are ~isfinite, then the whole vector will
%   be set to this value
% interp_method: interpolation passed to interp1
%
% vals = all isfinite elements remain unchanged, all ~isfinite values will
%   interpolated with from the isfinite values using this scheme:
%   1. values on the end are interpolated using nearest neighbor
%   2. values in the middle will be linearly interpolated
%   3. if no isfinite values exist, the whole vector is set to default_val
%
% Author: John Paden

if ~exist('interp_method','var')
  interp_method = 'linear';
end

good_mask = isfinite(vals);

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
