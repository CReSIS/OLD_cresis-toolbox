function vals = interp_finite(vals,default_val,interp_fh,interp_mask_fh,extrap_mode)
% vals = interp_finite(vals,default_val,interp_fh,interp_mask_fh,extrap_mode)
%
% Inputs
% =========================================================================
% vals: a vector of numbers
%
% default_val: scalar. If the mask created by interp_mask_fh returns no
% valid data, then the whole vector will be set to this value. If this
% value
%
% interp_fh: interpolation function handle that takes three arguments
% similar to Matlab's interp1. To support legacy interface, if this
% function is a string, then the default is interp1 with the method
% specified in the string (see interp1.m for supported methods). Otherwise
% the default is interp1 with 'linear' method. Note that edge samples that
% would require extrapolation are always interpolated with nearest neighbor
% in extrap_mode==0. Common examples:
%   @interp1 % linear interpolation (default)
%   'linear' % linear interpolation with interp1 (default)
%   'nearest' % nearest neighbor interpolation with interp1
%   @(xi,yi,xq) interp1(xi,yi,xq,'spline') % spline interpolation
%   @gps_interp1 % gps_interp1 interpolation for polar data
%
% interp_mask_fh: function handle that returns true for values that will be
% used to interpolate the values that return false. Default is @isfinite.
% Common examples are:
%   @isfinite       % interpolate NaN, inf, and -inf (default)
%   @(x) ~isnan(x)  % interpolate only NaN
%   @(x) x~=inf     % interpolate only inf
%
% extrap_mode: string containing 'nearest' (default) or 'interp'
%
% Outputs
% =========================================================================
% vals: vector of the same size as vals. All elements which had a true mask
% remain unchanged, all elements that had a false mask are interpolated
% from the isfinite values using this scheme:
%   1. values on the end are interpolated using nearest neighbor
%   2. values in the middle will be interpolated with interp_fh
%   3. if no isfinite values exist, the whole vector is set to default_val
%
% Author: John Paden

if ~exist('interp_fh','var') || isempty(interp_fh)
  interp_fh = @interp1;
elseif ischar(interp_fh)
  interp_fh = @(x,y,z) interp1(x,y,z,interp_fh);
end

if ~exist('interp_mask_fh','var') || isempty(interp_mask_fh)
  interp_mask_fh = @isfinite;
end

if ~exist('extrap_mode','var') || isempty(extrap_mode)
  extrap_mode = 0;
elseif strcmpi(extrap_mode,'interp')
  extrap_mode = 1;
else
  extrap_mode = 0;
end

good_mask = interp_mask_fh(vals);

if extrap_mode == 0
  %% For bad values at the beginning and end, use nearest neighbor interpolation
  for col = 1:numel(vals)/size(vals,1)
    first_good = find(good_mask(:,col),2);
    if isempty(first_good)
      % The whole vector is ~isfinite, so set to default_val and return
      if ~exist('default_val','var') || isempty(default_val)
        error('No valid samples and no default_val was given.');
      end
      vals(:) = default_val;
      
    elseif length(first_good) == 1
      vals(:) = vals(first_good);
      
    elseif ~extrap_mode
      % Extrapolation uses nearest neighbor
      first_good = first_good(1);
      vals(1:first_good,col) = vals(first_good,col);
      good_mask(1:first_good,col) = 1;
      
      last_good = find(good_mask(:,col),1,'last');
      vals(last_good:end,col) = vals(last_good,col);
      good_mask(last_good:end,col) = 1;
      
      %% Use linear interpolation for everything in between
      vals(~good_mask(:,col),col) = interp_fh(find(good_mask(:,col)),vals(good_mask(:,col)),find(~good_mask(:,col)));
    end
  end
  
else
  %% Interpolation and Extrapolation use interpolation function
  vals(~good_mask) = interp_fh(find(good_mask),vals(good_mask),find(~good_mask));
  
end
