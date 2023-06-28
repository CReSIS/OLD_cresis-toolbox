function yi = interp1_robust(x,y,xi,varargin)
% yi = interp1_robust(x,y,xi,varargin)
%
% Wrapper for interp1 that handles cases where x is empty or length 1. If
% empty, NaN are returned. If length 1, then 'nearest' interpolation is
% used.
%
% Examples:
% interp1_robust([],[],[2.5 4])
% ans =
%    NaN   NaN
% interp1_robust([2],[5],[2.5 4])
% ans =
%      5     5
% interp1_robust([2 3],[5 6],[2.5 4],'linear','extrap')
% ans =
%    5.500000000000000   7.000000000000000
%
% Author: John Paden

if isempty(x)
  yi = nan(size(xi));
elseif length(x) == 1
  yi = repmat(y,size(xi));
else
  if isempty(varargin)
    yi = interp1(x,y,xi);
  else
    yi = interp1(x,y,xi,varargin{:});
  end
end
