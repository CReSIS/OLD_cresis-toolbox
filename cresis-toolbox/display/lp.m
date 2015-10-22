function data = lp(data,method)
% function data = lp(data,method)
%
% Log power plot
%
% data = n-d array of data
% method = Optional scalar. 1 does 10*log10, 2 does 20*log10(abs())
%   Default is:
%   If real data:
%     10*log10(data)
%   If complex data:
%     20*log10(abs(data))
%
% imagesc(lp(data));
%
% See also: inc_filt.m

if ~exist('method','var') || isempty(method)
  if isreal(data) && all(data(~isnan(data)) >= 0)
    % Data is probably magnitude or power (real and non-negative)
    method = 1;
  else
    % Data is either real w/ negative values or complex
    method = 2;
  end
end

if method == 1
  data = 10*log10(abs(data));
else
  data = 20*log10(abs(data));
end

return;

