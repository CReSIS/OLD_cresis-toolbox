function data = inc_filt(data,rows,cols,method)
% data = inc_filt(data,rows,cols,method)
%
% Incoherently averages data in row and column dimensions. Detects the
% type of data (complex or real) and then converts to log scale.
% It does not decimate the data.
%
% INPUTS:
%  method: optional parameter, default is 1 for real data, 2 for complex
%    1: takes abs(data)
%    2: takes abs(data)^2
%
% Example:
%
% imagesc(incfilt(data,4,10))
%
% imagesc(incfilt(data,4,10,2))
%
% Author: John Paden
%
% See also: lp.m


if ~exist('method','var') || isempty(method)
  if isreal(data) && all(data(:) > 0)
    method = 1;
  else
    method = 2;
  end
end

if method == 1
  data = 10*log10(filter2(ones(rows,cols)/(rows*cols),abs(data)));
else
  data = 10*log10(filter2(ones(rows,cols)/(rows*cols),abs(data).^2));
end

return;

