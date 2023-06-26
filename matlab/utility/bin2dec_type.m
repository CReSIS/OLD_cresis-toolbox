function x=bin2dec_type(s, type_fh)
% x=bin2dec_type(s, type_fh)
%
% Fix for Matlab's bin2dec to work with uint64 and arbitrarily large
% numbers. Numbers larger than the type_fh will cause an overflow.
% Modifications are noted in the code below.
%
% type_fh: type function handle such as uint64.m, uint32.m, double.m, etc.
%
% Example:
%  x=bin2dec_type('1010111110000110100110101000110010111110100001001001000000000000', @uint64)
%
% Most of the code is from Matlab's bin2dec function. This function may not
% be used if you do not have a valid Matlab license.
%
%   Copyright 1984-2006 The MathWorks, Inc.

% handle input
if iscellstr(s)
  s = char(s);
end

if ~ischar(s)
  error(message('MATLAB:bin2dec:InvalidInputClass'));
end

if isempty(s)
  x = [];
  return,
end

% MODIFICATION: Check for type argument
if nargin<=1
  type_fh = @double;
end

% MODIFICATION: Remove size check
% if size(s,2)>52
%     error(message('MATLAB:bin2dec:InputOutOfRange'));
% end

% remove significant spaces
for i = 1:size(s,1)
  spacesHere = (s(i,:)==' '|s(i,:)==0);
  if any(spacesHere)
    stmp = s(i,:);                                  % copy this row
    nrOfZeros=sum(spacesHere);                      % number zeros to prepend
    stmp(spacesHere)='';                            % remove significant spaces
    s(i,:) = [repmat('0',1,nrOfZeros) stmp];        % prepend '0' to pad this row
  else
    continue;
  end
end

% check for illegal binary strings
if any(any(~(s == '0' | s == '1')))
  error(message('MATLAB:bin2dec:IllegalBinaryString'))
end

[m,n] = size(s);

% MODIFICATION: Convert to numbers using type_fh
v = type_fh(s - '0');
twos = type_fh(pow2(n-1:-1:0));
x = zeros(m,1,func2str(type_fh));
for s_idx=1:size(s,2)
  for r_idx = 1:size(s,1)
    if s(r_idx,s_idx) == '1'
      x(r_idx) = x(r_idx) + twos(1,s_idx);
    end
  end
end
