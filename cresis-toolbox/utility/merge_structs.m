function S_out = merge_structs(S_in,S_over)
% S_out = merge_structs(S_in,S_over)
%
% Merges two structures. S_over takes precedence.
%
% S_in = input struct
% S_over = struct which will override any fields in S_in
%
% S_out = merged structure
%
% Example: See bottom of file
%
% Authors: Huan Zhao, John Paden

if isempty(S_over)
  S_out = S_in;
  return;
end
if isempty(S_in)
  S_out = S_over;
  return;
end

S_out = S_in;
outputfields = fieldnames(S_out);
overfields = fieldnames(S_over);

% Traverse the structure as a tree where leafs are copied/overwritten
% and branches/children (type == struct) are dealt with recursively.

for x = 1:length(overfields)
  a = 0 ;
  for y = 1:length(outputfields)
    if strcmp(overfields{x},outputfields{y})
      % If the override has a field that IS present in S_in
      if isstruct(getfield(S_over,overfields{x})) ...
          && length(getfield(S_over,overfields{x})) == 1
        % Recursively call the function for children structs
        S_out.(outputfields{y}) = merge_structs(getfield(S_out, outputfields{y}),getfield(S_over, overfields{x}));
      else
        % Just overwrite the field if:
        % 1. it is not a struct (i.e. a leaf)
        % 2. if it is a struct array where array lengths do not match
        S_out = rmfield(S_out, outputfields{y});
        for idx = 1:length(S_over)
          S_out(idx).(outputfields{y}) = getfield(S_over(idx), overfields{x});
        end
      end
      a=1;
      break;
    end
    
  end
  if a==0
    % If the override has a field that IS NOT present in S_in
    for idx = 1:length(S_over)
      S_out(idx).(overfields{x}) = getfield(S_over(idx), overfields{x});
    end
  end
end

return ;

% ==================================================================
% ==================================================================
% Examples
% ==================================================================
% ==================================================================

S_in.A = 1;
S_in.B = 2;
S_in.C.D = 4;
S_in.C.D2 = 4.2;
S_over.A = 5;
S_over.C.D = 42;
S_over.C.E = 6;
S_over.F = 7;
S_out = merge_structs(S_in,S_over)

