function [struct_size,type] = struct_array_size(struct,struct_name)
% Returns the maximum dimension array required to fit every element.
% N-dim struct arrays are flattened into 1xN
% Array will be 1:1 match in size only when every struct entry is the
% same size. 
% Array type will be the last type encountered.
%
% % Example of struct that will not fill the array
% s1.s2(1).s3(1).f1 = uint8([2 3 4]);
% s1.s2(2).s3(1).f1 = double([6 7 8 9]);
% s1.s2(2).s3(2).f1 = [4 5];
% 
% struct_array_size(s1,'s2.s3.f1')
% 
% % Example of struct that will fill the array
% s1.s2(1).s3(1).f1 = [2 3 4 2];
% s1.s2(1).s3(2).f1 = [4 5 6 7];
% s1.s2(2).s3(1).f1 = [8 9 10 11];
% s1.s2(2).s3(2).f1 = [12 13 14 15];
% 
% struct_array_size(s1,'s2.s3.f1')
%
% Author: John Paden
%
% See also: struct_array.m

C = textscan(struct_name,'%s','Delimiter','.');
field_names = C{1};
[struct_size,type] = struct_array_size_recurse(struct,field_names);

return

function [struct_size,type] = struct_array_size_recurse(struct,field_names)

if length(field_names) > 1
  max_size = [];
  for idx = 1:numel(struct.(field_names{1}))
    [struct_size,type] = struct_array_size_recurse(struct.(field_names{1})(idx),field_names(2:end));
    if isempty(max_size)
      max_size = struct_size;
    else
      max_size = max(max_size,struct_size);
    end
  end
  struct_size = [numel(struct.(field_names{1})) max_size];
else
  type = class(struct.(field_names{1}));
  struct_size = size(struct.(field_names{1}));
end

return



