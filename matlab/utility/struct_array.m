function M = struct_array(struct,struct_name,fill_type,fill_val)
% Returns the maximum dimension array required to fit every element.
% N-dim struct arrays are flattened into 1xN
% Array will be 1:1 match in size only when every struct entry is the
% same size. 
% Array type will be the last type encountered.
%
% % Example of struct that will not fill the array
% clear s1;
% s1.s2(1).s3(1).f1 = uint8([2 3 4]);
% s1.s2(2).s3(1).f1 = double([6 7 8 9]);
% s1.s2(2).s3(2).f1 = [4 5];
% 
% M = struct_array(s1,'s2.s3.f1','double',NaN)
% 
% % Example of struct that will fill the array
% clear s1;
% s1.s2(1).s3(1).f1 = [2 3 4 2];
% s1.s2(1).s3(2).f1 = [4 5 6 7];
% s1.s2(2).s3(1).f1 = [8 9 10 11];
% s1.s2(2).s3(2).f1 = [12 13 14 15];
% 
% M = struct_array(s1,'s2.s3.f1','double',NaN)
%
% Author: John Paden
%
% See also: struct_array.m

struct_size = struct_array_size(struct,struct_name);

C = textscan(struct_name,'%s','Delimiter','.');
field_names = C{1};

state.fill_type = fill_type;
state.fill_type = fill_val;
state.fill_type = fill_type;
M = struct_array_recurse(struct,field_names,fill_type,fill_val,struct_size);

return

function M = struct_array_recurse(struct,field_names,fill_type,fill_val,struct_size)

if length(field_names) > 1
  M = fill_val * ones(struct_size, fill_type);
  for idx = 1:numel(struct.(field_names{1}))
    val = struct_array_recurse(struct.(field_names{1})(idx),field_names(2:end),fill_type,fill_val,struct_size(2:end));
    val_size = size(val);
    S.type = '()';
    S.subs = {idx};
    for val_dim = 1:length(val_size)
      S.subs{val_dim+1} = 1:size(val,val_dim);
    end
    M = subsasgn(M,S,val);
  end
else
  M = feval(fill_type,struct.(field_names{1}));
end

return

