function added = struct_add_field(current,new_field,value,force_overwrite,no_size_change)
% added = struct_add_field(current,new_field,value,force_overwrite,no_size_change)
%
% Function for adding a field to a structure without changing the size of
% the structure (even for empty structures) or putting anything into the
% new field unless required. Structure which contain structures must have a
% size of 1 though. A value or values can be specified too.
%
% current: current structure before field is added
% new_field: string containing the field name. Can also include periods
%   which will cause it to create a substructure (see examples).
% value: cell array of values to assign to the field
% force_overwrite: do not write a value into a field that already exists
%   unless this is set to true (default is false)
% no_size_change: value field is ignored... this mode tries to only make
%   empty structs (default is false)
%
% added: structure with the new field added
%
% S = struct_add_field([],'A');
% S = struct_add_field(S,'B');
% S = struct_add_field(S,'B.C.D');
%
% S = struct_add_field([],'B.C.D');
% S = struct_add_field(S,'B.C.E');
% S = struct_add_field(S,'F');
%
% S = struct_add_field([],'B.C.D',{[]});
% S = struct_add_field(S,'B.C.E',{[]});
% S = struct_add_field(S,'F',{[]});
%
% S = struct_add_field([],'B.C.D',{1,2});
% S(1) = struct_add_field(S(1),'B.C.E',{3,4});
% S(2) = struct_add_field([],'B.C.D',{5,6});
% S(2) = struct_add_field(S(2),'B.C.E',{7,8});
% S = struct_add_field(S,'F',{9,10});
%
% Author: John Paden

if ~exist('force_overwrite','var')
  force_overwrite = false;
end
if ~exist('no_size_change','var')
  no_size_change = false;
end
if ~exist('value','var')
  value = [];
end

period_idx = find(new_field == '.',1);
end_name_idx = find(new_field == '.' | new_field == '{',1);
if isempty(end_name_idx)
  new_field_name = new_field;
else
  new_field_name = new_field(1:end_name_idx-1);
end

added = current;
if ~isempty(period_idx)
  if force_overwrite  || isempty(added) || ~isfield(added,new_field_name)
    added(1).(new_field_name) = struct([]);
  end
  added(1).(new_field_name) = struct_add_field(added(1).(new_field_name), ...
    new_field(period_idx+1:end),value,force_overwrite,no_size_change);
else
  if no_size_change
    % Trick to add a new field to the struct and keep the struct empty
    if ~isstruct(added) || isempty(fieldnames(added))
      added = struct(new_field_name,{});
    elseif force_overwrite || ~isfield(added,new_field_name)
      [added.(new_field_name)] = deal([]);
    end
  else
    if force_overwrite || ~isfield(added,new_field_name)
      if isempty(added)
        added = struct(new_field_name,value);
      else
        [added.(new_field_name)] = deal(value{:});
      end
    end
  end
end

end
