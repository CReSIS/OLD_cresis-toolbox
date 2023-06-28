function result = compare_structs(struct1,struct2,print_differences,prestring)
% result = compare_structs(struct1,struct2,print_differences,prestring)
%
% Returns 1 if structs differ, 0 otherwise
%
% Cell arrays and inline functions cannot be compared
%
% struct1: First struct to compare
% struct2: Second struct to compare
% print_differences: double, 0 does not print, 1 prints differences, 2 prints
%   differences and fields which could not be compared
% prestring: DO NOT USE (used for internal recursion)
%
% Author: John Paden

if ~exist('print_differences','var')
  print_differences = 0;
end
if ~exist('prestring','var')
  prestring = '';
end

result = 0;

fn1 = fieldnames(struct1);
fn2 = fieldnames(struct2);
fn2_mask = zeros(size(fn2));

for idx = 1:length(fn1)
  match_idx = strmatch(fn1{idx},fn2,'exact');
  if isempty(match_idx)
    if print_differences
      fprintf('Struct 2 does not have field %s%s\n', prestring, fn1{idx})
    end
    result = 1;
  else
    % Both structs have this field
    fn2_mask(match_idx) = 1;
    if iscell(struct1.(fn1{idx}))
      if print_differences >= 2
        fprintf('Field %s%s is cell array, skipping\n', prestring, fn1{idx});
      end
      continue;
    end
    if isa(struct1.(fn1{idx}),'function_handle')
      if isa(struct2.(fn1{idx}),'function_handle')
        if strcmp(func2str(struct1.(fn1{idx})),func2str(struct2.(fn1{idx})))
          continue;
        else
          if print_differences
            fprintf('Field %s%s does not match\n', prestring, fn1{idx});
          end
        end
      else
        if print_differences
          fprintf('Field %s%s is function handle for struct1 and not struct2\n', prestring, fn1{idx});
        end
      end
      continue;
    end
    if isa(struct1.(fn1{idx}),'inline')
      if print_differences >= 2
        fprintf('Field %s%s is an inline function handle, skipping\n', prestring, fn1{idx});
      end
      continue;
    end
    % This field is not a cell
    if isstruct(struct1.(fn1{idx}))
      if isstruct(struct2.(fn1{idx}))
        if length(struct1.(fn1{idx})) ~= length(struct2.(fn1{idx}))
          if print_differences
            fprintf('Field %s%s is struct with different length\n', prestring, fn1{idx});
          end
          result = 1;
        end
        for struct_idx = 1:min(length(struct1.(fn1{idx})),length(struct2.(fn1{idx})))
          new_result = compare_structs(struct1.(fn1{idx})(struct_idx),struct2.(fn1{idx})(struct_idx),print_differences,[prestring fn1{idx} '.']);
          if new_result
            result = 1;
          end
        end
      end
    elseif isstruct(struct2.(fn1{idx}))
      if print_differences
        fprintf('Field %s%s is struct for struct2, but not struct1\n', prestring, fn1{idx});
      end
      result = 1;
    else
      s1 = struct1.(fn1{idx})(:);
      s2 = struct2.(fn1{idx})(:);
      if ndims(struct1.(fn1{idx})) ~= ndims(struct2.(fn1{idx})) ...
          || any(size(struct1.(fn1{idx})) ~= size(struct2.(fn1{idx}))) ...
          || any(isnan(s1) ~= isnan(s2)) ...
          || any(s1(~isnan(s1)) ~= s2(~isnan(s1)))
        if print_differences
          fprintf('Field %s%s does not match\n', prestring, fn1{idx});
          if print_differences >= 3
            struct1_val = struct1.(fn1{idx})
            struct2_val = struct2.(fn2{match_idx})
          end
        end
        result = 1;
      end
    end
  end
end

for idx = 1:length(fn2)
  if fn2_mask(idx) == 0
    if print_differences
      fprintf('Struct 1 does not have field %s%s\n', prestring, fn2{idx})
    end
    result = 1;
  end
end

return
