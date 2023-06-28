function out = print_struct(in,mode)
% out = print_struct(in,mode)
%
% Create a string for output to the matlab command window that shows the
% contents on a structure. Only numeric and char fields are supported.
%
% in: structure
% mode: bit mask representing:
%   0x0001: 0: uses spaces, 1: uses tabs to separate fields
%   0x0002: 0: regular structure array, 1: flat structure where each field
%                                          is the same length
% out: string which can be printed via fprintf(out) for display
%
% Examples:
% in = struct('field1',2,'field2',3,'field3','asdf')
% in = struct('field1',{2 3},'field2',{3 123123},'field3',{'asdf', 'Longer field than first'})
%
% Author: John Paden

if ~exist('mode','var')
  mode = 0;
end

fields = fieldnames(in);
num_fields = numel(fields);

out = '';
field_strs = {};
for field_idx = 1:num_fields
  field_strs{1,field_idx} = fields{field_idx};
end
if bitand(mode,2)
  % Single struct element with arrays in each field
  num_elements = length(in(1).(fields{1}));
else
  % Typical struct array
  num_elements = length(in);
end
for struct_idx = 1:num_elements
  for field_idx = 1:num_fields
    if bitand(mode,2)
      field = in(1).(fields{field_idx})(struct_idx);
    else
      field = in(struct_idx).(fields{field_idx});
    end
    if ischar(field)
      field(field==10|field==13) = 0;
      field_strs{struct_idx+1,field_idx} = field;
    elseif isnumeric(field)
      if isinteger(field)
        field_strs{struct_idx+1,field_idx} = sprintf('%d',field);
      else
        field_strs{struct_idx+1,field_idx} = sprintf('%g',field);
      end
    end
  end
end
field_lens = max(cellfun(@length,field_strs));

row = 1;
for col = 1:size(field_strs,2)
  if bitand(mode,1)
    format_str = sprintf('%%s\t', field_lens(col));
  else
    format_str = sprintf('<strong>%%%ds</strong> ', field_lens(col));
  end
  new_out = sprintf(format_str,field_strs{row,col});
  out(end+(1:numel(new_out))) = new_out;
end
out(end+1) = 13;
for row = 2:size(field_strs,1)
  for col = 1:size(field_strs,2)
  if bitand(mode,1)
    format_str = sprintf('%%s\t', field_lens(col));
  else
    format_str = sprintf('%%%ds ', field_lens(col));
  end
    new_out = sprintf(format_str,field_strs{row,col});
    out(end+(1:numel(new_out))) = new_out;
  end
  out(end+1) = 13;
end
