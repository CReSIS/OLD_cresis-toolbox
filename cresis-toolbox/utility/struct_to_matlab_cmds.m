function cmd_strs = struct_to_matlab_cmds(in_struct, out_fn)
% cmd_strs = struct_to_matlab_cmds(in_struct, out_fn)
%
% Creates a set of Matlab commands that will recreate a cell, struct,
% or object (only standard class types are supported).
%
% in_struct: Input structure to be converted to equivalent matlab commands
% out_fn: File name of output file containing the generated matlab commands
%   (optional)
%
% cmd_strs = All the commands in a string (commands terminated with \n).
%   If no output arguments, commands are printed to the screen.
%   Commands will be very long strings for big arrays and may exceed the
%   maximum line length allowed by the Matlab command window.
%
% Author: Sam Buchanan, John Paden

global gRadar;

if ~isstruct(in_struct) && ~isobject(in_struct) && ~iscell(in_struct)
  fprintf('Please supply an input of struct, object, or cell.\n');
  return;
end

if nargin > 1 && ~isempty(nargin(2))
  % Print to a file
  out.type(1) = 1;
  out.fn = out_fn;
  [out.fid,msg] = fopen(out_fn,'w');
  if out.fid < 0
    error('Could not open %s for writing: %s', out_fn, msg);
  end
else
  % Do nothing
  out.type(1) = 0;
  out.fn = '';
end

if nargout < 1
  if out.type(1) == 0;
    % Print to standard out
    out.type(2) = 1;
    else
    % Do nothing
    out.type(2) = 0;
            end
else
  % Print to a string
  out.type(2) = 2;
  out.cmd_strs = '';
end

% Function that actually does the work
struct_to_matlab_cmds_recurse(in_struct,out,inputname(1));

if out.type(1) == 1
  fclose(out.fid);
end

if out.type(2) == 2
  cmd_strs = out.cmd_strs;
end

end


function struct_to_matlab_cmds_recurse(in_struct,out,parent)
% struct_to_matlab_cmds_recurse(in_struct,out,parent)
%
% Hidden function that actually does the work for struct_to_matlab_cmds

if ~iscell(in_struct)
  fields = fieldnames(in_struct);
  num_fields = length(fields);
else
  num_fields = 1;
end

cmd_strs = '';

for struct_idx = 1:numel(in_struct)
  
  for field_idx = 1:num_fields
    if ~iscell(in_struct)
      field = fields{field_idx};
      field_class = class(in_struct(struct_idx).(field));
    else
      field_class = class(in_struct{struct_idx});
  end
    if any(strcmpi(field_class,{'struct','object','cell'}))
      if ~iscell(in_struct)
        struct_to_matlab_cmds_recurse(in_struct(struct_idx).(field),out,sprintf('%s(%d).(''%s'')', parent, struct_idx, field));
      else
        struct_to_matlab_cmds_recurse(in_struct{struct_idx},out,sprintf('%s{%d}', parent, struct_idx));
  end
    elseif any(strcmpi(field_class,{'logical','char','uint8','int8','uint16','int16','uint32','int32','uint64','int64','single','double'}))
      % Eventually should break these commands into multiple lines for
      % large arrays.
      if strcmpi(field_class,'double')
        format_str = '%.17g ';
      elseif strcmpi(field_class,'single')
        format_str = '%.10g ';
      else
        format_str = '%g ';
      end
      if ~iscell(in_struct)
        size_str = strcat('[', sprintf('%d ', size(in_struct(struct_idx).(field))), ']');
        cmd = cat(2, sprintf('%s(%d).(''%s'') = reshape(%s([', parent, struct_idx, field, field_class), ...
          sprintf(format_str, in_struct(struct_idx).(field)), ...
          sprintf(']),%s);\n', size_str));
      else
        size_str = strcat('[', sprintf('%d ', size(in_struct{struct_idx})), ']');
        cmd = cat(2, sprintf('%s{%d} = reshape(%s([', parent, struct_idx, field_class), ...
          sprintf(format_str, in_struct{struct_idx}), ...
          sprintf(']),%s);\n', size_str));
      end
      out = struct_to_matlab_cmds_print(cmd,out);
    elseif strcmpi(field_class,'function_handle')
      if ~iscell(in_struct)
        cmd = sprintf('%s(%d).(''%s'') = @%s;\n', ...
          parent, struct_idx, field, func2str(in_struct(struct_idx).(field)));
      end
      out = struct_to_matlab_cmds_print(cmd,out);
      else
      warning('No support for %s class', field_class);
      end
    end
  
end

end

function out = struct_to_matlab_cmds_print(str,out)
% Hidden function for struct_to_matlab_cmds that outputs the commands
% in the desired format

if out.type(1) == 1
  fwrite(out.fid,str,'char');
end
if out.type(2) == 1
  fprintf('%s',str);
elseif out.type(2) == 2
  cmd_strs = [cmd_strs fprintf('%s',str)];
end

end
