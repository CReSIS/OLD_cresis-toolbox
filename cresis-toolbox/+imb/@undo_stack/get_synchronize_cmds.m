function [cmds_list,cmds_direction] = get_synchronize_cmds(obj)
% [cmds_list,cmds_direction] = get_synchronize_cmds(obj)
%
% This command should be called by each document when a synchronize event
% occurs so that each document can execute the cmds.
%
% cmds_list = cell vector of commands
% cmds_direction = 'undo' or 'redo'

if obj.pointer > obj.last_pointer
  cmds_direction = 'redo';
  cmds_list = obj.stack(obj.last_pointer+1:obj.pointer);
else
  cmds_direction = 'undo';
  cmds_list = obj.stack(obj.last_pointer:-1:obj.pointer+1);
end

end
