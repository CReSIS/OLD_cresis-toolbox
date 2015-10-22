function cmds_list = get_save_cmds(obj,remove_cmds)
% cmds_list = get_save_cmds(obj,remove_cmds)
%
% This command returns all the commands that need to be saved.  Removes
% the saved commands if remove_commands set to true.
%
% cmds_list = cell vector of commands
% cmds_direction is assumed to be 'redo'

cmds_list = obj.stack(1:obj.pointer);

if remove_cmds
  obj.stack = obj.stack(obj.pointer+1:end);
  if isempty(obj.stack)
    obj.pointer = 0;
  else
    obj.pointer = 1;
  end
  obj.last_pointer = NaN;
end

end
