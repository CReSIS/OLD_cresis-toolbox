function cmds_synchronize(obj,varargin)
% cmds_synchronize(obj,varargin)
%
% Listener to undo stack synchronize_event

% Get the list of commands from the undo stack that need to be run to
% synchronize the echowin
[cmds_list,cmds_direction] = obj.undo_stack.get_synchronize_cmds();

% Execute the commands
obj.cmds_execute(cmds_list,cmds_direction);

% Check to see if any modifications have been done relative to the last
% save.
if obj.undo_stack.ismodified()
  set(obj.left_panel.savePB,'String','(S)ave Layer*');
else
  set(obj.left_panel.savePB,'String','(S)ave Layer');
end

return;
