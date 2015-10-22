function cmds_set_undo_stack(obj,undo_stack)
% cmds_set_undo_stack(obj,undo_stack)
%
% Attach or detach undo stack. Also adds or deletes listeners.
%
% Attach:
% undo_stack = undo stack handle class
% Detach:
% undo_stack = empty ([])

if isempty(undo_stack)
  % Detach the current undo stack
  if isa(obj.undo_stack,'imb.undo_stack')
    try
      obj.undo_stack.remove_document(obj);
    end
    try
      delete(obj.undo_stack_synchronize_listener);
    end
    try
      delete(obj.undo_stack_save_listener);
    end
    obj.undo_stack = [];
  end
else
  obj.undo_stack = undo_stack;
  
  obj.undo_stack_save_listener ...
    = addlistener(obj.undo_stack,'save_event',@obj.cmds_save);
  obj.undo_stack_synchronize_listener ...
    = addlistener(obj.undo_stack,'synchronize_event',@obj.cmds_synchronize);
  
  % An undo stack already exists for this system-segment pair, so attach this echowin
  % to it. Since there may be commands in the undo stack already, we will
  % run these commands so that the new echowin is synced up with the stack.
  cmds_list = obj.undo_stack.attach_document(obj);
  obj.cmds_execute(cmds_list,'redo');
end

end
