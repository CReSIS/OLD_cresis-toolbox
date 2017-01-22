
function push(obj,cmds)
% Pushes a set of commands onto the stack

obj.last_pointer = obj.pointer;
obj.pointer = obj.pointer + 1;
obj.stack{obj.pointer} = cmds;

% Since the stack may contain cmds for indices > obj.pointer (e.g. because
% the user chose to undo some commands and then pushed this command), we
% must remove these commands from the stack.
obj.stack = obj.stack(1:obj.pointer);

% Synch all the documents that are linked to current undo_stack
notify(obj,'synchronize_event');

end