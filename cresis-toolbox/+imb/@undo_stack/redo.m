function redo(obj)
% Repushes a set of commands onto the stack

if obj.pointer < length(obj.stack)
  obj.last_pointer = obj.pointer;
  obj.pointer = obj.pointer + 1;
  
  % Synch all the documents that are linked to current undo_stack
  notify(obj,'synchronize_event');
end

end
