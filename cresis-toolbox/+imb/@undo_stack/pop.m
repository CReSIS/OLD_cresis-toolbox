function pop(obj)
% Pops a set of commands off the stack

if obj.pointer > 0
  obj.last_pointer = obj.pointer;
  obj.pointer = obj.pointer - 1;

  % Synch all the documents that are linked to current undo_stack
  notify(obj,'synchronize_event');
end

end
