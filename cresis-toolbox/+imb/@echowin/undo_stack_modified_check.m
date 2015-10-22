function cancel_operation = undo_stack_modified_check(obj)

% Check to see if undo stack got created and, if it did, if there are
% any items in the undostack (i.e. operations that have not been saved)
if ~obj.undo_stack.ismodified()
  cancel_operation = false;
else
  prompt = questdlg('There are unsaved layers. You do not need to save them now, but would you like to?',...
    'Unsaved Layers','Yes','No','Cancel','Yes');
  switch prompt
    case 'Yes'
      % save
      obj.savePB_callback();
      cancel_operation = false;
    case 'No'
      cancel_operation = false;
    case 'Cancel'
      % Go back
      cancel_operation = true;
    case ''
      % Go back
      cancel_operation = true;
  end
end

return;
