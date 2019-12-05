function close_win(obj,varargin)
% Go throught undo_stack list and check if there is unsaved stacks
saved = true;
for idx = 1:length(obj.undo_stack_list)
  if obj.undo_stack_list(idx).ismodified()
    saved = false;
    break;
  end
end
if saved
  % close window
  close = true;
else
  prompt = questdlg('There are unsaved layers, do you want to save them?',...
    'Unsaved','Yes','No','Cancel','Yes');
  switch prompt
    case 'Yes'
      % save first, then close window
      for idx = 1:length(obj.undo_stack_list)
        if ~isempty(obj.undo_stack_list(idx).stack)
          imb.save_undo_stack(obj.undo_stack_list(idx));
        end
      end
      close = true;
    case 'No'
      % Do not save, close window
      close = true;
    case 'Cancel'
      % Go back
      close = false;
    case ''
      % Go back
      close = false;
  end
end

if close
  try
    obj.save_default_params();
  catch ME
    warning('Failed to save default parameters file %s', obj.default_params.picker_param_fn);
  end
  
  % Delete the object
  delete(obj);
end
