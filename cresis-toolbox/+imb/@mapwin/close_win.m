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
    cur_unit = get(obj.h_fig,'Units');
    set(obj.h_fig,'Units','pixels')
    mapwin_pos = get(obj.h_fig,'Position');
    set(obj.h_fig,'Unit',cur_unit);
    obj.default_params.mapwin.x = mapwin_pos(1);
    obj.default_params.mapwin.y = mapwin_pos(2);
    obj.default_params.mapwin.w = mapwin_pos(3);
    obj.default_params.mapwin.h = mapwin_pos(4);
    
    obj.default_params.prefwin = obj.map_pref.default_params;
    
    if length(obj.echowin_list) >= 1
      idx = 1;
      cur_unit = get(obj.echowin_list(idx).h_fig,'Units');
      set(obj.echowin_list(idx).h_fig,'Units','pixels')
      echowin_pos = get(obj.echowin_list(idx).h_fig,'Position');
      set(obj.echowin_list(idx).h_fig,'Unit',cur_unit);
      obj.default_params.echowin.x = echowin_pos(1);
      obj.default_params.echowin.y = echowin_pos(2);
      obj.default_params.echowin.w = echowin_pos(3);
      obj.default_params.echowin.h = echowin_pos(4);
    end
    
    default_params = obj.default_params;
    save(obj.default_params.picker_param_fn,'-append','-struct','default_params','mapwin','prefwin','echowin');
  catch ME
    warning('Failed to save default parameters file %s', obj.default_params.picker_param_fn);
  end
  
  % Delete the object
  delete(obj);
end

return;
