function close_echowin(obj,h_obj,event)

%% Find this echowindow in the echowin list
idx = find(obj.echowin_list == h_obj);

%% Save current parameters as defaults for next time an echowin is opened
try
  cur_unit = get(obj.echowin_list(idx).h_fig,'Units');
  set(obj.echowin_list(idx).h_fig,'Units','pixels')
  echowin_pos = get(obj.echowin_list(idx).h_fig,'Position');
  set(obj.echowin_list(idx).h_fig,'Unit',cur_unit);
  obj.default_params.echowin.x = echowin_pos(1);
  obj.default_params.echowin.y = echowin_pos(2);
  obj.default_params.echowin.w = echowin_pos(3);
  obj.default_params.echowin.h = echowin_pos(4);
  default_params = obj.default_params;
  save(obj.default_params.picker_param_fn,'-append','-struct','default_params','echowin');
catch ME
  warning('Failed to save default parameters file %s', obj.default_params.picker_param_fn);
end

%% Delete Echowindow
try
  delete(obj.echowin_list(idx));
end
obj.echowin_list = obj.echowin_list([1:idx-1 idx+1:end]);

%% Update echo window popup menu
picker_window_selection = get(obj.top_panel.picker_windowPM,'Value');
if picker_window_selection > idx
  set(obj.top_panel.picker_windowPM,'Value',picker_window_selection - 1);
end
menu_string = get(obj.top_panel.picker_windowPM,'String');
menu_string = menu_string([1:idx idx+2:end]);
set(obj.top_panel.picker_windowPM,'String',menu_string);

%% Delete map objects associated with this echo window
try
  delete(obj.echowin_maps(idx).h_line);
end
try
  delete(obj.echowin_maps(idx).h_cursor);
end
try
  delete(obj.echowin_maps(idx).h_text);
end
obj.echowin_maps = obj.echowin_maps([1:idx-1 idx+1:end]);

return;
