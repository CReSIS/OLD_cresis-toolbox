function close_echowin(obj,h_obj,event)

%% Find this echowindow in the echowin list
idx = find(obj.echowin_list == h_obj);

%% Save current parameters as defaults for next time an echowin is opened
try
  obj.save_default_params();
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
