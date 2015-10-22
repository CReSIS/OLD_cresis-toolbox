function open_crossover_echowin(obj,source,event)
% mapwin.open_crossover_echowin(obj,source,event)

% Set the currently selected frame 
obj.cur_sel = source.get_crossover();

% Force to open new window
set(obj.top_panel.picker_windowPM,'Value',1);

% Load the currently selected frame
loadPB_callback(obj);

end
