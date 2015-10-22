function create_ui(obj,img_title)
% column_browser.create_ui(obj,img_title)
%
% Support function for column_browser class. Create user interface for
% column browser figure
%

%--------------------------------------------------------------------------
% Create column browser window
%--------------------------------------------------------------------------
obj.h_fig = figure;
set(obj.h_fig,'DockControls','off')
set(obj.h_fig,'NumberTitle','off');
set(obj.h_fig,'Name',img_title);
set(obj.h_fig,'ToolBar','none');
set(obj.h_fig,'MenuBar','none');
set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
zoom_setup(obj.h_fig);

obj.h_gui.h_axes = axes('Parent',obj.h_fig);
hold(obj.h_gui.h_axes,'on');
grid(obj.h_gui.h_axes,'on');
obj.h_gui.h_plot = [];
obj.h_gui.h_marker = [];
obj.update_plot();

set(obj.h_fig,'WindowButtonUpFcn',@obj.button_up);
set(obj.h_fig,'WindowButtonDownFcn',@obj.button_down);
set(obj.h_fig,'WindowScrollWheelFcn',@obj.button_scroll);
set(obj.h_fig,'WindowKeyPressFcn',@obj.key_press);



return;
