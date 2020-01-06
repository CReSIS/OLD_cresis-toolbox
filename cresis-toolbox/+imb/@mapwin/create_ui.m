function create_ui(obj)
% create_ui(obj)
%
% Create user interface for mapwin class. Called only from the constructor.

%==========================================================================
%% Create map figure
%==========================================================================
figure(obj.h_fig);
set(obj.h_fig,'Position',[obj.default_params.mapwin.x obj.default_params.mapwin.y obj.default_params.mapwin.w obj.default_params.mapwin.h]);
set(obj.h_fig,'DockControls','off')
set(obj.h_fig,'NumberTitle','off');
if strcmpi(class(obj.h_fig),'double')
  set(obj.h_fig,'Name',sprintf('%d: map',obj.h_fig));
else
  set(obj.h_fig,'Name',sprintf('%d: map',obj.h_fig.Number));
end
set(obj.h_fig,'ToolBar','none');
set(obj.h_fig,'MenuBar','none'); 
set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
set(obj.h_fig,'Interruptible','off');
set(obj.h_fig,'Pointer','custom');
zoom_pointer = [NaN   NaN   NaN   NaN     1     1     1     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
  NaN   NaN     1     1   NaN     2   NaN     2     1     1   NaN   NaN   NaN   NaN   NaN   NaN
  NaN     1     2   NaN     2     1     1   NaN     2   NaN     1   NaN   NaN   NaN   NaN   NaN
  NaN     1   NaN     2   NaN     1     1     2   NaN     2     1   NaN   NaN   NaN   NaN   NaN
  1   NaN     2   NaN     2     1     1   NaN     2   NaN     2     1   NaN   NaN   NaN   NaN
  1     2     1     1     1     1     1     1     1     1   NaN     1   NaN   NaN   NaN   NaN
  1   NaN     1     1     1     1     1     1     1     1     2     1   NaN   NaN   NaN   NaN
  1     2   NaN     2   NaN     1     1     2   NaN     2   NaN     1   NaN   NaN   NaN   NaN
  NaN     1     2   NaN     2     1     1   NaN     2   NaN     1   NaN   NaN   NaN   NaN   NaN
  NaN     1   NaN     2   NaN     1     1     2   NaN     2     1     2   NaN   NaN   NaN   NaN
  NaN   NaN     1     1     2   NaN     2   NaN     1     1     1     1     2   NaN   NaN   NaN
  NaN   NaN   NaN   NaN     1     1     1     1   NaN     2     1     1     1     2   NaN   NaN
  NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     1     1     2   NaN
  NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     1     1     2
  NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     1     1
  NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     2];
search_icon = [
  0   0   0   0     1     1     1     1   0   0   0   0   0   0   0   0
  0   0   1   1     0     0     0     1   1   1   0   0   0   0   0   0
  0   1   1   0     0     0     0     0   0   0   1   0   0   0   0   0
  0   1   0   0     0     0     0     0   0   2   1   0   0   0   0   0
  1   0   0   0     0     0     0     0   0   0   2   1   0   0   0   0
  1   2   0   0     0     0     0     0   0   0   0   1   0   0   0   0
  1   0   0   0     0     0     0     0   0   0   2   1   0   0   0   0
  0   1   2   0     0     0     0     0   0   0   1   0   0   0   0   0
  0   1   0   0     0     0     0     0   0   2   1   0   0   0   0   0
  0   1   0   0     0     0     0     0   0   2   1   2   0   0   0   0
  0   0   1   1     2     0     2     0   1   1   1   1   2   0   0   0
  0   0   0   0     1     1     1     1   0     2     1     1     1     2   0   0
  0   0   0   0   0   0   0   0   0   0     2     1     1     1     2   0
  0   0   0   0   0   0   0   0   0   0   0     2     1     1     1     2
  0   0   0   0   0   0   0   0   0   0   0   0     2     1     1     1
  0   0   0   0   0   0   0   0   0   0   0   0   0     2     1     2];
search_icon(search_icon==0) = NaN;
search_icon(search_icon==1) = 0;
search_icon(search_icon==2) = 0.5;
search_icon = repmat(double(search_icon),[1 1 3]);
set(obj.h_fig,'PointerShapeCData',zoom_pointer);
set(obj.h_fig,'PointerShapeHotSpot',[6 6])

%==========================================================================
%% Create main table content + map axes
%==========================================================================

% ----top_panel
obj.top_panel.handle = uipanel('Parent',obj.h_fig);

%---- map panel
obj.map_panel.handle = uipanel('Parent',obj.h_fig);
obj.map_panel.h_axes = axes('Parent',obj.map_panel.handle,'Visible','off');
colormap(obj.map_panel.h_axes,gray(256));
hold(obj.map_panel.h_axes,'on');
obj.map_panel.h_image = imagesc(0,'Parent',obj.map_panel.h_axes,'Visible','off');
xlabel(obj.map_panel.h_axes,'X (km)');
ylabel(obj.map_panel.h_axes,'Y (km)');
obj.map_panel.h_flightline = plot(1,1,'Parent',obj.map_panel.h_axes,'Color','blue','LineWidth',1);
set(obj.map_panel.h_flightline,'XData',[],'YData',[]);
obj.map_panel.h_cur_sel = plot(1,1,'Parent',obj.map_panel.h_axes,'Color','red','LineWidth',4);
set(obj.map_panel.h_cur_sel,'XData',[],'YData',[]);
set(obj.map_panel.h_axes,'Visible','off');

%----status panel (status bar)
obj.status_panel.handle = uipanel('Parent',obj.h_fig);

%==========================================================================
%% Setup main table
%==========================================================================
obj.table.ui=obj.h_fig;
obj.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.table.height_margin = NaN*zeros(30,30);
obj.table.false_width = NaN*zeros(30,30);
obj.table.false_height = NaN*zeros(30,30);
obj.table.offset = [0 0];

row = 1; col = 1;
obj.table.handles{row,col}        = obj.top_panel.handle;
obj.table.width(row,col)          = inf;
obj.table.height(row,col)         = 25;
obj.table.width_margin(row,col)   = 0;
obj.table.height_margin(row,col)  = 0;

row = 2; col = 1;
obj.table.handles{row,col}        = obj.map_panel.handle;
obj.table.width(row,col)          = inf;
obj.table.height(row,col)         = inf;
obj.table.width_margin(row,col)   = 0;
obj.table.height_margin(row,col)  = 0;

row = 3; col = 1;
obj.table.handles{row,col}        = obj.status_panel.handle;
obj.table.width(row,col)          = inf;
obj.table.height(row,col)         = 12;
obj.table.width_margin(row,col)   = 0;
obj.table.height_margin(row,col)  = 0;

clear row col
table_draw(obj.table);


%==========================================================================
%% Create top_panel table contents
%==========================================================================

%----search input text box
obj.top_panel.searchTB = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.searchTB,'Style','edit');
set(obj.top_panel.searchTB,'String','');
if isunix
  set(obj.top_panel.searchTB,'FontSize',9);
else
  set(obj.top_panel.searchTB,'FontSize',8);
end
set(obj.top_panel.searchTB,'Callback',@obj.search_callback);
set(obj.top_panel.searchTB,'TooltipString','Enter frame ID to search for here');

%----search push buttion
obj.top_panel.searchPB = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.searchPB,'Style','PushButton');
set(obj.top_panel.searchPB,'String','');
set(obj.top_panel.searchPB,'CData',search_icon);
set(obj.top_panel.searchPB,'Callback',@obj.search_callback);
set(obj.top_panel.searchPB,'TooltipString','Search database for frame ID');

%----GIS push button
obj.top_panel.gisPB = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.gisPB,'Style','PushButton');
set(obj.top_panel.gisPB,'String','GIS');
set(obj.top_panel.gisPB,'Callback',@obj.gisPB_callback);
set(obj.top_panel.gisPB,'TooltipString','Open GIS window to load raster and vector data');

%----preference push button
obj.top_panel.preferencePB = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.preferencePB,'Style','PushButton');
set(obj.top_panel.preferencePB,'String','Prefs');
set(obj.top_panel.preferencePB,'Callback',@obj.prefPB_callback);
set(obj.top_panel.preferencePB,'TooltipString','Open preference window');

%----Track checkbox
obj.top_panel.trackCB = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.trackCB,'Style','CheckBox');
set(obj.top_panel.trackCB,'String','Track Echogram');
set(obj.top_panel.trackCB,'Value',true);
set(obj.top_panel.trackCB,'TooltipString','Map tracks movement in echogram');

%---- flight label 
obj.top_panel.flightLabel = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.flightLabel,'Style','text');
set(obj.top_panel.flightLabel,'String','');
set(obj.top_panel.flightLabel,'TooltipString','Selected frame ID');
obj.top_panel.flightCM = uicontextmenu('Parent',obj.h_fig);
% Define the context menu items and install their callbacks
uimenu(obj.top_panel.flightCM, 'Label', 'Copy', 'Callback', @obj.flightCM_callback);
set(obj.top_panel.flightLabel,'UIContextMenu',obj.top_panel.flightCM);

%----picker window pop up menu
obj.top_panel.picker_windowPM = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.picker_windowPM,'String',{'New Window'});
set(obj.top_panel.picker_windowPM,'Value',1);
set(obj.top_panel.picker_windowPM,'Style','popupmenu');
set(obj.top_panel.picker_windowPM,'HorizontalAlignment','Center');
set(obj.top_panel.picker_windowPM,'FontName','fixed');
set(obj.top_panel.picker_windowPM,'TooltipString','Specify window to load echogram into');

%----load push button
obj.top_panel.loadPB = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.loadPB,'Style','PushButton');
set(obj.top_panel.loadPB,'String','Load'); 
set(obj.top_panel.loadPB,'Callback',@obj.loadPB_callback);
set(obj.top_panel.loadPB,'TooltipString','Load selected echogram'); 


%==========================================================================
%% Set up top_panel table
%==========================================================================
obj.top_panel.table.ui = obj.top_panel.handle;
obj.top_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.top_panel.table.height_margin = NaN*zeros(30,30);
obj.top_panel.table.false_width = NaN*zeros(30,30);
obj.top_panel.table.false_height = NaN*zeros(30,30);
obj.top_panel.table.offset = [0 0];

row = 1; col = 1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.searchTB;
obj.top_panel.table.width(row,col)     = 74;
obj.top_panel.table.height(row,col)    = 20;
obj.top_panel.table.width_margin(row,col) = 0;
obj.top_panel.table.height_margin(row,col) = 0;

row = 1; col = col+1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.searchPB;
obj.top_panel.table.width(row,col)     = 28;
obj.top_panel.table.height(row,col)    = 20;
obj.top_panel.table.width_margin(row,col) = 0;
obj.top_panel.table.height_margin(row,col) = 0;

row = 1; col = col+1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.gisPB;
obj.top_panel.table.width(row,col)     = 28;
obj.top_panel.table.height(row,col)    = 20;
obj.top_panel.table.width_margin(row,col) = 3;
obj.top_panel.table.height_margin(row,col) = 0;

row = 1; col = col+1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.preferencePB;
obj.top_panel.table.width(row,col)     = 33;
obj.top_panel.table.height(row,col)    = 20;
obj.top_panel.table.width_margin(row,col) = 3;
obj.top_panel.table.height_margin(row,col) = 0;

row = 1; col = col+1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.trackCB;
obj.top_panel.table.width(row,col)     = 100;
obj.top_panel.table.height(row,col)    = 20;
obj.top_panel.table.width_margin(row,col) = 3;
obj.top_panel.table.height_margin(row,col) = 0;

row = 1; col = col+1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.flightLabel;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = 20;
obj.top_panel.table.width_margin(row,col) = 0;
obj.top_panel.table.height_margin(row,col) = 0;

row = 1; col = col+1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.picker_windowPM;
obj.top_panel.table.width(row,col)     = 100;
obj.top_panel.table.height(row,col)    = 20;
obj.top_panel.table.width_margin(row,col) = 3;
obj.top_panel.table.height_margin(row,col) = 0;

row = 1; col = col+1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.loadPB;
obj.top_panel.table.width(row,col)     = 45;
obj.top_panel.table.height(row,col)    = 20;
obj.top_panel.table.width_margin(row,col) = 3;
obj.top_panel.table.height_margin(row,col) = 0;

% Draw table
table_draw(obj.top_panel.table);

%==========================================================================
%% Create status_panel components
%==========================================================================

% create statusbar's elements:mouse coordinate info
% mouse coordinate info display
obj.status_panel.mouseCoordText = uicontrol('parent',obj.status_panel.handle);
set(obj.status_panel.mouseCoordText,'Style','text');
%set(obj.status_panel.mouseCoordText,'FontName','courier');
if isunix
  set(obj.status_panel.mouseCoordText,'FontSize',9)
else
  set(obj.status_panel.mouseCoordText,'FontSize',8)
end
set(obj.status_panel.mouseCoordText,'HorizontalAlignment','right');
set(obj.status_panel.mouseCoordText,'String','');

% %----mapwin context menu
% obj.status_panel.mapCM= uicontextmenu('parent',obj.h_fig);
% % Define the context menu items and install their callbacks
% obj.status_panel.mapCM_item2 = uimenu(obj.status_panel.mapCM, 'Label', 'Crossovers');%, 'Callback', @obj.sourceLB_callback);
% % attach context menu to axes
% set(obj.map_panel.h_axes,'uicontextmenu',obj.status_panel.mapCM);
% set(obj.status_panel.mouseCoordText,'uicontextmenu',obj.status_panel.mapCM);

%==========================================================================
%% Create status_panel's table
%==========================================================================
obj.status_panel.table.ui = obj.status_panel.handle;
obj.status_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.status_panel.table.height_margin = NaN*zeros(30,30);
obj.status_panel.table.false_width = NaN*zeros(30,30);
obj.status_panel.table.false_height = NaN*zeros(30,30);
obj.status_panel.table.offset = [0 0];

obj.status_panel.table.ui = obj.status_panel.handle;
obj.status_panel.table.offset = [0 2];
row = 1; col = 1;
obj.status_panel.table.handles{row,col}       = obj.status_panel.mouseCoordText;
obj.status_panel.table.width(row,col)         = inf;
obj.status_panel.table.height(row,col)        = inf;
obj.status_panel.table.width_margin(row,col)  = 0;
obj.status_panel.table.height_margin(row,col) = 0;
clear row col
table_draw(obj.status_panel.table);

clear row col

%==========================================================================
%% Connect map figure callbacks
%==========================================================================
set(obj.h_fig,'WindowButtonUpFcn',@obj.button_up)
set(obj.h_fig,'WindowButtonDownFcn',@obj.button_down);
set(obj.h_fig,'WindowScrollWheelFcn',@obj.button_scroll);
set(obj.h_fig,'WindowButtonMotionFcn',@obj.button_motion);
set(obj.h_fig,'WindowKeyPressFcn',@obj.key_press);
set(obj.h_fig,'WindowKeyReleaseFcn',@obj.key_release);
