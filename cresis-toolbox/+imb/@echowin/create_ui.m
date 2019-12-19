function create_ui(obj)
% echowin.create_ui(obj)
%
% Support function for imb.echowin class. Create user interface for
% echogram window
%

%--------------------------------------------------------------------------
% Create pick window
%--------------------------------------------------------------------------
set(obj.h_fig,'Position',[obj.default_params.x obj.default_params.y obj.default_params.w obj.default_params.h]);
set(obj.h_fig,'DockControls','off')
set(obj.h_fig,'NumberTitle','off');
if strcmpi(class(obj.h_fig),'double')
  set(obj.h_fig,'Name',sprintf('%d: pick',obj.h_fig));
else
  set(obj.h_fig,'Name',sprintf('%d: pick',obj.h_fig.Number));
end
set(obj.h_fig,'ToolBar','none');
set(obj.h_fig,'MenuBar','none');
set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
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
set(obj.h_fig,'PointerShapeCData',zoom_pointer);
set(obj.h_fig,'PointerShapeHotSpot',[6 6])

%==============================================================
% figure (pick) table content
% ----left_panel
obj.left_panel.handle = uipanel('Parent',obj.h_fig);
%----right_panel
obj.right_panel.handle = uipanel('Parent',obj.h_fig);
%---------------------------------------------------------------
% echowin table setup
obj.table.ui=obj.h_fig;
obj.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.table.height_margin = NaN*zeros(30,30);
obj.table.false_width = NaN*zeros(30,30);
obj.table.false_height = NaN*zeros(30,30);
obj.table.offset = [0 0];
row = 1; col = 1;
obj.table.handles{row,col}   = obj.left_panel.handle;
obj.table.width(row,col)     = 180;
obj.table.height(row,col)    = inf;
obj.table.width_margin(row,col)=0;
obj.table.height_margin(row,col)=0;
row = 1; col = 2;
obj.table.handles{row,col}   = obj.right_panel.handle;
obj.table.width(row,col)     = inf;
obj.table.height(row,col)    = inf;
obj.table.width_margin(row,col)=0;
obj.table.height_margin(row,col)=0;
clear row col
table_draw(obj.table);

%================================================================
%================================================================
%% Left_panel table contents
%================================================================
%================================================================

obj.tool_list = {};

obj.tool_list{end+1} = imb.picktool_interpolate;
addlistener(obj.tool_list{end},'hide_param',@obj.toolparam_close_callback);
% Left click: Enters a point based on parameters
%   Find max in range (specify range line/bin extent to search)
%   Find leading edge in range
%   Recomputes interp if point is within last interpolation range
% Left click and drag: Interpolates between manual points based on paramaters
%   Deletes all previous automated points in range
%   Interpolation tools: linear, spline, max-track, leading-edge-track
% Right click: Set cursor point
% Right click and drag: Delete all points in range
% Scroll: Zooms in/out
% Ctrl + any click: Select layer
% Ctrl + any click and drag: Zoom
% Any double click: Nothing
% Ctrl + double click: Zoom reset

obj.tool_list{end+1} = imb.picktool_quality;
addlistener(obj.tool_list{end},'hide_param',@obj.toolparam_close_callback);
% Left click: Nothing
% Left click and drag: Sets all points contains in draw to the currently
%   selected quality level
% Right click: Set cursor point
% Right click and drag: Nothing
% Scroll: Zooms in/out
% Ctrl + any click: Select layer
% Ctrl + any click and drag: Zoom
% Any double click: Nothing
% Ctrl + double click: Zoom reset

obj.tool_list{end+1} = imb.picktool_snake;
addlistener(obj.tool_list{end},'hide_param',@obj.toolparam_close_callback);
% Left click: Enters a point based on parameters
%   Find max in range (specify range line/bin extent to search)
%   Recomputes snake if point is within last range
% Left click and drag: Snakes between manual points based on paramaters
%   Deletes all previous automated points in range
%   Snake tools: basic, crandall, panton
% Right click: Set cursor point
% Right click and drag: Delete all points in range
% Scroll: Zooms in/out
% Ctrl + any click: Select layer
% Ctrl + any click and drag: Zoom
% Any double click: Nothing
% Ctrl + double click: Zoom reset

obj.tool_list{end+1} = imb.picktool_browse;
addlistener(obj.tool_list{end},'hide_param',@obj.toolparam_close_callback);
% Left click: Open A-scope window
% Left click and drag: Nothing
% Right click: Set cursor point
% Right click and drag: Nothing
% Scroll: Zooms in/out
% Ctrl + any click: Select layer
% Ctrl + any click and drag: Zoom
% Any double click: Nothing
% Ctrl + double click: Zoom reset

obj.tool_list{end+1} = imb.picktool_convert;
addlistener(obj.tool_list{end},'hide_param',@obj.toolparam_close_callback);
% Left click and drag: Converts selected layers to the specified layer
%   Deletes all previous points in range
% Right click: Set cursor point
% Right click and drag: Nothing
% Scroll: Zooms in/out
% Ctrl + any click: Select layer
% Ctrl + any click and drag: Zoom
% Any double click: Nothing
% Ctrl + double click: Zoom reset

obj.tool_list{end+1} = imb.picktool_viterbi;
addlistener(obj.tool_list{end},'hide_param',@obj.toolparam_close_callback);
% Left click: Enters a manual point based on parameters
%   Find max in range (specify range line/bin extent to search)
% Left click and drag: Runs tomo.viterbi algorithm on selected data.
% Right click: Set cursor point
% Right click and drag: Delete all points in range
% Scroll: Zooms in/out
% Ctrl + any click: Select layer
% Ctrl + any click and drag: Zoom
% Any double click: Nothing
% Ctrl + double click: Zoom reset

% obj.tool_list{end+1} = imb.picktool_landmark;
% Left click: Create landmark point
% Left click and drag: Create landmark region
% Right click: Set cursor point
% Right click and drag: 
% Scroll: Zooms in/out
% Ctrl + any click: Select layer
% Ctrl + any click and drag: Zoom
% Any double click: Nothing
% Ctrl + double click: Zoom reset

%% Build tool popup menu
obj.left_panel.toolPM = uicontrol('Parent',obj.left_panel.handle);
menuString = {};
for idx = 1:length(obj.tool_list)
  menuString{idx} = obj.tool_list{idx}.tool_name;
end
set(obj.left_panel.toolPM,'String',menuString);
clear menuString;
set(obj.left_panel.toolPM,'Value',1);
tmp = obj.tool_list{1}; obj.left_click = @tmp.left_click;
tmp = obj.tool_list{1}; obj.left_click_and_drag = @tmp.left_click_and_drag;
tmp = obj.tool_list{1}; obj.right_click_and_drag = @tmp.right_click_and_drag;

set(obj.left_panel.toolPM,'Style','popupmenu');
set(obj.left_panel.toolPM,'HorizontalAlignment','Center');
set(obj.left_panel.toolPM,'FontName','fixed');
set(obj.left_panel.toolPM,'Callback',@obj.toolPM_callback);
set(obj.left_panel.toolPM,'TooltipString','Select the active tool');

%--Tool Params Window Pushbutton
obj.left_panel.paramPB = uicontrol('Parent',obj.left_panel.handle);
set(obj.left_panel.paramPB,'Style','PushButton');
set(obj.left_panel.paramPB,'String','Tool Params');
set(obj.left_panel.paramPB,'Callback',@obj.paramPB_callback);
set(obj.left_panel.paramPB,'TooltipString','Open the tool parameters window');

%--Quality Popup Menu
obj.left_panel.qualityPM = uicontrol('Parent',obj.left_panel.handle);
menuString = {};
menuString{1} = 'good';
menuString{2} = 'moderate';
menuString{3} = 'derived';
set(obj.left_panel.qualityPM,'String',menuString);
set(obj.left_panel.qualityPM,'Value',1);
clear menuString;
set(obj.left_panel.qualityPM,'Style','popupmenu');
set(obj.left_panel.qualityPM,'HorizontalAlignment','Center');
set(obj.left_panel.qualityPM,'FontName','fixed');
set(obj.left_panel.qualityPM,'callback',@obj.quality_menu_callback);
set(obj.left_panel.qualityPM,'TooltipString','Set the active quality level');

%--Image Processing Window Pushbutton
if strcmpi(class(obj.h_fig),'double')
  obj.left_panel.imagewin = imagewin(sprintf('%d: Image Params', obj.h_fig), -1, true);
else
  obj.left_panel.imagewin = imagewin(sprintf('%d: Image Params', obj.h_fig.Number), -1, true);
end
obj.left_panel.imagewin.set_auto_caxis(true);
obj.left_panel.imagewin.img_fn_save_fh = @obj.imagewin_fn_callback;
obj.left_panel.imagewin.save_mat_fh = @obj.imagewin_save_mat_callback;
obj.left_panel.imagewin.img_fn_open_fh = @obj.imagewin_fn_callback;
obj.left_panel.imagePB = uicontrol('Parent',obj.left_panel.handle);
set(obj.left_panel.imagePB,'Style','PushButton');
set(obj.left_panel.imagePB,'String','Image Params');
set(obj.left_panel.imagePB,'Callback',@obj.toggle_imagewin_visibility);
set(obj.left_panel.imagePB,'TooltipString','Open image processing window');

%--Y axis choice
obj.left_panel.yaxisPM = uicontrol('Parent',obj.left_panel.handle);
set(obj.left_panel.yaxisPM,'Style','PopupMenu');
set(obj.left_panel.yaxisPM,'String',{'TWTT','WGS84','Range','Range Bin','Surface Flat'});
set(obj.left_panel.yaxisPM,'Value',1);
set(obj.left_panel.yaxisPM,'Callback',@obj.yaxisPM_callback);
set(obj.left_panel.yaxisPM,'TooltipString','Set the y-axis units');

%--X axis choice
obj.left_panel.xaxisPM = uicontrol('Parent',obj.left_panel.handle);
set(obj.left_panel.xaxisPM,'Style','PopupMenu');
set(obj.left_panel.xaxisPM,'String',{'Range Line','Along Track','GPS Time'});
set(obj.left_panel.xaxisPM,'Value',1);
set(obj.left_panel.xaxisPM,'Callback',@obj.xaxisPM_callback);
set(obj.left_panel.xaxisPM,'TooltipString','Set the x-axis units');

%--Max frames selection
obj.left_panel.framesPM = uicontrol('Parent',obj.left_panel.handle);
set(obj.left_panel.framesPM,'Style','popupmenu');
set(obj.left_panel.framesPM,'HorizontalAlignment','Center');
set(obj.left_panel.framesPM,'FontName','fixed');
set(obj.left_panel.framesPM,'String',{'1 Frame','2 Frames','3 Frames','4 Frames','5 Frames','6 Frames','7 Frames','8 Frames','9 Frames','10 Frames'});
set(obj.left_panel.framesPM,'Value',obj.default_params.max_frames);
set(obj.left_panel.framesPM,'Callback',@obj.framesPM_callback);
set(obj.left_panel.framesPM,'TooltipString','Set the number of frames to load/buffer');

%--Crossovers
if strcmpi(class(obj.h_fig),'double')
  obj.eg.crossovers.gui = imb.crossover(sprintf('%d: Crossovers',obj.h_fig), true);
else
  obj.eg.crossovers.gui = imb.crossover(sprintf('%d: Crossovers',obj.h_fig.Number), true);
end
addlistener(obj.eg.crossovers.gui,'update_event',@obj.set_visibility);
addlistener(obj.eg.crossovers.gui,'open_crossover_event',@obj.open_crossover);
addlistener(obj.eg.crossovers.gui,'update_cursor',@obj.cursor_crossover);
addlistener(obj.eg.crossovers.gui,'refresh_crossovers_event',@obj.load_crossovers);
obj.left_panel.crossoverPB = uicontrol('Parent',obj.left_panel.handle);
set(obj.left_panel.crossoverPB,'Style','PushButton');
set(obj.left_panel.crossoverPB,'String','Crossovers');
set(obj.left_panel.crossoverPB,'Callback',@obj.crossoverPB_callback);
set(obj.left_panel.crossoverPB,'TooltipString','Open flightline crossover browse window');

%--save
obj.left_panel.savePB = uicontrol('Parent',obj.left_panel.handle);
set(obj.left_panel.savePB,'Style','PushButton');
set(obj.left_panel.savePB,'String','(S)ave Layer');
set(obj.left_panel.savePB,'Callback',@obj.savePB_callback);
set(obj.left_panel.savePB,'TooltipString','Save layers to database');

% top part table
obj.left_panel.topTable.ui = []; % Parent is a table container
obj.left_panel.topTable.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.left_panel.topTable.height_margin = NaN*zeros(30,30);
obj.left_panel.topTable.false_width = NaN*zeros(30,30);
obj.left_panel.topTable.false_height = NaN*zeros(30,30);

row = 1; col = 1;
obj.left_panel.topTable.handles{row,col}   = obj.left_panel.toolPM;
obj.left_panel.topTable.width(row,col)     = inf;
obj.left_panel.topTable.height(row,col)    = 25;

col = col +1;
obj.left_panel.topTable.handles{row,col}   = obj.left_panel.paramPB;
obj.left_panel.topTable.width(row,col)     = inf;
obj.left_panel.topTable.height(row,col)    = 25;

row = row + 1; col = 1;
obj.left_panel.topTable.handles{row,col}   = obj.left_panel.qualityPM;
obj.left_panel.topTable.width(row,col)     = inf;
obj.left_panel.topTable.height(row,col)    = 25;

col = col + 1;
obj.left_panel.topTable.handles{row,col}   = obj.left_panel.imagePB;
obj.left_panel.topTable.width(row,col)     = inf;
obj.left_panel.topTable.height(row,col)    = 25;

row = row + 1; col = 1;
obj.left_panel.topTable.handles{row,col}   = obj.left_panel.yaxisPM;
obj.left_panel.topTable.width(row,col)     = inf;
obj.left_panel.topTable.height(row,col)    = 25;

col = 2;
obj.left_panel.topTable.handles{row,col}   = obj.left_panel.xaxisPM;
obj.left_panel.topTable.width(row,col)     = inf;
obj.left_panel.topTable.height(row,col)    = 25;

row = row + 1; col = 1;
obj.left_panel.topTable.handles{row,col}   = obj.left_panel.framesPM;
obj.left_panel.topTable.width(row,col)     = inf;
obj.left_panel.topTable.height(row,col)    = 25;

col = 2;
obj.left_panel.topTable.handles{row,col}   = obj.left_panel.crossoverPB;
obj.left_panel.topTable.width(row,col)     = inf;
obj.left_panel.topTable.height(row,col)    = 25;

row = row + 1; col = 1;
obj.left_panel.topTable.handles{row,col}   = obj.left_panel.savePB;
obj.left_panel.topTable.width(row,col)     = inf;
obj.left_panel.topTable.height(row,col)    = 25;

clear row col
%-----------------------

%---- layers list box
obj.left_panel.layerLB_panel_handle = uipanel('Parent',obj.left_panel.handle);
set(obj.left_panel.layerLB_panel_handle,'bordertype','none');
obj.layerLB_init();

%----frames list box
obj.left_panel.frameLB = uicontrol('Parent',obj.left_panel.handle);
set(obj.left_panel.frameLB,'Style','listbox');
set(obj.left_panel.frameLB,'HorizontalAlignment','Center');
set(obj.left_panel.frameLB,'FontName','fixed');
set(obj.left_panel.frameLB,'Callback',@obj.frameLB_callback);
obj.left_panel.frameCM = uicontextmenu('Parent',obj.h_fig);
% Define the context menu items and install their callbacks
uimenu(obj.left_panel.frameCM, 'Label', 'Copy', 'Callback', @obj.frameCM_callback);
set(obj.left_panel.frameLB,'UIContextMenu',obj.left_panel.frameCM);
set(obj.left_panel.frameLB,'TooltipString','Load frame');

%---- Source data Listbox
obj.left_panel.sourceLB = uicontrol('Parent',obj.left_panel.handle);
set(obj.left_panel.sourceLB,'Style','listbox');
set(obj.left_panel.sourceLB,'HorizontalAlignment','Center');
set(obj.left_panel.sourceLB,'FontName','fixed');
set(obj.left_panel.sourceLB,'Callback',@obj.sourceLB_callback);
set(obj.left_panel.sourceLB,'TooltipString','Echogram data source to load');
obj.left_panel.sourceCM = uicontextmenu('Parent',obj.h_fig);
% Define the context menu items and install their callbacks
uimenu(obj.left_panel.sourceCM, 'Label', 'Add', 'Callback', @obj.sourceCM_callback);
uimenu(obj.left_panel.sourceCM, 'Label', 'Remove', 'Callback', @obj.sourceCM_callback);
uimenu(obj.left_panel.sourceCM, 'Label', 'Refresh', 'Callback', @obj.sourceCM_callback);
set(obj.left_panel.sourceLB,'uicontextmenu',obj.left_panel.sourceCM);

%----------------------------------------------------------
% set up left_panel table
obj.left_panel.table.ui = obj.left_panel.handle;
obj.left_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.left_panel.table.height_margin = NaN*zeros(30,30);
obj.left_panel.table.false_width = NaN*zeros(30,30);
obj.left_panel.table.false_height = NaN*zeros(30,30);
obj.left_panel.table.offset = [0 0];

row = 1; col = 1;
obj.left_panel.table.handles{row,col}   = obj.left_panel.topTable;
obj.left_panel.table.width(row,col)     = inf;
obj.left_panel.table.height(row,col)    = 130;
obj.left_panel.table.width_margin(row,col) = 3;
obj.left_panel.table.height_margin(row,col) = 3;

row = 2; col = 1;
obj.left_panel.table.handles{row,col}   = obj.left_panel.layerLB_panel_handle;
obj.left_panel.table.width(row,col)     = inf;
obj.left_panel.table.height(row,col)    = 71;
obj.left_panel.table.width_margin(row,col) = 3;
obj.left_panel.table.height_margin(row,col) = 3;

row = 3; col = 1;
obj.left_panel.table.handles{row,col}   = obj.left_panel.frameLB;
obj.left_panel.table.width(row,col)     = inf;
obj.left_panel.table.height(row,col)    = inf;
obj.left_panel.table.width_margin(row,col) = 3;
obj.left_panel.table.height_margin(row,col) = 3;

row = 4; col = 1;
obj.left_panel.table.handles{row,col}   = obj.left_panel.sourceLB;
obj.left_panel.table.width(row,col)     = inf;
obj.left_panel.table.height(row,col)    = 75;
obj.left_panel.table.width_margin(row,col) = 3;
obj.left_panel.table.height_margin(row,col) = 9;
clear row col

% Draw table
table_draw(obj.left_panel.table);


%==========================================================
% right_panel contents: echogram and statusbar panel

% create axes panel and axes
obj.right_panel.axes.panel_h = uipanel('parent',obj.right_panel.handle);
obj.right_panel.axes.handle = axes('parent',obj.right_panel.axes.panel_h);
axis_pos = get(obj.right_panel.axes.handle,'Position');
% scale up the size of the axes (image part)
% about 87% width of panel and 87% height of panel
xy_diff = axis_pos(1)-axis_pos(2);
x_delta = axis_pos(1)-0.08;
y_delta = axis_pos(2)-0.10+xy_diff;
new_pos = [axis_pos(1)-x_delta axis_pos(2)-y_delta axis_pos(3)+2*x_delta axis_pos(4)+2*y_delta];
set(obj.right_panel.axes.handle,'Position',new_pos);
set(obj.right_panel.axes.handle,'activepositionproperty','outerposition');
set(obj.h_fig,'Colormap',1-gray(256));

% Create imagesc handle
hold(obj.right_panel.axes.handle,'on');
obj.eg.h_image = imagesc([],[],[],'parent',obj.right_panel.axes.handle);
obj.left_panel.imagewin.set_img(obj.eg.h_image);

% Create cursor plot handle
obj.cursor.h = plot(NaN,NaN,'k--','parent',obj.right_panel.axes.handle);

% create statusbar panel
obj.right_panel.status_panel.handle = uipanel('Parent',obj.right_panel.handle);

% create statusbar's elements: cmd info, cursor info, mouse coordinate info
% command info display
obj.right_panel.status_panel.statusText = uicontrol('parent',obj.right_panel.status_panel.handle);
set(obj.right_panel.status_panel.statusText,'Style','text');
%set(obj.right_panel.status_panel.statusText,'FontName','courier');
if isunix
  set(obj.right_panel.status_panel.statusText,'FontSize',9)
else
  set(obj.right_panel.status_panel.statusText,'FontSize',8)
end
set(obj.right_panel.status_panel.statusText,'HorizontalAlignment','left');
set(obj.right_panel.status_panel.statusText,'String','');

% mouse coordinate info display
obj.right_panel.status_panel.mouseCoordText = uicontrol('parent',obj.right_panel.status_panel.handle);
set(obj.right_panel.status_panel.mouseCoordText,'Style','text');
%set(obj.right_panel.status_panel.mouseCoordText,'FontName','courier');
if isunix
  set(obj.right_panel.status_panel.mouseCoordText,'FontSize',9)
else
  set(obj.right_panel.status_panel.mouseCoordText,'FontSize',8)
end
set(obj.right_panel.status_panel.mouseCoordText,'HorizontalAlignment','left');
set(obj.right_panel.status_panel.mouseCoordText,'String','');

%----echogram context menu
obj.right_panel.echoCM= uicontextmenu('parent',obj.h_fig);
% Define the context menu items and install their callbacks
obj.right_panel.echoCM_item1 = uimenu(obj.right_panel.echoCM, 'Label', 'Copy Cursor Info (Ctrl-C)', 'Callback',@obj.status_text_copy_callback);
set(obj.right_panel.status_panel.statusText,'uicontextmenu',obj.right_panel.echoCM);

% create statusbar's table
obj.right_panel.status_panel.table.ui = obj.right_panel.status_panel.handle;
obj.right_panel.status_panel.table.offset = [0 2];
row = 1; col = 1;
obj.right_panel.status_panel.table.handles{row,col}       = obj.right_panel.status_panel.statusText;
obj.right_panel.status_panel.table.width(row,col)         = inf;
obj.right_panel.status_panel.table.height(row,col)        = inf;
obj.right_panel.status_panel.table.width_margin(row,col)  = 0;
obj.right_panel.status_panel.table.height_margin(row,col) = 0;
row = 1; col = 2;
obj.right_panel.status_panel.table.handles{row,col}       = obj.right_panel.status_panel.mouseCoordText;
obj.right_panel.status_panel.table.width(row,col)         = 200;
obj.right_panel.status_panel.table.height(row,col)        = inf;
obj.right_panel.status_panel.table.width_margin(row,col)  = 0;
obj.right_panel.status_panel.table.height_margin(row,col) = 0;
clear row col
table_draw(obj.right_panel.status_panel.table);

% make the right panel's table 
obj.right_panel.table.ui=obj.right_panel.handle;
obj.right_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.right_panel.table.height_margin = NaN*zeros(30,30);
obj.right_panel.table.false_width = NaN*zeros(30,30);
obj.right_panel.table.false_height = NaN*zeros(30,30);
obj.right_panel.table.offset = [0 0];
row = 1; col = 1;
obj.right_panel.table.handles{row,col}   = obj.right_panel.axes.panel_h;
obj.right_panel.table.width(row,col)     = inf;
obj.right_panel.table.height(row,col)    = inf;
obj.right_panel.table.width_margin(row,col)=0;
obj.right_panel.table.height_margin(row,col)=0;
row = 2; col = 1;
obj.right_panel.table.handles{row,col}   = obj.right_panel.status_panel.handle;
obj.right_panel.table.width(row,col)     = inf;
obj.right_panel.table.height(row,col)    = 15;
obj.right_panel.table.width_margin(row,col)=0;
obj.right_panel.table.height_margin(row,col)=0;
obj.right_panel.table.false_height(row,col)=0;
clear row col
table_draw(obj.right_panel.table);

return;
