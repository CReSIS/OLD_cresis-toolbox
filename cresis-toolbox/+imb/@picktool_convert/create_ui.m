function create_ui(obj)

set(obj.h_fig,'DockControls','off')
set(obj.h_fig,'NumberTitle','off');
if strcmpi(class(obj.h_fig),'double')
  set(obj.h_fig,'Name',sprintf('%d: %s tool parameters', obj.h_fig, obj.tool_name_title));
else
  set(obj.h_fig,'Name',sprintf('%d: %s tool parameters', obj.h_fig.Number, obj.tool_name_title));
end
set(obj.h_fig,'ToolBar','none');
set(obj.h_fig,'MenuBar','none');
%set(obj.h_fig,'Resize','off');
set(obj.h_fig,'KeyPressFcn',@obj.key_press);
set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
set(obj.h_fig,'Position',[0 0 obj.w obj.h]);
%==========================================================================
% top panel
obj.top_panel.handle = uipanel('Parent',obj.h_fig);
set(obj.top_panel.handle,'HighlightColor',[0.8 0.8 0.8]);
set(obj.top_panel.handle,'ShadowColor',[0.6 0.6 0.6]);

obj.bottom_panel.handle = uipanel('Parent',obj.h_fig);
set(obj.bottom_panel.handle,'HighlightColor',[0.8 0.8 0.8]);
set(obj.bottom_panel.handle,'ShadowColor',[0.6 0.6 0.6]);
%--------------------------------------
% table
obj.table.ui=obj.h_fig;
row = 1; col = 1;
obj.table.handles{row,col}   = obj.top_panel.handle;
obj.table.width(row,col)     = inf;
obj.table.height(row,col)    = inf;
obj.table.width_margin(row,col) = 0;
obj.table.height_margin(row,col) = 0;
row = 2; col = 1;
obj.table.handles{row,col}   = obj.bottom_panel.handle;
obj.table.width(row,col)     = inf;
obj.table.height(row,col)    = inf;
obj.table.width_margin(row,col) = 0;
obj.table.height_margin(row,col) = 0;
clear row col
table_draw(obj.table);
%============================================================================================
% top panel table contents
%----Target layer type
obj.top_panel.target_layers_name = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.target_layers_name,'Style','text');
set(obj.top_panel.target_layers_name,'FontSize',10);
set(obj.top_panel.target_layers_name,'String','Source layer:');

%----Target layer input
obj.top_panel.target_layers_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.target_layers_TE,'Style','edit');
set(obj.top_panel.target_layers_TE,'String','1');
%---------------------------------------------------------------------------------------------
% set up top panel table
obj.top_panel.table.ui=obj.top_panel.handle;
obj.top_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.top_panel.table.height_margin = NaN*zeros(30,30);
obj.top_panel.table.false_width = NaN*zeros(30,30);
obj.top_panel.table.false_height = NaN*zeros(30,30);
obj.top_panel.table.offset = [0 0];


row = 1; col = 1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.target_layers_name;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 1; col = 2;
obj.top_panel.table.handles{row,col}   = obj.top_panel.target_layers_TE;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

clear row col

% Draw table
table_draw(obj.top_panel.table);

%-----Info message
obj.bottom_panel.info_label = uicontrol('Parent',obj.bottom_panel.handle);
set(obj.bottom_panel.info_label,'Style','text');
set(obj.bottom_panel.info_label,'String','Enter the number of the layers you want to write into the selected layers.');

%---------------------------------------------------------------------------------------------
% set up top panel table
obj.bottom_panel.table.ui=obj.bottom_panel.handle;
obj.bottom_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.bottom_panel.table.height_margin = NaN*zeros(30,30);
obj.bottom_panel.table.false_width = NaN*zeros(30,30);
obj.bottom_panel.table.false_height = NaN*zeros(30,30);
obj.bottom_panel.table.offset = [0 0];

row = 1; col = 1;
obj.bottom_panel.table.handles{row,col}   = obj.bottom_panel.info_label;
obj.bottom_panel.table.width(row,col)     = inf;
obj.bottom_panel.table.height(row,col)    = inf;
obj.bottom_panel.table.width_margin(row,col)= 1.5;
obj.bottom_panel.table.height_margin(row,col)=1.5;

clear row col

table_draw(obj.bottom_panel.table);


return

