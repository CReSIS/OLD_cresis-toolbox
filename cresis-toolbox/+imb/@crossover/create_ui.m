function create_ui(obj)
% crossover.create_ui(obj)
%
% Create the cross over figure and GUI

%% Create Figure

if obj.hide_only
  obj.h_fig = figure('Visible','off');
else
  obj.h_fig = figure;
end

set(obj.h_fig,'DockControls','off');
set(obj.h_fig,'NumberTitle','off');
set(obj.h_fig,'Name',obj.img_title);
set(obj.h_fig,'ToolBar','none');
set(obj.h_fig,'MenuBar','none'); 
set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
set(obj.h_fig,'Units','Points');
pos = get(obj.h_fig,'Position');
set(obj.h_fig,'Position',[pos(1:2) 200 240]);

%% Create GUI Objects

obj.h_gui.visibleCB = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.visibleCB,'Style','CheckBox');
set(obj.h_gui.visibleCB,'String','Show crossovers');
set(obj.h_gui.visibleCB,'Callback',@obj.update_callback);

obj.h_gui.resetPB = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.resetPB,'Style','Pushbutton');
set(obj.h_gui.resetPB,'String','Reset');
set(obj.h_gui.resetPB,'Callback',@obj.reset);

obj.h_gui.sortPM = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.sortPM,'Style','popupmenu');
set(obj.h_gui.sortPM,'String',{'No sort','Error','Angle','Frame','Time'});
set(obj.h_gui.sortPM,'Value',1);
set(obj.h_gui.sortPM,'Callback',@obj.sortPM_callback);

obj.h_gui.sort_orderCB = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.sort_orderCB,'Style','CheckBox');
set(obj.h_gui.sort_orderCB,'String','Descending');
set(obj.h_gui.sort_orderCB,'Callback',@obj.sortPM_callback);

obj.h_gui.crossoverLB  = listbox_mask(obj.h_fig,'',{},[],[],1);
addlistener(obj.h_gui.crossoverLB,'list_changed',@obj.update_callback);

obj.h_gui.crossoverCM = uicontextmenu('Parent',obj.h_fig);
% Define the context menu items and install their callbacks
uimenu(obj.h_gui.crossoverCM, 'Label', 'Set Cursor', 'Callback', @obj.set_cursor);
uimenu(obj.h_gui.crossoverCM, 'Label', 'Refresh Crossovers', 'Callback', @obj.refresh_crossovers);
uimenu(obj.h_gui.crossoverCM, 'Label', 'Copy', 'Callback', @obj.copy_crossover);
set(obj.h_gui.crossoverLB.h_list,'uicontextmenu',obj.h_gui.crossoverCM);

obj.h_gui.nanCB = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.nanCB,'Style','CheckBox');
set(obj.h_gui.nanCB,'String','Show NaN');
set(obj.h_gui.nanCB,'Value',true);
set(obj.h_gui.nanCB,'Callback',@obj.update_callback);

obj.h_gui.angleText = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.angleText,'Style','Text');
set(obj.h_gui.angleText,'String','Min/Max Angle');
      
obj.h_gui.minAngleLE = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.minAngleLE,'Style','Edit');
set(obj.h_gui.minAngleLE,'String','0');
set(obj.h_gui.minAngleLE,'Callback',@obj.update_callback);
      
obj.h_gui.maxAngleLE = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.maxAngleLE,'Style','Edit');
set(obj.h_gui.maxAngleLE,'String','inf');
set(obj.h_gui.maxAngleLE,'Callback',@obj.update_callback);

obj.h_gui.errorText = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.errorText,'Style','Text');
set(obj.h_gui.errorText,'String','Min/Max Error');
      
obj.h_gui.minErrorLE = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.minErrorLE,'Style','Edit');
set(obj.h_gui.minErrorLE,'String','0');
set(obj.h_gui.minErrorLE,'Callback',@obj.update_callback);
      
obj.h_gui.maxErrorLE = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.maxErrorLE,'Style','Edit');
set(obj.h_gui.maxErrorLE,'String','inf');
set(obj.h_gui.maxErrorLE,'Callback',@obj.update_callback);

obj.h_gui.qualityText = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.qualityText,'Style','Text');
set(obj.h_gui.qualityText,'String','Quality');
      
obj.h_gui.qualityLE = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.qualityText,'Style','Edit');
set(obj.h_gui.qualityText,'String','[1 2 3]');
set(obj.h_gui.qualityText,'Callback',@obj.update_callback);

%% Create GUI Table
obj.h_gui.table.ui = obj.h_fig;
obj.h_gui.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.h_gui.table.height_margin = NaN*zeros(30,30);
obj.h_gui.table.false_width = NaN*zeros(30,30);
obj.h_gui.table.false_height = NaN*zeros(30,30);
obj.h_gui.table.offset = [0 0];

row = 1;
col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.visibleCB;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

col = 2;
obj.h_gui.table.handles{row,col}   = obj.h_gui.resetPB;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

row = row + 1;
col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.sortPM;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
col = 2;
obj.h_gui.table.handles{row,col}   = obj.h_gui.sort_orderCB;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;

row = row + 1;
col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.crossoverLB.h_list;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = inf;
col = 2;
obj.h_gui.table.width(row,col)     = 0;
obj.h_gui.table.height(row,col)    = 20;

row = row + 1;
col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.crossoverLB.h_LE;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
col = 2;
obj.h_gui.table.handles{row,col}   = obj.h_gui.nanCB;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;

row = row + 1;
col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.angleText;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
col = 2;
obj.h_gui.table.width(row,col)     = 0;
obj.h_gui.table.height(row,col)    = 20;

row = row + 1;
col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.minAngleLE;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
col = 2;
obj.h_gui.table.handles{row,col}   = obj.h_gui.maxAngleLE;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;

row = row + 1;
col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.errorText;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
col = 2;
obj.h_gui.table.width(row,col)     = 0;
obj.h_gui.table.height(row,col)    = 20;

row = row + 1;
col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.minErrorLE;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
col = 2;
obj.h_gui.table.handles{row,col}   = obj.h_gui.maxErrorLE;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;

row = row + 1;
col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.qualityText;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
col = 2;
obj.h_gui.table.handles{row,col}   = obj.h_gui.qualityLE;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
clear row col
table_draw(obj.h_gui.table);

end
