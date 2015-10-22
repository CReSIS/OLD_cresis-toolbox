function layerLB_init(obj)

% layerLB_int(obj)
%
% Layer listbox initialization function. Looks in the echowin object passed
% in for obj.left_panel.layerLB_panel_handle and plots the listbox into
% that panel. Might be better to determine this panel via a parameter
% passed in the function call...
%
% This function initializes only UI elements of the listbox. See function
% layerLB_setdata for setting the value of initial layer names as well as
% visibility and selection state variables.
%

% "constant" parameters
% number of rows that will simultaneously be displayed by the listbox
obj.left_panel.layer_panel.MAX_ROW = 5;

%==========================================================================
% top obj.left_panel.layer_panel
obj.left_panel.layer_panel.handle = uipanel('Parent',obj.left_panel.layerLB_panel_handle);
set(obj.left_panel.layer_panel.handle,'HighlightColor',[0.8 0.8 0.8]);
set(obj.left_panel.layer_panel.handle,'ShadowColor',[0.6 0.6 0.6]);
set(obj.left_panel.layer_panel.handle,'Units','pixels');
%set(obj.left_panel.layer_panel.handle,'Position',0,0,180,60);
%set(obj.left_panel.layer_panel.handle,'visible','off');

%==========================================================================
% slider
obj.left_panel.layer_panel.slider = uicontrol('Parent',obj.left_panel.layerLB_panel_handle);
set(obj.left_panel.layer_panel.slider,'Style','slider');
set(obj.left_panel.layer_panel.slider,'Enable','off');
set(obj.left_panel.layer_panel.slider,'Min',0-eps);
set(obj.left_panel.layer_panel.slider,'Max',0);
set(obj.left_panel.layer_panel.slider,'Value',0);
set(obj.left_panel.layer_panel.slider,'Callback',@obj.layerLB_slider_callback);

%--------------------------------------
% base table
table.ui=obj.left_panel.layerLB_panel_handle;

row = 1; col = 1;
table.handles{row,col}       = obj.left_panel.layer_panel.handle;
table.width(row,col)         = inf;
table.height(row,col)        = inf;
table.width_margin(row,col)  = 0;
table.height_margin(row,col) = 0;

row = 1; col = 2;
table.handles{row,col}       = obj.left_panel.layer_panel.slider;
table.width(row,col)         = 12.5;
table.height(row,col)        = inf;
table.width_margin(row,col)  = 0;
table.height_margin(row,col) = 0;

clear row col
table_draw(table);

% configure base table

obj.left_panel.layer_panel.table.ui=obj.left_panel.layer_panel.handle;
obj.left_panel.layer_panel.table.handles = [];
obj.left_panel.layer_panel.table.offset = [0 0];

% build the table of visible layers
for idx = 1:obj.left_panel.layer_panel.MAX_ROW
  % create obj.left_panel.layer_panel components
  %----ELEMENT descriptor
  element = uicontrol('Parent',obj.left_panel.layer_panel.handle);
  set(element,'Style','text');
  if isunix
    set(element,'FontSize',10)
    set(element,'FontName','fixed');
  else
    set(element,'FontSize',8)
    set(element,'FontName','courier');
  end
  set(element,'HorizontalAlignment','left');
  
  %----ENABLE checkbox
  %----All enabled on initialization
  cbox = uicontrol('Parent',obj.left_panel.layer_panel.handle);
  set(cbox,'Style','checkbox');
  set(cbox,'Callback',@obj.layerLB_check_callback);
  
  %----EDITING radio button
  %----All disabled on initialization
  rbut = uicontrol('Parent',obj.left_panel.layer_panel.handle);
  set(rbut,'Style','radiobutton');
  set(rbut,'Callback',@obj.layerLB_radio_callback);
  
  % create table
  row = idx; col = 1;
  obj.left_panel.layer_panel.table.handles{row,col}   = cbox;
  obj.left_panel.layer_panel.table.width(row,col)     = 20;
  obj.left_panel.layer_panel.table.height(row,col)    = 12;
  obj.left_panel.layer_panel.table.width_margin(row,col) = 2;
  obj.left_panel.layer_panel.table.height_margin(row,col) = 1;
  
  col = 2;
  obj.left_panel.layer_panel.table.handles{row,col}   = rbut;
  obj.left_panel.layer_panel.table.width(row,col)     = 10;
  obj.left_panel.layer_panel.table.height(row,col)    = 12;
  obj.left_panel.layer_panel.table.width_margin(row,col) = 0;
  obj.left_panel.layer_panel.table.height_margin(row,col) = 1;
  
  col = 3;
  obj.left_panel.layer_panel.table.handles{row,col}   = element;
  obj.left_panel.layer_panel.table.width(row,col)     = inf;
  obj.left_panel.layer_panel.table.height(row,col)    = 12;
  obj.left_panel.layer_panel.table.width_margin(row,col) = 5;
  obj.left_panel.layer_panel.table.height_margin(row,col) = 1.5;
  
  set(element,'Visible','off');
  set(cbox,'Visible','off');
  set(rbut,'Visible','off');
end

% Add a row of NaN handles to avoid table_draw messing up alignment
row = obj.left_panel.layer_panel.MAX_ROW+1; col = 1;
obj.left_panel.layer_panel.table.handles{row,col}   = NaN;
obj.left_panel.layer_panel.table.width(row,col)     = 25;
obj.left_panel.layer_panel.table.height(row,col)    = 12;
obj.left_panel.layer_panel.table.width_margin(row,col) = 5;
obj.left_panel.layer_panel.table.height_margin(row,col) = 1;

col = 2;
obj.left_panel.layer_panel.table.handles{row,col}   = NaN;
obj.left_panel.layer_panel.table.width(row,col)     = 15;
obj.left_panel.layer_panel.table.height(row,col)    = 12;
obj.left_panel.layer_panel.table.width_margin(row,col) = 0;
obj.left_panel.layer_panel.table.height_margin(row,col) = 1;

col = 3;
obj.left_panel.layer_panel.table.handles{row,col}   = NaN;
obj.left_panel.layer_panel.table.width(row,col)     = 140;
obj.left_panel.layer_panel.table.height(row,col)    = 12;
obj.left_panel.layer_panel.table.width_margin(row,col) = 5;
obj.left_panel.layer_panel.table.height_margin(row,col) = 1.5;

clear row col
table_draw(obj.left_panel.layer_panel.table);
