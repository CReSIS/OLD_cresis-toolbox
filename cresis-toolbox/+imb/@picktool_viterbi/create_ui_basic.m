function create_ui_basic(obj,xpos,ypos)

% create_ui_basic(obj,xpos,ypos)
%
% Creates components for the HMM param window's UI when the Viterbi
% tool is selected. Plots the window at xpos,ypos.
%

set(obj.h_fig,'visible','off');

% set default position (changed when window accessed)
set(obj.h_fig,'Units','Pixels');
set(obj.h_fig,'Position',[xpos ypos obj.w obj.h]);
% show top panel 
% set(obj.top_panel.handle,'visible','on');
% set(obj.bottom_panel.handle,'visible','off');

if ~obj.first_time
  figure(obj.h_fig);
  clf;
  obj.table = [];
end

%==========================================================================
% top panel
obj.top_panel.handle = uipanel('Parent',obj.h_fig);
set(obj.top_panel.handle,'HighlightColor',[0.8 0.8 0.8]);
set(obj.top_panel.handle,'ShadowColor',[0.6 0.6 0.6]);
%set(obj.top_panel.handle,'visible','off');

%--------------------------------------
% table
obj.table.ui=obj.h_fig;

obj.table.handles{1,1}        = obj.top_panel.handle;
obj.table.width(1,1)          = inf;
obj.table.height(1,1)         = inf;
obj.table.width_margin(1,1)   = 0;
obj.table.height_margin(1,1)  = 0;

table_draw(obj.table);

%============================================================================================
% top panel table contents

%----HMM tool list box
obj.top_panel.tool_PM = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.tool_PM,'Style','popupmenu');
set(obj.top_panel.tool_PM,'String',{'basic'});
set(obj.top_panel.tool_PM,'Value',1)
set(obj.top_panel.tool_PM,'Callback',@obj.toolPM_callback);

%-----mode label
obj.top_panel.mode_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.mode_label,'Style','text');
set(obj.top_panel.mode_label,'String','Mode');
%----insert range
obj.top_panel.insert_range_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.insert_range_label,'Style','text');
set(obj.top_panel.insert_range_label,'String','Max point range:');
%----insert pt search range box
obj.top_panel.insert_range_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.insert_range_TE,'Style','edit');
set(obj.top_panel.insert_range_TE,'String',obj.in_rng_sv);
%----column restriction label
obj.top_panel.column_restriction_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.column_restriction_label,'Style','text');
set(obj.top_panel.column_restriction_label,'String','Column tracking restriction:');
%----column restriction cbox
obj.top_panel.column_restriction_cbox = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.column_restriction_cbox,'Style','checkbox');
set(obj.top_panel.column_restriction_cbox,'Value', 1);
%----top suppression label
obj.top_panel.top_sup_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.top_sup_label,'Style','text');
set(obj.top_panel.top_sup_label,'String',sprintf('Top\nsuppression:'));
%----top suppression cbox
obj.top_panel.top_sup_cbox = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.top_sup_cbox,'Style','checkbox');
set(obj.top_panel.top_sup_cbox,'Value', 1);
%----multiple suppression label
obj.top_panel.mult_sup_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.mult_sup_label,'Style','text');
set(obj.top_panel.mult_sup_label,'String','Multiple suppression:');
%----multiple suppression cbox
obj.top_panel.mult_sup_cbox = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.mult_sup_cbox,'Style','checkbox');
set(obj.top_panel.mult_sup_cbox,'Value', 1);
%%
%---------------------------------------------------------------------------------------------
rows = 5;  % Update with number of rows and columns
cols = 2;
% set up top panel table
default_dimensions = NaN*zeros(rows,cols);
obj.top_panel.table.ui=obj.top_panel.handle;
obj.top_panel.table.width_margin = default_dimensions;
obj.top_panel.table.height_margin = default_dimensions;
obj.top_panel.table.false_width = default_dimensions;
obj.top_panel.table.false_height = default_dimensions;
obj.top_panel.table.offset = [0 0];

obj.top_panel.table.width = ones(rows, cols) * inf;
obj.top_panel.table.height = ones(rows, cols) * inf;
obj.top_panel.table.width_margin = ones(rows, cols) * 1.5;
obj.top_panel.table.height_margin = ones(rows, cols) * 1.5;

%% Mode
obj.top_panel.table.handles{1,1}   = obj.top_panel.mode_label;
obj.top_panel.table.handles{1,2}   = obj.top_panel.tool_PM;
%% Insert Range
obj.top_panel.table.handles{2,1}   = obj.top_panel.insert_range_label;
obj.top_panel.table.handles{2,2}   = obj.top_panel.insert_range_TE;
%% Column restriction
obj.top_panel.table.handles{3,1}   = obj.top_panel.column_restriction_label;
obj.top_panel.table.handles{3,2}   = obj.top_panel.column_restriction_cbox;
%% Top suppression
obj.top_panel.table.handles{4,1}   = obj.top_panel.top_sup_label;
obj.top_panel.table.handles{4,2}   = obj.top_panel.top_sup_cbox;
%% Multiple suppression
obj.top_panel.table.handles{5,1}   = obj.top_panel.mult_sup_label;
obj.top_panel.table.handles{5,2}   = obj.top_panel.mult_sup_cbox;
% TODO[reece]: Add number inputs for all weights
clear rows cols

% Draw table
table_draw(obj.top_panel.table);

if obj.first_time
  obj.first_time = false;
else
  set(obj.h_fig,'visible','on');
end

return;