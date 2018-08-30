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

row = 1; col = 1;
obj.table.handles{row,col}   = obj.top_panel.handle;
obj.table.width(row,col)     = inf;
obj.table.height(row,col)    = inf;
obj.table.width_margin(row,col) = 0;
obj.table.height_margin(row,col) = 0;

clear row col
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
%----Viterbi search range name
obj.top_panel.Viterbi_range_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.Viterbi_range_label,'Style','text');
set(obj.top_panel.Viterbi_range_label,'String','HMM detection range:');
%----Viterbi search range box
obj.top_panel.Viterbi_range_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.Viterbi_range_TE,'Style','edit');
set(obj.top_panel.Viterbi_range_TE,'String',obj.sn_rng_sv);
%----reinterp mode enable label
obj.top_panel.reinterp_mode_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.reinterp_mode_label,'Style','text');
set(obj.top_panel.reinterp_mode_label,'String','Reinterpolation (beta):');
%----reinterp mode enable cbox
obj.top_panel.reinterp_mode_cbox = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.reinterp_mode_cbox,'Style','checkbox');
%----column restriction label
obj.top_panel.column_restriction_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.column_restriction_label,'Style','text');
set(obj.top_panel.column_restriction_label,'String','Column tracking restriction:');
%----column restriction cbox
obj.top_panel.column_restriction_cbox = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.column_restriction_cbox,'Style','checkbox');
set(obj.top_panel.column_restriction_cbox,'Value', 1);
%%
%---------------------------------------------------------------------------------------------
% set up top panel table
obj.top_panel.table.ui=obj.top_panel.handle;
obj.top_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.top_panel.table.height_margin = NaN*zeros(30,30);
obj.top_panel.table.false_width = NaN*zeros(30,30);
obj.top_panel.table.false_height = NaN*zeros(30,30);
obj.top_panel.table.offset = [0 0];

row = 1; col = 1; 
obj.top_panel.table.handles{row,col}   = obj.top_panel.mode_label;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 1; col = 2;
obj.top_panel.table.handles{row,col}   = obj.top_panel.tool_PM;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 2; col = 1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.insert_range_label;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 2; col = 2; 
obj.top_panel.table.handles{row,col}   = obj.top_panel.insert_range_TE;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 3; col = 1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.Viterbi_range_label;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 3; col = 2; 
obj.top_panel.table.handles{row,col}   = obj.top_panel.Viterbi_range_TE;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 4; col = 1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.reinterp_mode_label;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 4; col = 2; 
obj.top_panel.table.handles{row,col}   = obj.top_panel.reinterp_mode_cbox;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 5; col = 1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.column_restriction_label;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 5; col = 2; 
obj.top_panel.table.handles{row,col}   = obj.top_panel.column_restriction_cbox;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

clear row col

% Draw table
table_draw(obj.top_panel.table);

if obj.first_time
  obj.first_time = false;
else
  set(obj.h_fig,'visible','on');
end

return;