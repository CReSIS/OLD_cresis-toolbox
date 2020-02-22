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
%----quality output label
obj.top_panel.quality_output_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.quality_output_label,'Style','text');
set(obj.top_panel.quality_output_label,'String',sprintf('Quality\noutput:'));
%----quality output cbox
obj.top_panel.quality_output_cbox = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.quality_output_cbox,'Style','checkbox');
set(obj.top_panel.quality_output_cbox,'Value', 0);
%----quality threshold
obj.top_panel.quality_threshold_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.quality_threshold_label,'Style','text');
set(obj.top_panel.quality_threshold_label,'String',sprintf('Quality\nthreshold:'));
%----quality threshold box
obj.top_panel.quality_threshold_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.quality_threshold_TE,'Style','edit');
set(obj.top_panel.quality_threshold_TE,'String','-20');
%----smoothness weight label
obj.top_panel.smoothness_weight_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.smoothness_weight_label,'Style','text');
set(obj.top_panel.smoothness_weight_label,'String',sprintf('Smoothness\nweight:'));
%----smoothness weight box
obj.top_panel.smoothness_weight_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.smoothness_weight_TE,'Style','edit');
set(obj.top_panel.smoothness_weight_TE,'String','3');
%----smoothness variance label
obj.top_panel.smoothness_variance_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.smoothness_variance_label,'Style','text');
set(obj.top_panel.smoothness_variance_label,'String',sprintf('Smoothness\nvariance:'));
%----smoothness variance box
obj.top_panel.smoothness_variance_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.smoothness_variance_TE,'Style','edit');
set(obj.top_panel.smoothness_variance_TE,'String','Inf');
%----repulsion label
obj.top_panel.repulsion_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.repulsion_label,'Style','text');
set(obj.top_panel.repulsion_label,'String',sprintf('Repulsion:'));
%----repulsion box
obj.top_panel.repulsion_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.repulsion_TE,'Style','edit');
set(obj.top_panel.repulsion_TE,'String','150000');
%----ice bin threshold label
obj.top_panel.icebinthr_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.icebinthr_label,'Style','text');
set(obj.top_panel.icebinthr_label,'String',sprintf('Ice bin\nthreshold:'));
%----ice bin threshold box
obj.top_panel.icebinthr_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.icebinthr_TE,'Style','edit');
set(obj.top_panel.icebinthr_TE,'String','10');
%%
%---------------------------------------------------------------------------------------------
% set up top panel table
obj.top_panel.table.ui=obj.top_panel.handle;
obj.top_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.top_panel.table.height_margin = NaN*zeros(30,30);
obj.top_panel.table.false_width = NaN*zeros(30,30);
obj.top_panel.table.false_height = NaN*zeros(30,30);
obj.top_panel.table.offset = [0 0];

rows = 13;
cols = 2;
obj.top_panel.table.width = ones(rows, cols) * inf;
obj.top_panel.table.height = ones(rows, cols) * inf;
obj.top_panel.table.width_margin = ones(rows, cols) * 1.5;
obj.top_panel.table.height_margin = ones(rows, cols) * 1.5;

% TODO[reece]: Remove unused and add weights
obj.top_panel.table.handles{1,1}   = obj.top_panel.mode_label;
obj.top_panel.table.handles{1,2}   = obj.top_panel.tool_PM;
obj.top_panel.table.handles{2,1}   = obj.top_panel.insert_range_label;
obj.top_panel.table.handles{2,2}   = obj.top_panel.insert_range_TE;
obj.top_panel.table.handles{3,1}   = obj.top_panel.Viterbi_range_label;
obj.top_panel.table.handles{3,2}   = obj.top_panel.Viterbi_range_TE;
obj.top_panel.table.handles{4,1}   = obj.top_panel.reinterp_mode_label;
obj.top_panel.table.handles{4,2}   = obj.top_panel.reinterp_mode_cbox;
%% Column restriction
obj.top_panel.table.handles{5,1}   = obj.top_panel.column_restriction_label;
obj.top_panel.table.handles{5,2}   = obj.top_panel.column_restriction_cbox;
%% Top suppression
obj.top_panel.table.handles{6,1}   = obj.top_panel.top_sup_label;
obj.top_panel.table.handles{6,2}   = obj.top_panel.top_sup_cbox;
%% Multiple suppression
obj.top_panel.table.handles{7,1}   = obj.top_panel.mult_sup_label;
obj.top_panel.table.handles{7,2}   = obj.top_panel.mult_sup_cbox;
%% Quality output
obj.top_panel.table.handles{8,1}   = obj.top_panel.quality_output_label;
obj.top_panel.table.handles{8,2}   = obj.top_panel.quality_output_cbox;
%% Quality threshold
obj.top_panel.table.handles{9,1}   = obj.top_panel.quality_threshold_label;
obj.top_panel.table.handles{9,2}   = obj.top_panel.quality_threshold_TE;
%% Smoothness weight
obj.top_panel.table.handles{10,1}  = obj.top_panel.smoothness_weight_label;
obj.top_panel.table.handles{10,2}  = obj.top_panel.smoothness_weight_TE;
%% Smoothness variance
obj.top_panel.table.handles{11,1}  = obj.top_panel.smoothness_variance_label;
obj.top_panel.table.handles{11,2}  = obj.top_panel.smoothness_variance_TE;
%% Repulsion
obj.top_panel.table.handles{12,1}  = obj.top_panel.repulsion_label;
obj.top_panel.table.handles{12,2}  = obj.top_panel.repulsion_TE;
%% Ice bin threshold
obj.top_panel.table.handles{13,1}  = obj.top_panel.icebinthr_label;
obj.top_panel.table.handles{13,2}  = obj.top_panel.icebinthr_TE;

clear rows cols

% Draw table
table_draw(obj.top_panel.table);

if obj.first_time
  obj.first_time = false;
else
  set(obj.h_fig,'visible','on');
end

return;