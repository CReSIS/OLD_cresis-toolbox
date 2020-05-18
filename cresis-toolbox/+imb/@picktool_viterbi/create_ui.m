function create_ui(obj)
% create_ui(obj)
%
% Creates components for the window's UI when the echogram
% window is first opened.
%

set(obj.h_fig,'DockControls','off')
set(obj.h_fig,'NumberTitle','off');
if strcmpi(class(obj.h_fig),'double')
  set(obj.h_fig,'Name',sprintf('%d: %s tool parameters', obj.h_fig, obj.tool_name_title));
else
  set(obj.h_fig,'Name',sprintf('%d: %s tool parameters', obj.h_fig.Number, obj.tool_name_title));
end
set(obj.h_fig,'ToolBar','none');
set(obj.h_fig,'MenuBar','none'); 
set(obj.h_fig,'KeyPressFcn',@obj.key_press);
set(obj.h_fig,'CloseRequestFcn',@obj.close_win);

set(obj.h_fig,'visible','off');

% set default position (changed when window accessed)
set(obj.h_fig,'Units','Pixels');
set(obj.h_fig,'Position',[0 0 obj.w obj.h]);

%==========================================================================
% top panel
obj.top_panel.handle = uipanel('Parent',obj.h_fig);
set(obj.top_panel.handle,'HighlightColor',[0.8 0.8 0.8]);
set(obj.top_panel.handle,'ShadowColor',[0.6 0.6 0.6]);

%--------------------------------------
% table
obj.table.ui=obj.h_fig;

obj.table.handles{1,1}        = obj.top_panel.handle;
obj.table.width(1,1)          = inf;
obj.table.height(1,1)         = inf;
obj.table.width_margin(1,1)   = 0;
obj.table.height_margin(1,1)  = 0;

table_draw(obj.table);

%==========================================================================
% top panel table contents

%----Mode dropdown
tooltip = 'Switch Viterbi functionality';
obj.top_panel.tool_PM = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.tool_PM,'Style','popupmenu');
set(obj.top_panel.tool_PM,'String',{'basic'});
set(obj.top_panel.tool_PM,'Value',1)
set(obj.top_panel.tool_PM,'TooltipString', tooltip);
%-----mode label
obj.top_panel.mode_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.mode_label,'Style','text');
set(obj.top_panel.mode_label,'String','Mode');
set(obj.top_panel.mode_label,'TooltipString', tooltip);

%----insert range
tooltip = 'Viterbi will search +/- this many bins for the peak intensity on insert';
obj.top_panel.insert_range_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.insert_range_label,'Style','text');
set(obj.top_panel.insert_range_label,'String','Max point range:');
set(obj.top_panel.insert_range_label,'TooltipString', tooltip);
%----insert pt search range box
obj.top_panel.insert_range_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.insert_range_TE,'Style','edit');
set(obj.top_panel.insert_range_TE,'String','5');
set(obj.top_panel.insert_range_TE,'TooltipString', tooltip);

%----Horizontal Bounding
tooltip = 'How to bound the input and output of Viterbi horizontally.';
obj.top_panel.hori_bound_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.hori_bound_label,'Style','text');
set(obj.top_panel.hori_bound_label,'String',sprintf('Horizontal Bounding:'));
set(obj.top_panel.hori_bound_label,'TooltipString', tooltip);
%----Horizontal Bounding radio buttons
obj.top_panel.hori_bound_bg = uibuttongroup('Parent',obj.top_panel.handle);  % bg = button_group
obj.top_panel.r_echo = uicontrol('Parent', obj.top_panel.hori_bound_bg);
set(obj.top_panel.r_echo,'Style','radiobutton');
set(obj.top_panel.r_echo,'String','Entire Echogram');
set(obj.top_panel.r_echo,'Position',[0 30 200 15]);
set(obj.top_panel.r_echo,'Value', 0);
set(obj.top_panel.r_echo,'TooltipString', 'Pass in the entire echogram (no bounding).');
obj.top_panel.r_sel = uicontrol('Parent', obj.top_panel.hori_bound_bg);
set(obj.top_panel.r_sel,'Style','radiobutton');
set(obj.top_panel.r_sel,'String','Selection Box');
set(obj.top_panel.r_sel,'Position',[0 15 200 15]);
set(obj.top_panel.r_sel,'Value', 0);
set(obj.top_panel.r_sel,'TooltipString', 'Use echogram within horizontal bounds of selection box.');
obj.top_panel.r_extr = uicontrol('Parent', obj.top_panel.hori_bound_bg);
set(obj.top_panel.r_extr,'Style','radiobutton');
set(obj.top_panel.r_extr,'String','Extreme Groundtruth');
set(obj.top_panel.r_extr,'Position',[0 0 200 15]);
set(obj.top_panel.r_extr,'Value', 1);
set(obj.top_panel.r_extr,'TooltipString', 'Use echogram within outer-most (horizontally) ground truth points within selection box.');

%----Vertical Bounding
tooltip = 'How to bound the input and output of Viterbi vertically.';
obj.top_panel.vert_bound_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.vert_bound_label,'Style','text');
set(obj.top_panel.vert_bound_label,'String',sprintf('Vertical Bounding:'));
set(obj.top_panel.vert_bound_label,'TooltipString', tooltip);
%----Vertical Bounding radio buttons
obj.top_panel.vert_bound_bg = uibuttongroup('Parent',obj.top_panel.handle);  % bg = button_group
obj.top_panel.r_echo_vert = uicontrol('Parent', obj.top_panel.vert_bound_bg);
set(obj.top_panel.r_echo_vert,'Style','radiobutton');
set(obj.top_panel.r_echo_vert,'String','Entire Echogram');
set(obj.top_panel.r_echo_vert,'Position',[0 15 200 15]);
set(obj.top_panel.r_echo_vert,'Value', 0);
set(obj.top_panel.r_echo_vert,'TooltipString', 'Pass in the entire echogram (no vertical bounding).');
obj.top_panel.r_sel_vert = uicontrol('Parent', obj.top_panel.vert_bound_bg);
set(obj.top_panel.r_sel_vert,'Style','radiobutton');
set(obj.top_panel.r_sel_vert,'String','Selection Box');
set(obj.top_panel.r_sel_vert,'Position',[0 0 200 15]);
set(obj.top_panel.r_sel_vert,'Value', 1);
set(obj.top_panel.r_sel_vert,'TooltipString', 'Use echogram within vertical bounds of selection box.');

%----multiple weight label
tooltip = 'Amount by which to repel surface multiples if suppression enabled. Greater value = greater avoidance.';
obj.top_panel.mult_weight_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.mult_weight_label,'Style','text');
set(obj.top_panel.mult_weight_label,'String','Multiple Repulsion:');
set(obj.top_panel.mult_weight_label,'TooltipString', tooltip);
%----multiple weight box
obj.top_panel.mult_weight_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.mult_weight_TE,'Style','edit');
set(obj.top_panel.mult_weight_TE,'String', '100');
set(obj.top_panel.mult_weight_TE,'TooltipString', tooltip);

%----multiple weight decay label
tooltip = 'Multiply repulsion of each subsequent multiple by this amount to reduce suppression of faded multiples. Smaller = faster repulsion decay.';
obj.top_panel.mult_weight_decay_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.mult_weight_decay_label,'Style','text');
set(obj.top_panel.mult_weight_decay_label,'String','Multiple Decay:');
set(obj.top_panel.mult_weight_decay_label,'TooltipString', tooltip);
%----multiple weight decay box
obj.top_panel.mult_weight_decay_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.mult_weight_decay_TE,'Style','edit');
set(obj.top_panel.mult_weight_decay_TE,'String', '0');
set(obj.top_panel.mult_weight_decay_TE,'TooltipString', tooltip);

%----multiple weight local decay label
tooltip = 'Multiply the multiple suppression repulsion by this amount for every subsequent bin past the multiple. Smaller = faster repulsion decay.';
obj.top_panel.mult_weight_local_decay_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.mult_weight_local_decay_label,'Style','text');
set(obj.top_panel.mult_weight_local_decay_label,'String','Multiple Local Decay:');
set(obj.top_panel.mult_weight_local_decay_label,'TooltipString', tooltip);
%----multiple weight local decay box
obj.top_panel.mult_weight_local_decay_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.mult_weight_local_decay_TE,'Style','edit');
set(obj.top_panel.mult_weight_local_decay_TE,'String', '0.8');
set(obj.top_panel.mult_weight_local_decay_TE,'TooltipString', tooltip);

%----along track weight label
tooltip = 'The weight by which to multiply the binary cost. Greater weight = smoother';
obj.top_panel.along_track_weight_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.along_track_weight_label,'Style','text');
set(obj.top_panel.along_track_weight_label,'String','Along Track weight:');
set(obj.top_panel.along_track_weight_label,'TooltipString', tooltip);
%----along track weight box
obj.top_panel.along_track_weight_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.along_track_weight_TE,'Style','edit');
set(obj.top_panel.along_track_weight_TE,'String', '1');
set(obj.top_panel.along_track_weight_TE,'TooltipString', tooltip);

%----image magnitude weight label
tooltip = 'The weight by which to multiply the image magnitude cost. Greater weight = prefer greater image magnitude';
obj.top_panel.image_mag_weight_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.image_mag_weight_label,'Style','text');
set(obj.top_panel.image_mag_weight_label,'String','Image Weight:');
set(obj.top_panel.image_mag_weight_label,'TooltipString', tooltip);
%----image magnitude weight box
obj.top_panel.image_mag_weight_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.image_mag_weight_TE,'Style','edit');
set(obj.top_panel.image_mag_weight_TE,'String', '1');
set(obj.top_panel.image_mag_weight_TE,'TooltipString', tooltip);

%----gt weight label
tooltip = 'The weight by which to multiply the ground truth cost. Greater weight = prefer closer to ground truth';
obj.top_panel.ground_truth_weight_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.ground_truth_weight_label,'Style','text');
set(obj.top_panel.ground_truth_weight_label,'String','Ground Truth Weight:');
set(obj.top_panel.ground_truth_weight_label,'TooltipString', tooltip);
%----gt weight box
obj.top_panel.ground_truth_weight_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.ground_truth_weight_TE,'Style','edit');
set(obj.top_panel.ground_truth_weight_TE,'String', '1');
set(obj.top_panel.ground_truth_weight_TE,'TooltipString', tooltip);

%----gt cutoff label
tooltip = 'Points must be chosen within this many rangebins of a ground truth point when present. -1 for any distance allowed.';
obj.top_panel.ground_truth_cutoff_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.ground_truth_cutoff_label,'Style','text');
set(obj.top_panel.ground_truth_cutoff_label,'String','Ground Truth Cutoff:');
set(obj.top_panel.ground_truth_cutoff_label,'TooltipString', tooltip);
%----gt cutoff box
obj.top_panel.ground_truth_cutoff_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.ground_truth_cutoff_TE,'Style','edit');
set(obj.top_panel.ground_truth_cutoff_TE,'String', '5');
set(obj.top_panel.ground_truth_cutoff_TE,'TooltipString', tooltip);
%%
%---------------------------------------------------------------------------------------------
% set up top panel table
cols = 2;
rows = 14;  % Just keep this larger or equal to actual number of rows.

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

row = 0;

%% Insert Range
row = row + 1;
obj.top_panel.table.handles{row,1}   = obj.top_panel.insert_range_label;
obj.top_panel.table.handles{row,2}   = obj.top_panel.insert_range_TE;
%% Horizontal Bound
row = row + 1;
obj.top_panel.table.handles{row,1}   = obj.top_panel.hori_bound_label;
obj.top_panel.table.handles{row,2}   = obj.top_panel.hori_bound_bg;
obj.top_panel.table.height(row, :)   = 40;
%% Vertical Bound
row = row + 1;
obj.top_panel.table.handles{row,1}   = obj.top_panel.vert_bound_label;
obj.top_panel.table.handles{row,2}   = obj.top_panel.vert_bound_bg;
obj.top_panel.table.height(row, :)   = 30;
%% Multiple Weight
row = row + 1;
obj.top_panel.table.handles{row,1}   = obj.top_panel.mult_weight_label;
obj.top_panel.table.handles{row,2}   = obj.top_panel.mult_weight_TE;
%% Multiple Weight Decay
row = row + 1;
obj.top_panel.table.handles{row,1}   = obj.top_panel.mult_weight_decay_label;
obj.top_panel.table.handles{row,2}   = obj.top_panel.mult_weight_decay_TE;
%% Multiple Weight Local Decay 
row = row + 1;
obj.top_panel.table.handles{row,1}   = obj.top_panel.mult_weight_local_decay_label;
obj.top_panel.table.handles{row,2}   = obj.top_panel.mult_weight_local_decay_TE;
%% Along-track Weight
row = row + 1;
obj.top_panel.table.handles{row,1}  = obj.top_panel.along_track_weight_label;
obj.top_panel.table.handles{row,2}  = obj.top_panel.along_track_weight_TE;
%% Image magnitude weight
row = row + 1;
obj.top_panel.table.handles{row,1}  = obj.top_panel.image_mag_weight_label;
obj.top_panel.table.handles{row,2}  = obj.top_panel.image_mag_weight_TE;
%% gt weight
row = row + 1;
obj.top_panel.table.handles{row,1}  = obj.top_panel.ground_truth_weight_label;
obj.top_panel.table.handles{row,2}  = obj.top_panel.ground_truth_weight_TE;
%% gt cutoff
row = row + 1;
obj.top_panel.table.handles{row,1}  = obj.top_panel.ground_truth_cutoff_label;
obj.top_panel.table.handles{row,2}  = obj.top_panel.ground_truth_cutoff_TE;

if row > rows
  warning('Viterbi create_ui does not have default values for new rows. Update rows variable to match number of rows present.');
end

clear row cols

% Draw table
table_draw(obj.top_panel.table);
