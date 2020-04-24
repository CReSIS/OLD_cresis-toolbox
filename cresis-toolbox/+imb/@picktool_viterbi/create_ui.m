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

%----column restriction label
tooltip = 'Crop echogram input horizontally to values between extreme ground truth points';
obj.top_panel.column_restriction_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.column_restriction_label,'Style','text');
set(obj.top_panel.column_restriction_label,'String','Column restriction:');
set(obj.top_panel.column_restriction_label,'TooltipString', tooltip);
%----column restriction cbox
obj.top_panel.column_restriction_cbox = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.column_restriction_cbox,'Style','checkbox');
set(obj.top_panel.column_restriction_cbox,'Value', 1);
set(obj.top_panel.column_restriction_cbox,'TooltipString', tooltip);

%----layers label
tooltip = 'List of layers to repulse or attract the viterbi layer. Enter as a vector. The first entry is the top and the viterbi layer may not exceed.';
obj.top_panel.layers_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.layers_label,'Style','text');
set(obj.top_panel.layers_label,'String','Layers:');
set(obj.top_panel.layers_label,'TooltipString', tooltip);
%----layers box
obj.top_panel.layers_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.layers_TE,'Style','edit');
set(obj.top_panel.layers_TE,'String', '[1]');
set(obj.top_panel.layers_TE,'TooltipString', tooltip);

%----layer weight label
tooltip = 'List of layer weights. Larger positive values cause repulsion. Larger negative values cause attraction.';
obj.top_panel.layers_weight_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.layers_weight_label,'Style','text');
set(obj.top_panel.layers_weight_label,'String','Layer Weights:');
set(obj.top_panel.layers_weight_label,'TooltipString', tooltip);
%----layer  weight box
obj.top_panel.layers_weight_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.layers_weight_TE,'Style','edit');
set(obj.top_panel.layers_weight_TE,'String', '[1000]');
set(obj.top_panel.layers_weight_TE,'TooltipString', tooltip);

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

%----surface slope label
tooltip = 'Use the slope of the surface layer as the expected slope of the target layer';
obj.top_panel.surf_slope_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.surf_slope_label,'Style','text');
set(obj.top_panel.surf_slope_label,'String',sprintf('Calc slope from surf:'));
set(obj.top_panel.surf_slope_label,'TooltipString', tooltip);
%----surface slope cbox
obj.top_panel.surf_slope_cbox = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.surf_slope_cbox,'Style','checkbox');
set(obj.top_panel.surf_slope_cbox,'Value', 1);
set(obj.top_panel.surf_slope_cbox,'TooltipString', tooltip);

%----max slope label
tooltip = 'The maximum allowed slope of the target layer. -1 for no max.';
obj.top_panel.max_slope_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.max_slope_label,'Style','text');
set(obj.top_panel.max_slope_label,'String','Max Slope:');
set(obj.top_panel.max_slope_label,'TooltipString', tooltip);
%----max slope box
obj.top_panel.max_slope_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.max_slope_TE,'Style','edit');
set(obj.top_panel.max_slope_TE,'String', '-1');
set(obj.top_panel.max_slope_TE,'TooltipString', tooltip);

%----transition weight label
tooltip = 'The weight by which to multiply the binary cost. Greater weight = smoother';
obj.top_panel.transition_weight_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.transition_weight_label,'Style','text');
set(obj.top_panel.transition_weight_label,'String','Transition weight:');
set(obj.top_panel.transition_weight_label,'TooltipString', tooltip);
%----transition weight box
obj.top_panel.transition_weight_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.transition_weight_TE,'Style','edit');
set(obj.top_panel.transition_weight_TE,'String', '1');
set(obj.top_panel.transition_weight_TE,'TooltipString', tooltip);

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
cols = 2;
row = 0;

%% Mode
row = row + 1;
obj.top_panel.table.handles{row,1}   = obj.top_panel.mode_label;
obj.top_panel.table.handles{row,2}   = obj.top_panel.tool_PM;
%% Insert Range
row = row + 1;
obj.top_panel.table.handles{row,1}   = obj.top_panel.insert_range_label;
obj.top_panel.table.handles{row,2}   = obj.top_panel.insert_range_TE;
%% Column restriction
row = row + 1;
obj.top_panel.table.handles{row,1}   = obj.top_panel.column_restriction_label;
obj.top_panel.table.handles{row,2}   = obj.top_panel.column_restriction_cbox;
%% Layers
row = row + 1;
obj.top_panel.table.handles{row,1}   = obj.top_panel.layers_label;
obj.top_panel.table.handles{row,2}   = obj.top_panel.layers_TE;
%% Layers Weight
row = row + 1;
obj.top_panel.table.handles{row,1}   = obj.top_panel.layers_weight_label;
obj.top_panel.table.handles{row,2}   = obj.top_panel.layers_weight_TE;
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
%% Transition Slope from Surface
row = row + 1;
obj.top_panel.table.handles{row,1}   = obj.top_panel.surf_slope_label;
obj.top_panel.table.handles{row,2}   = obj.top_panel.surf_slope_cbox;
%% Max Slope 
row = row + 1;
obj.top_panel.table.handles{row,1}  = obj.top_panel.max_slope_label;
obj.top_panel.table.handles{row,2}  = obj.top_panel.max_slope_TE;
%% Transition Weight
row = row + 1;
obj.top_panel.table.handles{row,1}  = obj.top_panel.transition_weight_label;
obj.top_panel.table.handles{row,2}  = obj.top_panel.transition_weight_TE;
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

% set up top panel table
default_dimensions = NaN*zeros(row,cols);
obj.top_panel.table.ui=obj.top_panel.handle;
obj.top_panel.table.width_margin = default_dimensions;
obj.top_panel.table.height_margin = default_dimensions;
obj.top_panel.table.false_width = default_dimensions;
obj.top_panel.table.false_height = default_dimensions;
obj.top_panel.table.offset = [0 0];

obj.top_panel.table.width = ones(row, cols) * inf;
obj.top_panel.table.height = ones(row, cols) * inf;
obj.top_panel.table.width_margin = ones(row, cols) * 1.5;
obj.top_panel.table.height_margin = ones(row, cols) * 1.5;

clear row cols

% Draw table
table_draw(obj.top_panel.table);
