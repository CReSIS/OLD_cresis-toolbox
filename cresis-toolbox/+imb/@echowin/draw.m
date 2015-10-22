function draw(obj,param)
% draw(obj,param)
%
% Member function of imb.echowin class.
%
% Loads a new dataset into the pick window from the map window.

obj.busy_mode = true;
set(obj.h_fig,'Pointer','watch');
obj.status_text_set(sprintf('(%s) Drawing...', datestr(now,'HH:MM:SS')),'replace');
drawnow;
fprintf('START ECHOWIN DRAW (%s)\n',datestr(now,'HH:MM:SS'));

obj.eg.sources = param.sources;
obj.eg.layers = param.layers;
obj.eg.cur_sel = param.cur_sel;
obj.eg.system = param.system;
% get the projection matrix for the current location
obj.eg.projmat = imb.get_proj_info(obj.eg.cur_sel.location);

%% Drop all currently loaded data
obj.eg.data = [];
obj.eg.gps_time = [];
obj.eg.latitude = [];
obj.eg.longitude = [];
obj.eg.elevation = [];
obj.eg.surface = [];
obj.eg.frame_idxs = [];

%% Get all the frames for this segment
ops_param = struct('properties',[]);
ops_param.properties.location = obj.eg.cur_sel.location;
ops_param.properties.segment_id = obj.eg.cur_sel.segment_id;
[status,data] = opsGetSegmentInfo(obj.eg.system,ops_param);
obj.eg.frame_names = {};
obj.eg.start_gps_time = [];
obj.eg.stop_gps_time = [];
for idx = 1:length(data.properties.frame)
  obj.eg.frame_names{idx} = data.properties.frame{idx};
  obj.eg.start_gps_time(idx) = double(data.properties.start_gps_time(idx));
  obj.eg.stop_gps_time(idx) = double(data.properties.stop_gps_time(idx));
end
% Sort the list of frames in case database isn't sorted
[obj.eg.frame_names sorted_idxs] = sort(obj.eg.frame_names);
obj.eg.start_gps_time = obj.eg.start_gps_time(sorted_idxs);
obj.eg.stop_gps_time = obj.eg.stop_gps_time(sorted_idxs);
obj.eg.cur_sel.day_seg = data.properties.segment;

%% Update file listing of all source files
obj.update_source_fns_existence();

%% Load echogram data

set(obj.left_panel.sourceLB,'Value',1);
set(obj.left_panel.sourceLB,'String',obj.eg.sources);
%set(obj.left_panel.sourceLB,'Value',find(strcmp(obj.eg.sources,existing_sources(first_idx))));

frame_idxs = find(strcmp(obj.eg.cur_sel.frame_name,obj.eg.frame_names));
% Set current listbox frame to this frame since it is the one that was
% specifically asked for (load_echogram uses this value)
set(obj.left_panel.frameLB,'Value',frame_idxs);
% But load some other frames around this one for more responsive GUI behavior
desire_frame_idxs = frame_idxs(1) ...
  : min(frame_idxs(1)+obj.default_params.max_frames-1,length(obj.eg.frame_names));
clipped = 3;
load_new_data = 1;
x_min = -inf;
x_max = inf;
y_min = -inf;
y_max = inf;

% Update "Frames" listbox
set(obj.left_panel.frameLB,'String',obj.eg.frame_names);

[x_min,x_max,y_min,y_max] = obj.load_echogram(desire_frame_idxs,clipped,x_min,x_max,y_min,y_max);

%% Create variables for the layer list box
% ===================================================================
obj.eg.layer_id = obj.eg.layers.lyr_id;
obj.eg.layer.layer_names = obj.eg.layers.lyr_name;

obj.eg.layers.surface = NaN;
layerLB_param.bottom = NaN;
% find indices of surface and bottom layers
for idx = 1:length(obj.eg.layers.lyr_name)
  if strcmp(obj.eg.layers.lyr_name{idx},'surface')
    obj.eg.layers.surface = idx;
  elseif strcmp(obj.eg.layers.lyr_name{idx},'bottom')
    layerLB_param.bottom = idx;
  end
end
layerLB_param.surface = obj.eg.layers.surface;

layerLB_param.names = cell(1,length(obj.eg.layer_id));
for idx = 1:length(obj.eg.layer_id)
  layerLB_param.names{idx} = sprintf('(%d) %s',idx,obj.eg.layers.lyr_name{idx});
end

obj.layerLB_setdata(layerLB_param);

%% Load flightlines, layers, crossovers and place the cursor
% ===================================================================
obj.load_flightline();
obj.load_layers_init();
obj.load_layers();
obj.plot_echogram(x_min,x_max,y_min,y_max);
obj.plot_layers();
obj.crossovers_en = obj.eg.crossovers.gui.crossovers_en();
obj.load_crossovers();
obj.plot_cursors();

%% Set layer and cross over visibility
obj.set_visibility();

%% Plot cursor (if valid)
xlimits = xlim(obj.right_panel.axes.handle);
start_gps = obj.eg.image_gps_time(find(obj.eg.image_xaxis>=xlimits(1),1));
stop_gps = obj.eg.image_gps_time(find(obj.eg.image_xaxis<=xlimits(2),1,'last'));

if ~isempty(obj.cursor.gps_time) ...
    && obj.cursor.gps_time >= start_gps && obj.cursor.gps_time <= stop_gps
  obj.plot_cursors();
  str = obj.status_text_cursor();
  obj.status_text_set(str,'replace');
end

%% Plot new selection on flight path
notify(obj,'update_echowin_flightline');

% Set echowin figure callback properties (these are set here rather than
% create_ui since most of the do not work before data is loaded for the
% first time)
set(obj.h_fig,'WindowButtonDownFcn',@obj.button_down);
set(obj.h_fig,'WindowButtonMotionFcn',@obj.button_motion);
set(obj.h_fig,'WindowButtonUpFcn',@obj.button_up);
set(obj.h_fig,'WindowKeyPressFcn',@obj.key_press);
set(obj.h_fig,'WindowKeyReleaseFcn',@obj.key_release);
set(obj.h_fig,'WindowScrollWheelFcn',@obj.button_scroll);

obj.busy_mode = false;
if obj.zoom_mode
  set(obj.h_fig,'Pointer','custom');
else
  set(obj.h_fig,'Pointer','Arrow');
end

obj.status_text_set(sprintf(' done. (%s)', datestr(now,'HH:MM:SS')),'append');  

fprintf(' DONE (%s)\n',datestr(now,'HH:MM:SS'));

return;




