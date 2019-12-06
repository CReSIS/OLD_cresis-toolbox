function cmds_execute(obj,cmds_list,cmds_direction)
% cmds_execute(obj,cmds_list,cmds_direction)
%
% Executes tool commands

if strcmpi(cmds_direction,'redo')
  %% Redo commands
  for cmd_idx = 1:length(cmds_list)
    for sub_idx = 1:length(cmds_list{cmd_idx})
      if strcmpi(cmds_list{cmd_idx}(sub_idx).redo_cmd,'insert')
        cmds_execute_insert(obj,cmds_list{cmd_idx}(sub_idx).redo_args);
      elseif strcmpi(cmds_list{cmd_idx}(sub_idx).redo_cmd,'delete')
        cmds_execute_delete(obj,cmds_list{cmd_idx}(sub_idx).redo_args);
      end
    end
  end
  
else
  %% Undo commands (cmds_direction == 'undo')
  for cmd_idx = 1:length(cmds_list)
    for sub_idx = 1:length(cmds_list{cmd_idx})
      if strcmpi(cmds_list{cmd_idx}(sub_idx).undo_cmd,'insert')
        cmds_execute_insert(obj,cmds_list{cmd_idx}(sub_idx).undo_args);
      elseif strcmpi(cmds_list{cmd_idx}(sub_idx).undo_cmd,'delete')
        cmds_execute_delete(obj,cmds_list{cmd_idx}(sub_idx).undo_args);
      end
    end
  end
  
end

obj.update_layer_plots();

end

function args = cmds_execute_insert(obj,args)
% Execute the insert command on obj.eg.layer
% args{1} = database layer IDs
% args{2} = point IDs
% args{3} = twtt
% args{4} = layer type
% args{5} = layer quality

%% Convert layer ID's to layer indices
layer_idx = find(args{1} == obj.eg.layer_id);

if isempty(layer_idx)
  % This echowin does not have this layer shown, so no updates required
  return;
end

%% Convert point ID's to point indices
point_idxs = [];
point_id_mask = logical(zeros(size(args{2})));
for point_id_idx = 1:length(args{2})
  point_id = args{2}(point_id_idx);
  new_point = find(point_id == obj.eg.map_id);
  if ~isempty(new_point)
    point_id_mask(point_id_idx) = true;
    point_idxs(end+1) = new_point;
  end
end
%% Apply insert
obj.eg.layer.y{layer_idx}(point_idxs) = args{3}(point_id_mask);
obj.eg.layer.type{layer_idx}(point_idxs) = args{4}(point_id_mask);
obj.eg.layer.qual{layer_idx}(point_idxs) = args{5}(point_id_mask);

%% Convert units from twtt to current units
yaxis_choice = get(obj.left_panel.yaxisPM,'Value');
if yaxis_choice == 1 % TWTT
  obj.eg.layer.y_curUnit{layer_idx}(point_idxs) = obj.eg.layer.y{layer_idx}(point_idxs) * 1e6;
  
elseif yaxis_choice == 2 % WGS_84 Elevation
  physical_constants;
  elevation = interp1(obj.eg.gps_time,...
    obj.eg.elevation,...
    obj.eg.layer.x{layer_idx},'linear');
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surface,...
    obj.eg.layer.x{layer_idx},'linear');
  for point_idx = point_idxs
    if isnan(obj.eg.layer.y{layer_idx}(point_idx))
      obj.eg.layer.y_curUnit{layer_idx}(point_idx) = NaN;
    else
      range = min(obj.eg.layer.y{layer_idx}(point_idx),surface(point_idx))*c/2 ...
        +max(0,obj.eg.layer.y{layer_idx}(point_idx)-surface(point_idx)) * c/(sqrt(er_ice)*2);
      obj.eg.layer.y_curUnit{layer_idx}(point_idx) = elevation(point_idx) - range;
    end
  end
  
elseif yaxis_choice == 3 % Depth/Range
  physical_constants;
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surface,...
    obj.eg.layer.x{layer_idx},'linear');
  for point_idx = point_idxs
    if isnan(obj.eg.layer.y{layer_idx}(point_idx))
      obj.eg.layer.y_curUnit{layer_idx}(point_idx) = NaN;
    else
      obj.eg.layer.y_curUnit{layer_idx}(point_idx) = min(obj.eg.layer.y{layer_idx}(point_idx),surface(point_idx))*c/2 ...
        +max(0,obj.eg.layer.y{layer_idx}(point_idx)-surface(point_idx)) * c/(sqrt(er_ice)*2);
    end
  end
  
elseif yaxis_choice == 4 % Range bin
  obj.eg.layer.y_curUnit{layer_idx}(point_idxs) = interp1(obj.eg.time,...
    1:length(obj.eg.time),...
    obj.eg.layer.y{layer_idx}(point_idxs),'linear');
end

end

function args = cmds_execute_delete(obj,args)
% Execute the delete command on obj.eg.layer and updates the layer handles
% obj.layer_h and obj.quality_h.
%
% args{1} = database layer IDs
% args{2} = [x_min x_max y_min y_max] with x in gps-time and y in twtt
%   units

% obj.eg.layer.x{layer_idx} % gps-time
% obj.eg.layer.y{layer_idx} % twtt
% obj.eg.layer.qual{layer_idx} % integer 1-3
% obj.eg.layer.type{layer_idx} % this is either 1 (manual) or 2 (auto)
% obj.eg.layer.x_curUnit{layer_idx} % current x-axis units (e.g. along-track km)
% obj.eg.layer.y_curUnit{layer_idx} % current y-axis units (e.g. WGS-84 m)


%% Convert layer ID's to layer indices
layer_idx = find(args{1} == obj.eg.layer_id);

if isempty(layer_idx)
  % This echowin does not have this layer shown, so no updates required
  return;
end

%% Determine which point indices need to be updated
point_idxs = find(obj.eg.layer.x{layer_idx} > args{2}(1) & obj.eg.layer.x{layer_idx} < args{2}(2) ...
  & obj.eg.layer.y{layer_idx} > args{2}(3) & obj.eg.layer.y{layer_idx} < args{2}(4));

obj.eg.layer.y{layer_idx}(point_idxs) = NaN;
obj.eg.layer.qual{layer_idx}(point_idxs) = 1;
obj.eg.layer.type{layer_idx}(point_idxs) = 1;
obj.eg.layer.y_curUnit{layer_idx}(point_idxs) = NaN;

end