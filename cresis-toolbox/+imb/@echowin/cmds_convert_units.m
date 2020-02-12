function cmds = cmds_convert_units(obj,cmds)
% cmds = convert_cmd_units(obj,cmds)
%
% Converts tool commands from current units to gps-time and twtt

for cmd_idx = 1:length(cmds)
  if strcmpi(cmds(cmd_idx).undo_cmd,'insert')
    cmds(cmd_idx).undo_args = cmds_convert_units_insert(obj,cmds(cmd_idx).undo_args);
  elseif strcmpi(cmds(cmd_idx).undo_cmd,'delete')
    cmds(cmd_idx).undo_args = cmds_convert_units_delete(obj,cmds(cmd_idx).undo_args);
  end
  if strcmpi(cmds(cmd_idx).redo_cmd,'insert')
    cmds(cmd_idx).redo_args = cmds_convert_units_insert(obj,cmds(cmd_idx).redo_args);
  elseif strcmpi(cmds(cmd_idx).redo_cmd,'delete')
    cmds(cmd_idx).redo_args = cmds_convert_units_delete(obj,cmds(cmd_idx).redo_args);
  end
end

end

function args = cmds_convert_units_insert(obj,args)
% args = cmds_convert_units_insert(obj,args)
%
% Convert units of the insert command arguments
% args{1} = layer indices --> database layer IDs
% args{2} = point indices --> point IDs
% args{3} = layer image y-units --> current image units converted to twtt
% args{4} = layer type (no conversion required)
% args{5} = layer quality (no conversion required)

physical_constants;
vel_air = c/2;
slowness_air = 1/vel_air;
vel_ice = c/(2*sqrt(er_ice));
slowness_ice = 1/vel_ice;

% Convert y-axis range units
yaxis_choice = get(obj.left_panel.yaxisPM,'Value');

if yaxis_choice == 1 % TWTT
  args{3} = args{3}/(1e6);
  
elseif yaxis_choice == 2 % WGS_84 Elevation
  elevation = interp1(obj.eg.gps_time,obj.eg.elev,obj.eg.map_gps_time(args{2}),'linear','extrap');
  surface = interp1(obj.eg.gps_time,obj.eg.surf_twtt,obj.eg.map_gps_time(args{2}),'linear','extrap');
  
  for idx = 1:length(args{3})
    if (elevation(idx) - args{3}(idx)) < surface(idx)*vel_air % above surface
      args{3}(idx) = (elevation(idx)-args{3}(idx))*slowness_air;
    else
      args{3}(idx) = surface(idx) + (elevation(idx) - args{3}(idx) - surface(idx)*vel_air)*slowness_ice;
    end
  end
  
elseif yaxis_choice == 3 % Depth/Range
  surface = interp1(obj.eg.gps_time,obj.eg.surf_twtt,obj.eg.map_gps_time(args{2}),'linear','extrap');
  for idx = 1:length(args{3})
    if args{3}(idx) < surface(idx)*vel_air % above surface
      args{3}(idx) = args{3}(idx)*slowness_air;
    else
      args{3}(idx) = surface(idx) + (args{3}(idx) - surface(idx)*vel_air)*slowness_ice;
    end
  end
  
elseif yaxis_choice == 4 % Range bin
  args{3} = interp1(obj.eg.image_yaxis,obj.eg.time,...
    args{3},'linear');
  
elseif yaxis_choice == 5 % Surface flat
  surface = interp1(obj.eg.gps_time,obj.eg.surf_twtt,obj.eg.map_gps_time(args{2}),'linear','extrap');
  
  below_ice_mask = args{3} > 0;
  args{3}(below_ice_mask) = surface(below_ice_mask)+args{3}(below_ice_mask)*slowness_ice; % Below ice
  args{3}(~below_ice_mask) = surface(~below_ice_mask)+args{3}(~below_ice_mask)*slowness_air; % Above ice
  
end

% Change layer idxs for layer ids
args{1} = obj.eg.layers.lyr_id(args{1});
% Change point idxs for point ids
args{2} = obj.eg.map_id(args{2});

end

function args = cmds_convert_units_delete(obj,args)
% args = cmds_convert_units_delete(obj,args)
%
% Convert units of the delete command arguments
% args{1} = layer indices --> database layer IDs
% args{2} = [x_min x_max y_min y_max] in image units --> gps-time and twtt
%   units
% args{3} = first and last point idxs --> point IDs

physical_constants;
vel_air = c/2;
slowness_air = 1/vel_air;
vel_ice = c/(2*sqrt(er_ice));
slowness_ice = 1/vel_ice;

% Convert y-axis range units
yaxis_choice = get(obj.left_panel.yaxisPM,'Value');

if yaxis_choice == 1 % TWTT
  args{2}(3:4) = args{2}(3:4)/(1e6);
  
elseif yaxis_choice == 2 % WGS_84 Elevation
  elevation = interp1(obj.eg.gps_time,obj.eg.elev,obj.eg.map_gps_time(args{3}(1)),'linear','extrap');
  surface = interp1(obj.eg.gps_time,obj.eg.surf_twtt,obj.eg.map_gps_time(args{3}(1)),'linear','extrap');
  
  % We are given a rectangle to delete in WGS-84 coordinates. When we
  % translate this to twtt, in the general case we will end up with a
  % parallelogram. We approximate the parallelogram with a rectangle
  % by looking at just the left side elevation/surface (i.e. we ignore
  % the right side coordinates by just using "args{3}(1)" and not
  % "args{3}(2)")
  for idx = 3:4
    if (elevation - args{2}(idx)) < surface*vel_air % above surface
      args{2}(idx) = (elevation-args{2}(idx))*slowness_air;
    else
      args{2}(idx) = surface + (elevation - args{2}(idx) - surface*vel_air)*slowness_ice;
    end
  end
  % Resort the range because WGS-84 elevation is opposite sign of twtt
  args{2}(3:4) = sort(args{2}(3:4));
  
elseif yaxis_choice == 3 % Range
  surface = interp1(obj.eg.gps_time,obj.eg.surf_twtt,obj.eg.map_gps_time(args{3}(1)),'linear','extrap');
  
  % We are given a rectangle to delete in Range coordinates. When we
  % translate this to twtt, in the general case we will end up with a
  % parallelogram. We approximate the parallelogram with a rectangle
  % by looking at just the left side elevation/surface.
  for idx = 3:4
    if args{2}(idx) < surface*vel_air % above surface
      args{2}(idx) = args{2}(idx)*slowness_air;
    else
      args{2}(idx) = surface + (args{2}(idx) - surface*vel_air)*slowness_ice;
    end
  end
  
elseif yaxis_choice == 4 % Range bin
  args{2}(3:4) = interp1(obj.eg.image_yaxis,obj.eg.time,...
    args{2}(3:4),'linear');
  
elseif yaxis_choice == 5 % Surface flat
  surface = interp1(obj.eg.gps_time,obj.eg.surf_twtt,obj.eg.map_gps_time(args{3}(1)),'linear','extrap');
  
  below_ice_mask = args{2} > 0;
  % Only modify indices 3,4:
  args{2}(below_ice_mask & [0 0 1 1]) = surface+args{2}(below_ice_mask & [0 0 1 1])*slowness_ice; % Below ice
  args{2}(~below_ice_mask & [0 0 1 1]) = surface+args{2}(~below_ice_mask & [0 0 1 1])*slowness_air; % Above ice
end


% Change layer idxs for layer ids
args{1} = obj.eg.layers.lyr_id(args{1});

% Convert x-axis range units
args{2}(1:2) = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,args{2}(1:2),'linear','extrap');

% Change point idxs for point ids
args{3} = obj.eg.map_id(args{3});
end
