% script tomo.run_surfData_modify.m
%
% Example script for running surfData_modify.m. Demonstrates a few of the
% most common operations to be performed with surfData_modify.
%
% Authors: John Paden

% =====================================================================
%% User Settings
% =====================================================================

% params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'));
params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'20140401_03','post');
params.cmd.generic = 1;
params.cmd.frms = [37 39 43:48];

surfdata_source = 'surfData';

SURFDATA_MODIFY_EXAMPLE = 'add_new_layer';
if strcmpi(SURFDATA_MODIFY_EXAMPLE,'update_fields')
  %% Example for updating fields of a particular layer
  layers = [1];
  args = [];
  args{1} = 'control_layer'
  args{2} = 7;
  
elseif strcmpi(SURFDATA_MODIFY_EXAMPLE,'add_new_layer')
  %% Example for adding a new layer in
  echo_source = 'CSA_music';
  ops_layer_name = 'surface';
  new_layer = 7;
  args = [];
  args{end+1} = 'plot_name_values';
  args{end+1} = {'color','green','marker','^'};
  args{end+1} = 'name';
  args{end+1} = 'surface gt';
  args{end+1} = 'surf_layer';
  args{end+1} = [];
  args{end+1} = 'active_layer';
  args{end+1} = 1;
  args{end+1} = 'mask_layer';
  args{end+1} = [];
  args{end+1} = 'control_layer';
  args{end+1} = 7;
  args{end+1} = 'visible';
  args{end+1} = true;
end

% =====================================================================
%% Automated Section
% =====================================================================

global gRadar;

if strcmpi(SURFDATA_MODIFY_EXAMPLE,'update_fields')
  %% Load each of the day segments
  for param_idx = 1:length(params)
    param = params(param_idx);
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) ...
        || ischar(param.cmd.generic) || ~param.cmd.generic
      continue;
    end
    
    param = merge_structs(param,gRadar);
    
    fprintf('surfData_modify %s\n', param.day_seg);
    tomo.surfData_modify(param,surfdata_source,layers,args{:});
  end
  
elseif strcmpi(SURFDATA_MODIFY_EXAMPLE,'add_new_layer')
  %% Load each of the day segments
  for param_idx = 1:length(params)
    param = params(param_idx);
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) ...
        || ischar(param.cmd.generic) || ~param.cmd.generic
      continue;
    end
    
    param = merge_structs(param,gRadar);
    
    % Determine which frames to process
    load(ct_filename_support(param,param.records.frames_fn,'frames'));
    
    if isempty(param.cmd.frms)
      param.cmd.frms = 1:length(frames.frame_idxs);
    end
    % Remove frames that do not exist from param.cmd.frms list
    [valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
    if length(valid_frms) ~= length(param.cmd.frms)
      bad_mask = ones(size(param.cmd.frms));
      bad_mask(keep_idxs) = 0;
      warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
        param.cmd.frms(find(bad_mask,1)));
      param.cmd.frms = valid_frms;
    end
    
    %% Only do one frame at a time
    all_frms = param.cmd.frms;
    for frm = all_frms(:).'
      param.cmd.frms = frm;
      echo_fn_dir = ct_filename_out(param,echo_source);
      echo_fn = fullfile(echo_fn_dir,sprintf('Data_%s_%03d.mat',param.day_seg, frm));
      
      mdata = load(echo_fn,'Surface','Time','GPS_time','Latitude','Longitude','Elevation','twtt');
      
      % Query OPS for surface and bottom information
      param_load_layers = param;
      param_load_layers.cmd.frms = round([-1,0,1] + frm);
      
      layer_params = [];
      idx = 0;
      idx = idx + 1;
      layer_params(idx).name = layer_name;
      layer_params(idx).source = 'ops';
      layers = opsLoadLayers(param_load_layers,layer_params);
      
      % Interpolate surface and bottom information to mdata
      master = [];
      master.GPS_time = mdata.GPS_time;
      master.Latitude = mdata.Latitude;
      master.Longitude = mdata.Longitude;
      master.Elevation = mdata.Elevation;
      for lay_idx = 1:length(layer_params)
        ops_layer = [];
        ops_layer{1}.gps_time = layers(lay_idx).gps_time;
        
        ops_layer{1}.type = layers(lay_idx).type;
        ops_layer{1}.quality = layers(lay_idx).quality;
        ops_layer{1}.twtt = layers(lay_idx).twtt;
        ops_layer{1}.type(isnan(ops_layer{1}.type)) = 2;
        ops_layer{1}.quality(isnan(ops_layer{1}.quality)) = 1;
        lay = opsInterpLayersToMasterGPSTime(master,ops_layer,[300 60]);
        layers(lay_idx).twtt_ref = lay.layerData{1}.value{2}.data;
      end
      
      % Add surface information to surfData file
      Surface = layers(1).twtt_ref;
      args{end+1} = 'x';
      args{end+1} = repmat((1:size(mdata.twtt,1)).',[1 size(mdata.twtt,2)]);
      
      args{end+1} = 'y';
      surf_layer = NaN*zeros(size(mdata.twtt));
      surf_layer(floor(size(mdata.twtt,1)/2)+1,:) = interp1(mdata.Time, 1:length(mdata.Time), Surface);
      args{end+1} = surf_layer;
      
      fprintf('surfData_modify %s_%03d\n', param.day_seg, frm);
      tomo.surfData_modify(param,surfdata_source,new_layer,args{:});
    end
  end
  
end
