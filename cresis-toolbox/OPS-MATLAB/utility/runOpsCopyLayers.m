% script runOpsCopyLayers.m
%
% Example script for running opsCopyLayers.m. Demonstrates a few of the
% most common operations to be performed with opsCopyLayers.
%
% Authors: John Paden

% =====================================================================
%% User Settings
% =====================================================================

% Set the parameter spreadsheet 
if 1
  % Use param spreadsheet:
  params = read_param_xls(ct_filename_param('rds_param_2010_Greenland_DC8.xls'),'','post');
else
  % Use param spreadsheet, but override settings for specific segments
  params = read_param_xls(ct_filename_param('rds_param_2010_Greenland_DC8.xls'),'20100414_02','post');
  for param_idx = 1:length(params)
    params(param_idx).cmd.generic = 1;
    params(param_idx).cmd.frms = [];
  end
end

% Set the operation to run (just choose one operation)
if 0
  runOpsCopyLayers_operation = 'copy_layer_nonmatch_sys';
else
  runOpsCopyLayers_operation = 'copy_layer';
end

%% copy_layer: Copy layer from one location to another and optionally apply operation during copy
if strcmp(runOpsCopyLayers_operation,'copy_layer')
  copy_param = [];
  copy_param.layer_source.existence_check = false;
  copy_param.layer_dest.existence_check = false;

  % Set the layer name for the source (e.g. 'surface', 'bottom')
  copy_param.layer_source.name = 'surface';
  
  % Set the layer name for the destination (e.g. 'surface', 'bottom')
  copy_param.layer_dest.name = 'surface';

  % Set the source (choose one)
  if 1
    copy_param.layer_source.source = 'ops';
  elseif 0
    copy_param.layer_source.source = 'records';
  elseif 0
    copy_param.layer_source.source = 'echogram';
    % Set the echogram source if using echogram
    if 1
      copy_param.layer_source.echogram_source = 'qlook';
    elseif 0
      copy_param.layer_source.echogram_source = 'deconv';
    else
      copy_param.layer_source.echogram_source = 'standard';
    end
  elseif 0
    copy_param.layer_source.source = 'layerdata';
    copy_param.layer_source.layerdata_source = 'layerData';
  else
    copy_param.layer_source.source = 'lidar';
    copy_param.layer_source.lidar_source = 'awi';
  end

  if 1
    copy_param.copy_method = 'overwrite';
  elseif 0
    copy_param.copy_method = 'fillgaps';
  else
    copy_param.copy_method = 'merge';
  end
  
  if strcmpi(copy_param.layer_dest.name,'surface')
    copy_param.gaps_fill.method = 'interp_finite';
  else
    copy_param.gaps_fill.method = 'preserve_gaps';
    copy_param.gaps_fill.method_args = [40 20];
  end
  
  % Set the twtt offset (for positive offset layer shifts down)
  twtt_offset = 0;
  
  % Set the GPS time offset (for positive offset layer shifts right)
  gps_time_offset = 0;
  
  if twtt_offset ~= 0 || gps_time_offset ~= 0
    warning('You have set a nonzero twtt_offset(%.12g) or gps_time_offset(%.3g). Normally these are both zero. Please verify that this is what you want to do before running "dbcont" to continue.\n', twtt_offset, gps_time_offset);
    keyboard
    copy_param.eval.cmd = sprintf('source = interp1(gps_time+%.3g,source + %.12g,gps_time);',gps_time_offset,twtt_offset);
  end
  
  % Set overwrite quality level (e.g. []: do not overwrite, 1: good, 2: medium, 3: bad)
  quality = [];
  if ~isempty(quality) && any(quality == [1 2 3])
    warning('You have set quality to %d. Normally it should be []. Please verify that you want to overwrite the quality level before running "dbcont" to continue.\n', quality);
    keyboard
    copy_param.quality.mode = 'overwrite';
    copy_param.quality.quality = quality;
  end
  
  % Set the destination (choose one): it can be the same as the source
  if 1
    copy_param.layer_dest.source = 'ops';
  elseif 0
    copy_param.layer_dest.source = 'records';
  elseif 0
    copy_param.layer_dest.source = 'echogram';
    % Set the echogram source if using echogram
    if 1
      copy_param.layer_dest.echogram_source = 'qlook';
    elseif 0
      copy_param.layer_dest.echogram_source = 'deconv';
    else
      copy_param.layer_dest.echogram_source = 'standard';
    end
  elseif 0
    copy_param.layer_dest.source = 'layerdata';
    copy_param.layer_dest.layerdata_source = 'layerData';
  end
  
end
%% Copy surface from mcords records to snow layerdata  
if strcmp(runOpsCopyLayers_operation,'copy_layer_nonmatch_sys')
  
  % Load mcords records data
  load_params = read_param_xls(ct_filename_param('rds_param_2015_Greenland_Polar6.xls'),'20150913.*','post');  
  layer_params = []; idx = 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'records';
  global gRadar;
  gps_time = []; twtt = [];
  for param_idx = 1:length(load_params)
    param = load_params(param_idx);
    param = merge_structs(param,gRadar);
    layers = opsLoadLayers(param,layer_params);
    gps_time = [gps_time, layers.gps_time];
    twtt = [twtt, layers.twtt];
  end
  
  copy_param.layer_source.source = 'custom';
  copy_param.layer_source.gps_time = gps_time;
  copy_param.layer_source.twtt = twtt;
  
  copy_param.layer_dest.name = 'surface_rds';
  copy_param.layer_dest.source = 'layerdata';
  copy_param.layer_dest.layerdata_source = 'layerData';
  copy_param.layer_dest.existence_check = false;

  copy_param.copy_method = 'overwrite';
  copy_param.gaps_fill.method = 'interp_finite';

end

% =====================================================================
%% Automated Section
% =====================================================================

%% Load each of the day segments
global gRadar;
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  param = merge_structs(param,gRadar);
  fprintf('opsCopyLayers %s (%s)\n', param.day_seg, datestr(now));
  opsCopyLayers(param,copy_param);
  fprintf('  Complete (%s)\n', datestr(now));
end
