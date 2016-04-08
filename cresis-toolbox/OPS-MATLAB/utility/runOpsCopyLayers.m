% script runOpsCopyLayers.m
%
% Example script for running opsCopyLayers.m. Demonstrates a few of the
% most common operations to be performed with opsCopyLayers.
%
% Authors: John Paden

% =====================================================================
%% User Settings
% =====================================================================

params = read_param_xls(ct_filename_param('rds_param_2013_Antarctica_P3.xls'),'','post');
% params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'','post');
% params = read_param_xls(ct_filename_param('rds_param_2015_Greenland_Polar6.xls'),'','post');
% params = read_param_xls(ct_filename_param('rds_param_2016_Greenland_Polar6.xls'),'','post');
% params = read_param_xls(ct_filename_param('snow_param_2009_Greenland_P3.xls'),'','post');
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_DC8.xls'),'','post');
% params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'),'','post');
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),'','post');
% params = read_param_xls(ct_filename_param('snow_param_2009_Antarctica_DC8.xls'),'','post');
% params = read_param_xls(ct_filename_param('snow_param_2010_Antarctica_DC8.xls'),'','post');
% params = read_param_xls(ct_filename_param('snow_param_2011_Antarctica_DC8.xls'),'','post');
% params = read_param_xls(ct_filename_param('snow_param_2012_Antarctica_DC8.xls'),'','post');
% params = read_param_xls(ct_filename_param('snow_param_2015_Greenland_Polar6.xls'),'','post');

copy_param = [];

%% Example Operations (just choose one)

if 1
  %% Copy surface from echogram qlook to records (and lots of other example sources/destinations)
  copy_param.layer_source.name = 'surface';
  copy_param.layer_source.existence_check = false;
  
  copy_param.layer_source.source = 'ops';
%   copy_param.layer_source.source = 'echogram';
%   copy_param.layer_source.echogram_source = 'qlook';
%   copy_param.layer_source.source = 'echogram';
%   copy_param.layer_source.echogram_source = 'deconv';
%   copy_param.layer_source.source = 'GIMP';
%   copy_param.layer_source.source = 'layerdata';
%   copy_param.layer_source.layerdata_source = 'layerData';
%   copy_param.layer_source.source = 'lidar';
%   copy_param.layer_source.lidar_source = 'awi';
  

  copy_param.layer_dest.name = 'surface';
  copy_param.layer_dest.existence_check = false;
  copy_param.copy_method = 'overwrite';
  copy_param.gaps_fill.method = 'interp_finite';
  
%   copy_param.layer_dest.source = 'records';
  copy_param.layer_dest.source = 'layerdata';
  copy_param.layer_dest.layerdata_source = 'layerData';
%   copy_param.layer_dest.source = 'echogram';
%   copy_param.layer_dest.echogram_source = 'qlook';
%   copy_param.layer_dest.source = 'echogram';
%   copy_param.layer_dest.echogram_source = 'deconv';
%   copy_param.layer_dest.source = 'echogram';
%   copy_param.layer_dest.echogram_source = 'standard';

elseif 0
  %% Copy surface from layerdata to records
  copy_param.layer_source.name = 'surface';
  copy_param.layer_source.source = 'layerdata';
  copy_param.layer_source.layerdata_source = 'layerData';
  
  copy_param.layer_dest.name = 'surface';
  copy_param.layer_dest.source = 'records';

  copy_param.copy_method = 'merge';
  copy_param.gaps_fill.method = 'interp_finite';
  
elseif 0
  %% Copy surface from ops to records
  copy_param.layer_source.name = 'surface';
  copy_param.layer_source.source = 'ops';
  
  copy_param.layer_dest.name = 'surface';
  copy_param.layer_dest.source = 'records';

  copy_param.copy_method = 'merge';
  copy_param.gaps_fill.method = 'interp_finite';
  
elseif 0
  %% Copy layer from ops to layerdata
%   layer_name = 'surface';
  layer_name = 'bottom';
  
  copy_param.layer_source.name = layer_name;
  copy_param.layer_source.source = 'ops';
  
  copy_param.layer_dest.name = layer_name;
  copy_param.layer_dest.source = 'layerdata';
  copy_param.layer_dest.layerdata_source = 'layerData';

  copy_param.copy_method = 'overwrite';
  copy_param.gaps_fill.method = 'preserve_gaps';
  
elseif 0
  %% Copy bottom from ops to echogram CSARP_post/mvdr
  copy_param.layer_source.name = 'bottom';
  copy_param.layer_source.source = 'ops';
  
  copy_param.layer_dest.name = 'bottom';
  copy_param.layer_dest.source = 'echogram';
  copy_param.layer_dest.echogram_source = 'standard';

  copy_param.copy_method = 'merge';
  copy_param.gaps_fill.method = 'preserve_gaps';
  copy_param.gaps_fill.method_args = [300 60];
    
elseif 0
  %% Copy "tomo_bottom2" layer from ops to echogram CSARP_post/mvdr
  copy_param.layer_source.name = 'tomo_bottom2';
  copy_param.layer_source.source = 'ops';
  
  copy_param.layer_dest.name = 'tomo_bottom2';
  copy_param.layer_dest.source = 'echogram';
  copy_param.layer_dest.existence_check = false;

  copy_param.copy_method = 'merge';
  copy_param.gaps_fill.method = 'preserve_gaps';
  copy_param.gaps_fill.method_args = [300 60];
  
elseif 0
  %% Copy "custom" layer from layerdata to "new" layer in ops
  copy_param.layer_source.name = 'custom';
  copy_param.layer_source.source = 'layerdata';
  copy_param.layer_source.layerdata_source = 'layerData';
  
  copy_param.layer_dest.name = 'new';
  copy_param.layer_dest.source = 'ops';
  copy_param.layer_dest.existence_check = false;

  copy_param.copy_method = 'merge';
  copy_param.gaps_fill.method = 'interp_finite';
  
elseif 0
  %% Copy surface from mcords records to snow layerdata  
  
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
  
elseif 0
  %% Copy surface from echogram deconv to echogram qlook
  copy_param.layer_source.name = 'surface';
  copy_param.layer_source.source = 'echogram';
  copy_param.layer_source.echogram_source = 'deconv';
  
  copy_param.layer_dest.name = 'surface';
  copy_param.layer_dest.source = 'echogram';
  copy_param.layer_dest.echogram_source = 'qlook';

  copy_param.copy_method = 'merge';
  copy_param.gaps_fill.method = 'interp_finite';
  
elseif 0
  %% Shift records surface
  copy_param.layer_source.name = 'surface';
  copy_param.layer_source.source = 'records';
  
  copy_param.layer_dest.name = 'surface';
  copy_param.layer_dest.source = 'records';
  
  twtt_offset = -2.424000000000000e-07;
  copy_param.eval.cmd = sprintf('source = source + %.12g;',twtt_offset);

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
  fprintf('opsCopyLayers %s\n', param.day_seg);
  opsCopyLayers(param,copy_param);
end