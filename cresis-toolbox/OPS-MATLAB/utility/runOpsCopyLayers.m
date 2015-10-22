% script runOpsCopyLayers.m
%
% Runs opsCopyLayers.m

%% User Settings
params = read_param_xls(ct_filename_param('rds_param_2009_Greenland_TO.xls'),'','post');

copy_param = [];

copy_param.layer_source.name = 'surface';
copy_param.layer_source.source = 'ops';
copy_param.layer_source.echogram_source = 'qlook';
copy_param.layer_source.layerdata_source = 'layerData';
copy_param.layer_source.existence_check = true;

copy_param.layer_dest.name = 'surface';
copy_param.layer_dest.source = 'records';
copy_param.layer_dest.echogram_source = 'qlook';
copy_param.layer_dest.layerdata_source = 'layerData';
copy_param.layer_dest.existence_check = true;

copy_param.copy_method = 'overwrite';
copy_param.gaps_fill.method = 'interp_finite';
  
%% Automated Section

%% Load each of the day segments
% =====================================================================
layers = {};
day_seg = {};
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