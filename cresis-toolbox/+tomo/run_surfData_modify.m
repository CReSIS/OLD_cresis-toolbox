% script tomo.run_surfData_modify.m
%
% Example script for running opsLoadLayers.m. Demonstrates a few of the
% most common operations to be performed with opsLoadLayers and supports
% interpolation of various layer sources for better comparison.
%
% Authors: John Paden

% =====================================================================
%% User Settings
% =====================================================================

params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'));

layers = [3 4];
args{1} = 'control_layer';
args{2} = 4;
args{3} = 'mask_layer';
args{4} = 3;

% =====================================================================
%% Automated Section
% =====================================================================

global gRadar;

%% Load each of the day segments
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) ...
      || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  param = merge_structs(param,gRadar);
  
  fprintf('surfData_modify %s\n', param.day_seg);
  tomo.surfData_modify(param,'',layers,args{:});
end
