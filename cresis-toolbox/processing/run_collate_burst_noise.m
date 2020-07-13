% script run_collate_burst_noise
%
% Runs collate_burst_noise
%
% Authors: John Paden

%% USER SETTINGS
% =========================================================================

param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2019_Antarctica_Ground.xls'),'20200105_02',{'analysis_noise','analysis'});

param_override.collate_burst_noise.in_path = 'analysis_burst';

param_override.collate_burst_noise.debug_plots = {'bn_plot'};
% param_override.collate_burst_noise.debug_plots = {'visible','bn_plot'}; % <== CHOOSE to debug

if 0
  % For debugging, use this to select a specific image and wf_adc to
  % collate instead of doing them all
  param_override.collate_burst_noise.wf_adcs{img} = 4;
  param_override.collate_burst_noise.imgs = 2;
end

cmd_method = 'generic';
rds_settings;

params.cmd.generic = true;

%% Automated Section
% =====================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Process each of the segments
ctrl_chain = {};
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
%   collate_burst_noise(param,param_override);
  collate_burst_noise
  
end
