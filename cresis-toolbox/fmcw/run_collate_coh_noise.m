% script run_collate_coh_noise
%
% Runs collate_coh_noise
%
% Authors: John Paden

%% USER SETTINGS
% =========================================================================

param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2016_Greenland_TOdtu.xls'),'','analysis');
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20161110_01');
% params = ct_set_params(params,'cmd.generic',1,'cmd.notes','^((?!do not process).)*$');

param_override.collate_coh_noise.in_dir = 'analysis';
param_override.collate_coh_noise.out_dir = 'analysis';
param_override.collate_coh_noise.cmd_idx = 1;
param_override.collate_coh_noise.debug_level = 0; % Set to 1 to view threshold and minimum samples result, set to 2 to also view Doppler plots

imgs = [1];
wf_adcs = [1];

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
  for img = imgs
    for wf_adc = wf_adcs
      param_override.collate_coh_noise.img = img;
      param_override.collate_coh_noise.wf_adc = wf_adc;
      collate_coh_noise(param,param_override);
    end
  end
  
end