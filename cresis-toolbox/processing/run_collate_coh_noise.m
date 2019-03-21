% script run_collate_coh_noise
%
% Runs collate_coh_noise
%
% Authors: John Paden

%% USER SETTINGS
% =========================================================================

param_override = [];

% params = read_param_xls(ct_filename_param('rds_param_2016_Greenland_TOdtu.xls'),'','analysis');
% params = read_param_xls(ct_filename_param('snow_param_2017_Greenland_P3.xls'),'',{'analysis_noise','analysis'});
params = read_param_xls(ct_filename_param('snow_param_2017_Arctic_Polar5.xls'),'',{'analysis_noise','analysis'});

params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20170330_01');

if 1
  param_override.collate_coh_noise.method = 'firdec';
  param_override.collate_coh_noise.firdec_fs = 1/7.5;
  param_override.collate_coh_noise.firdec_fcutoff = @(t) 1/30;
else
  param_override.collate_coh_noise.method = 'dft';
  param_override.collate_coh_noise.dft_corr_length = inf;
end
param_override.collate_coh_noise.in_dir = 'analysis';
param_override.collate_coh_noise.out_dir = 'analysis';

param_override.collate_coh_noise.min_samples = 1500;
param_override.collate_coh_noise.threshold_en = true;

% param_override.collate_coh_noise.debug_plots = {};
param_override.collate_coh_noise.debug_plots = {'visible','cn_plot','threshold_plot'};
% param_override.collate_coh_noise.debug_plots = {'cn_plot'};

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
%   collate_coh_noise(param,param_override);
  collate_coh_noise
  
end
