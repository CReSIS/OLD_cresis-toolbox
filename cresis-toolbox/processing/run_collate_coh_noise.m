% script run_collate_coh_noise
%
% Runs collate_coh_noise
%
% Authors: John Paden

%% USER SETTINGS
% =========================================================================

param_override = [];

% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),'',{'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),'',{'analysis_noise','analysis'});
params = read_param_xls(ct_filename_param('snow_param_2013_Greenland_P3.xls'),'',{'analysis_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2014_Greenland_P3.xls'),'',{'analysis_noise','analysis'});

% Enable a specific segment
params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20120330_04');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20130420_01');
params = ct_set_params(params,'cmd.generic',1,'day_seg','20130328_01'); %random test for reuse_debug in collate_coh_noise
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20140428_01');

if 1
  % Near-DC removal
  param_override.collate_coh_noise.method{1} = 'firdec';
  param_override.collate_coh_noise.firdec_fcutoff{1} = @(t) 1/30; % Update coherent noise estimate every 30 seconds
  param_override.collate_coh_noise.firdec_fs{1} = 1/7.5; % Should update about 4 times as often as the estimate: 30/4 = 7.5
else
  % DC removal when dft_corr_time set to inf
  param_override.collate_coh_noise.method = 'dft';
  param_override.collate_coh_noise.dft_corr_time = inf;
end
% param_override.collate_coh_noise.in_path = 'analysis_threshold'; % Enable during second pass
% param_override.collate_coh_noise.out_path = 'analysis_threshold'; % Enable during second pass

param_override.collate_coh_noise.cmd_idx = 1;
param_override.collate_coh_noise.debug_plots = {'visible','cn_plot','threshold_plot'}; % Debugging
param_override.collate_coh_noise.debug_plots = {'visible','cn_plot','threshold_plot','reuse'}; % Debugging
% param_override.collate_coh_noise.debug_plots = {'cn_plot','threshold_plot'}; % Typical setting when not debugging
% param_override.collate_coh_noise.debug_plots = {}; % Necessary if plots are too large for memory

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
