% script run_qlook
%
% Script for running qlook.m (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_qlook.m, qlook.m,
%   qlook_task.m

%% User Setup
% =====================================================================
param_override = [];

% params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'));
% params = read_param_xls(ct_filename_param('rds_param_2018_Antarctica_Ground.xls'));
params = read_param_xls(ct_filename_param('snow_param_2018_Greenland_P3.xls'));

% Example to run specific segments and frames by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.qlook',0);
params = ct_set_params(params,'cmd.qlook',1,'day_seg','20180320_01');
params = ct_set_params(params,'cmd.frms',[6:73]);

% param_override.records.file.base_dir = '/mnt/nfs_raid_ssdorange/2018_Greenland_P3/';

params = ct_set_params(params,'radar.wfs(1).coh_noise_method',''); % before analysis
params = ct_set_params(params,'radar.wfs(1).coh_noise_method','analysis'); % after analysis

params = ct_set_params(params,'radar.wfs(1).deconv.en',false);

param_override.qlook.out_path = 'qlook_noise';
param_override.qlook.surf.en = false;
% params = ct_set_params(params,'radar.wfs(1).t_ref', -0.000000071980808 -38.6e-6 );

dbstop if error;
% param_override.cluster.type = 'torque';
param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
% param_override.cluster.type = 'slurm';
% param_override.cluster.rerun_only = true;
% param_override.cluster.desired_time_per_job  = 240*60;
% param_override.cluster.cpu_time_mult  = 2;
% param_override.cluster.mem_mult  = 2;

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
  if param.cmd.qlook
    ctrl_chain{end+1} = qlook(param,param_override);
  end
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);
