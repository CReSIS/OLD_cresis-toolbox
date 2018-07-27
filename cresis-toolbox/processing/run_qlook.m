% script run_qlook
%
% Script for running qlook (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_qlook.m, qlook.m,
%   qlook_task.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('snow_param_2017_Greenland_P3.xls'),'20170311_02');

% Syntax for running a specific segment and frame by overriding parameter spreadsheet values
%params = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'),'20161024_05');
params = ct_set_params(params,'cmd.qlook',0);
% params = ct_set_params(params,'cmd.qlook',1,'day_seg','20170503_02');
% params = ct_set_params(params,'cmd.frms',[221]);
params = ct_set_params(params,'cmd.qlook',1,'day_seg','20170311_02');
% params = ct_set_params(params,'cmd.frms',[7]);
% params = ct_set_params(params,'qlook.presums',4);
% params = ct_set_params(params,'qlook.dec',1);
% params = ct_set_params(params,'qlook.out_path','qlook');
% params = ct_set_params(params,'radar.wfs(1).coh_noise_method','');

dbstop if error;
param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
%param_override.cluster.rerun_only = true;
param_override.cluster.desired_time_per_job  = 120*60;
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

% Potentially stop and inspect cluster_print_chain output to adjust
% cluster control parameters before running or to run the next lines on a
% different computer (the save/load functions are for this purpose).

return
%ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.desired_time_per_job',5*60);
%ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.cpu_time_mult',2);
%ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.mem_mult',2);

[ctrl_chain,chain_fn] = cluster_load_chain([],chain_id);
ctrl_chain = cluster_run(ctrl_chain);
