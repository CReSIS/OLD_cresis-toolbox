% script run_combine_wf_chan
%
% Script for running combine_wf_chan (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_combine_wf_chan.m, combine_wf_chan.m,
%   combine_wf_chan_task.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('accum_param_2018_Greenland_P3.xls'),'');

% Syntax for running a specific segment and frame by overriding parameter spreadsheet values
%params = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'),'20161024_05');
% params = ct_set_params(params,'cmd.combine_wf_chan',0);
% params = ct_set_params(params,'cmd.combine_wf_chan',1,'day_seg','20161102_03');
% params = ct_set_params(params,'cmd.frms',6);

params = ct_set_params(params,'cmd.combine_wf_chan',0);
params = ct_set_params(params,'cmd.combine_wf_chan',1,'day_seg','20180405_01');
params = ct_set_params(params,'cmd.frms',[115]);
% params = ct_set_params(params,'cmd.csarp',1,'day_seg','20161024_05');
% params = ct_set_params(params,'combine.in_path','out_cn');
% params = ct_set_params(params,'combine.array_path','out_cn');
% params = ct_set_params(params,'combine.out_path','standard_cn');

% params = ct_set_params(params,'cmd.combine_wf_chan',0);
% params = ct_set_params(params,'cmd.combine_wf_chan',1,'day_seg','20161107_02');
% params = ct_set_params(params,'cmd.combine',1,'day_seg','20161024_05');
% params = ct_set_params(params,'cmd.frms',2,'day_seg','20161107_02');
% params = ct_set_params(params,'combine.out_path','paden_standard');
% params = ct_set_params(params,'combine.in_path','paden_out');
% params = ct_set_params(params,'combine.array_path','paden_out');

dbstop if error;
% param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
param_override.cluster.type = 'debug';
%param_override.cluster.rerun_only = true;
param_override.cluster.desired_time_per_job  = 0*60;
% param_override.cluster.cpu_time_mult  = 2;
% param_override.cluster.mem_mult  = 2;
% param_override.cluster.max_jobs_active       = 1;
% param_override.cluster.qsub_submit_arguments = '-q debug -m n -l nodes=1:ppn=1:dcwan:dc2,pmem=%dmb,walltime=%d:00';

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
  if param.cmd.combine_wf_chan
    ctrl_chain{end+1} = combine_wf_chan(param,param_override);
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
