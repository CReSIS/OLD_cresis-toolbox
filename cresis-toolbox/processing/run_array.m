% script run_array
%
% Script for running array (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_array.m, array.m, load_sar_data.m,
% array_proc.m, array_task.m, array_combine_task.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2018_Greenland_P3.xls'),'');

% Syntax for running a specific segment and frame by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.array',0);
params = ct_set_params(params,'cmd.array',1,'day_seg','20180404_02');
params = ct_set_params(params,'cmd.frms',[]);
params = ct_set_params(params,'array.method','mvdr');
%params = ct_set_params(params,'array.DCM',[]);
params = ct_set_params(params,'array.bin_rng',-1:1);
params = ct_set_params(params,'array.rline_rng',-9:9);

params = ct_set_params(params,'array.method','music');
params = ct_set_params(params,'array.out_path','music_imgs4_Nsig2');
params = ct_set_params(params,'array.three_dim.en',1);
params = ct_set_params(params,'array.Nsv',64);
params = ct_set_params(params,'array.Nsig',2);
params = ct_set_params(params,'array.imgs',{[ones([7 1]),[6:12].'],[2*ones([7 1]),[6:12].'],[3*ones([7 1]),[6:12].'],[4*ones([7 1]),[6:12].']});
params = ct_set_params(params,'array.img_comb',[]);


for param_idx = 1:length(params)
  %param = params(param_idx);
  for wf = 1:length(params(param_idx).radar.wfs)
    params(param_idx).radar.wfs(wf).Tsys = [0.46 -4.66 0.14 -1.77 0 -2.63 -3.38 -69.66 -75.57 -75.45 -80.42 -80.49 -75.71 -77.69 -70.53]/1e9;
    params(param_idx).radar.wfs(wf).chan_equal_dB = [6.8 -0.6 3 0.1 0 3.5 3.9 7 3.3 4.8 6.1 6.2 4.6 3.1 6.2];
    params(param_idx).radar.wfs(wf).chan_equal_deg = [-166.2 -142.7 177 -95.9 0 -25.9 -86.5 -27.4 128.1 41.6 -46.8 43 90.7 121.3 31.6];
    params(param_idx).radar.wfs(wf).coh_noise_method = 'analysis';
    params(param_idx).radar.wfs(wf).deconv.en = 0;
    params(param_idx).radar.wfs(wf).deconv.fn = 'analysis';
  end
end

dbstop if error;
param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
%param_override.cluster.rerun_only = true;
param_override.cluster.desired_time_per_job  = 10*60;
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
  if param.cmd.array
    ctrl_chain{end+1} = array(param,param_override);
  end
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);
