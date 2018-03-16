% script test_get_heights
%
% Script for testing get_heights.
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_get_heights.m, get_heights.m,
%   get_heights_task.m

%% Test Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'),'20161024_05');

params.cmd.get_heights = 1;
params.cmd.frms = 1;

params.get_heights = [];
params.get_heights.out_path = '';
params.get_heights.block_size = 20000;
params.get_heights.img_comb = [3.0000e-06 -Inf 1.0000e-06 1.0000e-05 -Inf 3.0000e-06];
params.get_heights.imgs = {[ones([6 1]),(1:6).'],[2*ones([6 1]),(1:6).'],[3*ones([6 1]),(1:6).']};
params.get_heights.ft_wind = @hanning;
params.get_heights.lever_arm_fh = @lever_arm;
params.get_heights.B_filter = ones(20,1)./20;
params.get_heights.decimate_factor = 20;
params.get_heights.inc_ave = 10;
params.get_heights.surf.en = 1;
params.get_heights.surf.method = 'threshold';
params.get_heights.surf.noise_rng = [0 -50 -10];
params.get_heights.surf.min_bin = 1.5000e-06;
params.get_heights.surf.max_bin = [];
params.get_heights.surf.threshold = 9;
params.get_heights.surf.sidelobe = 15;
params.get_heights.surf.medfilt = 3;
params.get_heights.surf.max_diff = 1.0000e-07;
params.get_heights.surf.search_rng = [-5 5];

param_override.out_path = fullfile(gRadar.out_path,'tests');
param_override.cluster.type = 'debug';
param_override.cluster.rerun_only = false;

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
  if param.cmd.get_heights
    ctrl_chain{end+1} = get_heights(param,param_override);
  end
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);

[ctrl_chain,chain_fn] = cluster_load_chain([],chain_id);
ctrl_chain = cluster_run(ctrl_chain);
