% script test_combine_wf_chan
%
% Script for testing combine_wf_chan.
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_combine_wf_chan.m, combine_wf_chan.m,
%   combine_wf_chan_task.m

%% Test Setup
% =====================================================================
param = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'),'20161024_05');

param = ct_set_params(param,'cmd.frms',[1]);

param.combine = [];
param.combine.in_path = '';
param.combine.array_path = '';
param.combine.out_path = '';
param.combine.img_comb = [3.0000e-06 -Inf 1.0000e-06 1.0000e-05 -Inf 3.0000e-06];
param.combine.imgs = {[ones([6 1]),(1:6).'],[2*ones([6 1]),(1:6).'],[3*ones([6 1]),(1:6).']};
param.combine.method = 'standard';
param.combine.window = @hanning;
param.combine.bin_rng = [-1 0 1];
param.combine.rline_rng = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
param.combine.debug_level = 0;
param.combine.dbin = 1;
param.combine.dline = 6;
param.combine.DCM = [];
param.combine.three_dim.en = 0;
param.combine.Nsv = 1;
param.combine.theta_rng = [0 0];
param.combine.sv_fh = @array_proc_sv;
param.combine.diag_load = 0;
param.combine.Nsig = 2;

param.sched.type = 'no scheduler';
param.sched.rerun_only = false;

clear('param_override');
param_override.out_path = fullfile(gRadar.out_path,'tests');

%% Automated Section
% =====================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

combine_wf_chan(param,param_override);
