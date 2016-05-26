% script run_combine_wf_chan
%
% Script for running combine_wf_chan (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_combine_wf_chan.m, combine_wf_chan.m,
%   combine_wf_chan_task.m

% =====================================================================
% Debug Setup
% =====================================================================
param = read_param_xls(ct_filename_param('rds_param_2016_Greenland_G1XB.xls'),'20160413_01');

dbstop if error
param.cmd.frms = [];
param.cmd.combine = true;

clear('param_override');
param_override.sched.type = 'no scheduler';
param_override.sched.rerun_only = true;

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

combine_wf_chan(param,param_override);
