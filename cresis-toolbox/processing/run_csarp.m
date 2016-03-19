% script run_csarp
%
% Script for running csarp (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_csarp.m, csarp.m,
%   csarp_task.m

% =====================================================================
% Debug Setup
% =====================================================================
param = read_param_xls(ct_filename_param('rds_param_2014_Antarctica_DC8.xls'),'20141121_05');

clear('param_override');
param_override.sched.type = 'no scheduler';
param_override.sched.rerun_only = false;

% Input checking
if ~exist('param','var')
  error('A struct array of parameters must be passed in\n');
end
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

csarp(param,param_override);

return;

