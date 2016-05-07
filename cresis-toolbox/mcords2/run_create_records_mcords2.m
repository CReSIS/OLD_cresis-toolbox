% script run_create_records_mcords2
%
% Script for running create_records_mcords2 (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_create_records_mcords2.m, create_records_mcords2.m,
%   create_records_mcords2_sync.m

% =====================================================================
% Debug Setup
% =====================================================================
param = read_param_xls(ct_filename_param('rds_param_2016_Greenland_Polar6.xls'),'20160426_11');

clear('param_override');
dbstop if error
param_override.sched.type = 'no scheduler';
param_override.sched.rerun_only = true;

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

create_records_mcords2(param,param_override);

return;
