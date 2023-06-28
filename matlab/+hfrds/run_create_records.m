% script hfrds.run_create_records
%
% Script for running hfrds.create_records (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, hfrds.run_create_records.m, hfrds.create_records.m

% =====================================================================
% Debug Setup
% =====================================================================
param = read_param_xls(ct_filename_param('rds_param_2016_Greenland_G1XB.xls'),'20160413_01');

clear('param_override');
dbstop if error
param_override.sched.type = 'no scheduler';
param_override.sched.rerun_only = true;

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

hfrds.create_records(param,param_override);
