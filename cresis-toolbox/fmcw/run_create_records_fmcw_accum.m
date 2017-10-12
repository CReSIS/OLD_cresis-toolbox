% script run_create_records_fmcw_accum
%
% Script for running create_records_fmcw_accum (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_create_records_fmcw_accum.m, create_records_fmcw_accum.m,
%   create_records_fmcw_accum_sync.m, check_records.m

% =====================================================================
% Debug Setup
% =====================================================================
param = read_param_xls(ct_filename_param('snow_param_2016_Antarctica_DC8.xls'),'20161020_02');

clear('param_override');
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

create_records_fmcw_accum(param,param_override);

return;
