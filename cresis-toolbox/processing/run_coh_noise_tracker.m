% script run_coh_noise_tracker
%
% Script for running coh_noise_tracker
%
% Authors: John Paden
%
% See also: master.m, run_coh_noise_tracker.m coh_noise_tracker.m,
%   coh_noise_tracker_task.m

% =====================================================================
% Debug Setup
% =====================================================================
% param = read_param_xls(ct_filename_param('snow_param_2009_Antarctica_DC8.xls'),'20091107_05',{'analysis_spec','analysis'});
% param = read_param_xls(ct_filename_param('snow_param_2009_Antarctica_DC8.xls'),'20091118_01',{'analysis_coh_noise','analysis'});
param = read_param_xls(ct_filename_param('snow_param_2010_Greenland_P3.xls'),'20100507_01',{'analysis_coh_noise','analysis'});

clear('param_override');
dbstop if error;
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

coh_noise_tracker(param,param_override);

return
