% Script run_radiometric
% 
% Loads the "radiometric" worksheet from the parameter spreadsheet and then calls
% radiometric_calibration with this information.
%
% Authors: Theresa Stumpf, Anthony Hoch
% note: modified from run_post
%
% See also: radiometric_calibration, run_post.m

fprintf('\n\n========================================================\n');
fprintf('run radiometric\n');
fprintf('========================================================\n');

%% User Settings

% params = read_param_xls('/users/hocha/scripts/personal-params/mcrds_param_2008_Greenland_TO.xls',[],'radiometric');
params = read_param_xls('/users/tstumpf/scripts/matlab/rds_param_2011_Antarctica_TO.xls',[],'radiometric');

clear('param_override');
param_override = [];
% param_override.sched.type = 'no scheduler';
param_override.sched.cluster_size = inf;
%param_override.sched.submit_arguments    = '-b 240 -q pg@qm2 -e "/N/dcwan/scratch/jpaden/matlab_error" -l nodes=1:ppn=4:dcwan:gpfs,walltime=15:00';
%param_override.sched.stop_on_fail = true;



%% Automated Section
% =====================================================================
% Create param structure array
% =====================================================================
tic;
global gRadar;

% Input checking
if ~exist('params','var')
  error('Check Param Status: param spreadsheet did not load correctly\n');
  return
end
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% =====================================================================
% =====================================================================
% Process each of the days
% =====================================================================
% =====================================================================
for param_idx = 1:length(params)
  param = params(param_idx);
  cmd = param.cmd;
  if isfield(cmd,'generic') && cmd.generic
    radiometric_calibration(param,param_override);
  end  
end


return
