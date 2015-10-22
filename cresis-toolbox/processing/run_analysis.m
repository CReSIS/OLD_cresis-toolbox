% Script run_analysis
% 
% run_analysis description
fprintf('\n\n========================================================\n');
fprintf('run analysis\n');
fprintf('========================================================\n');

%% User Settings

params = read_param_xls('/users/paden/scripts/cresis-toolbox/params/mcords_param_2011_Greenland_P3.xls',[],'analysis');

clear('param_override');
param_override = [];
param_override.sched.type = 'no scheduler';
param_override.sched.cluster_size = 15;
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
  error('Use run_master: A struct array of parameters must be passed in\n');
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
    analysis(param,param_override);
  end  
end


