function doa_cluster(param,param_override)
% This function is called from Sim.crosstrack.
%
% Inputs:
%   param = struct with processing parameters
%

% clear('param_override');
% param_override = [];
% physical_constants;
%
% global gRadar;
%
% % param_override.sched.type = 'no scheduler';
% param_override.sched.type = 'custom_torque';
%
% param_override.sched.cluster_size = inf;
% param_override.sched.rerun_only   = false;
% param_override.sched.stop_on_fail = false;
%
% param_override.sched.num_submissions = 170;
% param_override.sched.group_size      = 1;
% param_override.sched.group_walltime  = 2*86400;
%
% param_override.sched.submit_arguments = '-l nodes=1:ppn=1,pmem=14000mb,walltime=48:00:00';
% % param_override.sched.submit_arguments = '-l nodes=1:ppn=1,pmem=8000mb,walltime=120:00';
%
% if exist('param_override','var')
%   param_override = merge_structs(gRadar,param_override);
% else
%   param_override = gRadar;
% end

param = merge_structs(param,param_override);

%% Initialize Torque setup
if strcmpi(param.sched.type,'custom_torque')
  global ctrl; % Make this global for convenience in debugging
  ctrl = torque_new_batch(param);
  fprintf('Torque batch: %s\n', ctrl.batch_dir);
  torque_compile('crosstrack.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
end

%% Create Tasks
task_param = param;

fh  = @crosstrack;
arg = {task_param};

if strcmp(param.sched.type,'custom_torque')
  create_task_param.conforming = true;
  %   create_task_param.notes = sprintf(' ');
  ctrl = torque_create_task(ctrl,fh,1,arg,create_task_param);
else
  fh(arg{:});
end

%% Wait for tasks to finish
if strcmpi(param.sched.type,'custom_torque')
  % Wait until all submitted jobs to complete
  ctrl = torque_rerun(ctrl);
  
  if ~all(ctrl.error_mask == 0)
    if ctrl.sched.stop_on_fail
      torque_cleanup(ctrl);
      error('Not all jobs completed, but out of retries (%s)', datestr(now));
    else
      warning('Not all jobs completed, but out of retries (%s)', datestr(now));
      keyboard;
    end
  else
    fprintf('Jobs completed (%s)\n\n', datestr(now));
  end
  
  torque_cleanup(ctrl);
end

end