function cluster_job(task_in_fn_dir,task_out_fn_dir,job_list)
%
% This M-file should be compiled:
%   mcc -m -C worker_task.m
% Then the worker_task and worker_task.ctf files should be moved
% to a spot in the PATH that comes before Matlab's worker executable.
% The worker shell script should also be placed in this same
% directory.
%
% Then when the distributed toolbox calls the worker command via the
% torque scheduler, it calls the worker shell script since it comes
% before the Matlab's worker executable. The worker shell script
% calls the compiled version of this function.
%
% demo_matlab_torque.m shows how to use the scheduler from Matlab
%
% In your bashrc file:
%  declare -x PATH="/N/u/jpaden/Quarry/bin:${PATH}"
% Compile the worker_task.m and put the two outputs into the folder.
% Put the shell script worker into this folder.

fprintf('%s: Start (%s)\n', mfilename, datestr(now));

if ~exist('task_in_fn_dir','var')
  task_in_fn_dir = getenv('INPUT_PATH');
end
if ~exist('task_out_fn_dir','var')
  task_out_fn_dir = getenv('OUTPUT_PATH');
end
if ~exist('job_list','var')
  job_list = getenv('JOB_LIST');
end
job_list = regexp(job_list, 'd', 'split');

for task_idx = 1:length(job_list)
  task_id = str2double(job_list{task_idx});
  fprintf('%s: Load task %d (%s)\n', mfilename, task_id, datestr(now));
  
  % Create in/out filenames
  static_in_fn = fullfile(task_in_fn_dir,'static.mat');
  dynamic_in_fn = fullfile(task_in_fn_dir,'dynamic.mat');
  out_fn = fullfile(task_out_fn_dir,sprintf('out_%d.mat',task_id));
  
  % Read input param struct
  %   param fields: 'taskfunction','argsin','num_args_out'
  sparam = load(static_in_fn);
  dparam_task_field = sprintf('dparam_%d',task_id);
  dparam = load(dynamic_in_fn,dparam_task_field);
  param = merge_structs(sparam.static_param,dparam.(dparam_task_field));
  
  % Creating command to evaluate
  if param.num_args_out == 0
    eval_cmd = sprintf('%s(',param.task_function);
  else
    eval_cmd = sprintf('[argsout{1:%i}] = %s(',param.num_args_out,param.task_function);
  end
  for argsin_idx = 1:length(param.argsin)
    if argsin_idx < length(param.argsin)
      eval_cmd = sprintf('%sparam.argsin{%i},', eval_cmd, argsin_idx);
    else
      eval_cmd = sprintf('%sparam.argsin{%i}', eval_cmd, argsin_idx);
    end
  end
  eval_cmd = sprintf('%s);', eval_cmd);
  
  % Evaluate command
  try
    argsout = {};
    fprintf('%s: Eval %s\n', mfilename, eval_cmd);
    eval(eval_cmd);
    fprintf('%s: Done Eval (%s)\n', mfilename, datestr(now));
    errorstruct = [];
    save(out_fn,'argsout','errorstruct');
  catch errorstruct
    fprintf('%s: Error\n  %s: %s (%s)\n', mfilename, errorstruct.identifier, errorstruct.message, datestr(now));
    for stack_idx = 1:length(errorstruct.stack)
      fprintf('  %s: %d\n', errorstruct.stack(stack_idx).name, errorstruct.stack(stack_idx).line);
    end
    save(out_fn,'argsout','errorstruct');
  end
end

return;
