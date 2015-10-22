function worker_task
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

fprintf('worker_task\n');

custom_torque = getenv('CUSTOM_TORQUE');
if ~isempty(custom_torque)
  fprintf('Custom torque variable found\n');
  task_in_fn_dir = getenv('INPUT_PATH');
  task_out_fn_dir = getenv('OUTPUT_PATH');
  job_list = getenv('JOB_LIST');
  job_list = regexp(job_list, 'd', 'split');
  task_in_fn = {};
  task_out_fn = {};
  for job_idx = 1:length(job_list)
    task_in_fn{job_idx} = fullfile(task_in_fn_dir,sprintf('in_%s.mat',job_list{job_idx}));
    task_out_fn{job_idx} = fullfile(task_out_fn_dir,sprintf('out_%s.mat',job_list{job_idx}));
  end
else
  mdce.storage_location = getenv('MDCE_STORAGE_LOCATION');
  mdce.job_location = getenv('MDCE_JOB_LOCATION');
  mdce.task_id = getenv('MDCE_TASK_ID');

  startInd = regexp(mdce.storage_location,'UNIX{') + length('UNIX{');
  stopInd = startInd + regexp(mdce.storage_location(startInd:end),'}') - 2;
  task_in_fn = {fullfile(mdce.storage_location(startInd:stopInd),mdce.job_location, ...
    sprintf('Task%s.in.mat',mdce.task_id))};
  task_out_fn = {fullfile(mdce.storage_location(startInd:stopInd),mdce.job_location, ...
    sprintf('Task%s.out.mat',mdce.task_id))};
end

for fn_idx = 1:length(task_in_fn)
  % task fields: 'taskfunction','argsin','num_args_out','param'
  fprintf('Loading input arguments and control parameters\n');
  task = load(task_in_fn{fn_idx});
  
  funch = func2str(task.taskfunction);
  
  fprintf('Creating command to evaluate\n');
  if task.num_args_out == 0
    eval_cmd = sprintf('funch(');
  else
    eval_cmd = sprintf('[argsout{1:%i}] = %s(',task.num_args_out,funch);
  end
  for argsin_idx = 1:length(task.argsin)
    if argsin_idx < length(task.argsin)
      eval_cmd = sprintf('%stask.argsin{%i},', eval_cmd, argsin_idx);
    else
      eval_cmd = sprintf('%stask.argsin{%i}', eval_cmd, argsin_idx);
    end
  end
  eval_cmd = sprintf('%s);', eval_cmd);
  
  try
    argsout = {};
    fprintf('Eval %s\n', eval_cmd);
    eval(eval_cmd);
    fprintf('Done eval\n');
    save(task_out_fn{fn_idx},'argsout');
  catch errorstruct
    fprintf('%s: %s\n', errorstruct.identifier, errorstruct.message);
    for stack_idx = 1:length(errorstruct.stack)
      fprintf('  %s: %d\n', errorstruct.stack(stack_idx).name, errorstruct.stack(stack_idx).line);
    end
    save(task_out_fn{fn_idx},'argsout','errorstruct');
  end
end

return;
