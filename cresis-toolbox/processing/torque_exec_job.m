function torque_exec_job(ctrl,job_ids,run_mode)
% torque_exec_job(ctrl,job_id,run_mode)
%
% Re-executes a specific job from the command line (runs it locally).
%
% Inputs:
% ctrl = ctrl structure returned from torque_new_batch
%  .sched = scheduler structure
%   .worker_fn = path to worker
%  .in_path_dir = input arguments directory
%  .out_path_dir = output arguments directory
% job_id = specific job ID, or vector of job IDs, to run
% run_mode = optional scalar integer indicating how to run the job:
%   1: Run job as if you were in no scheduler mode [default]
%   2: Run job through uncompiled worker_task function
%   3: Run job through compiled worker_task function
%
% Author: John Paden
%
% See also: torque_batch_list torque_cleanup torque_compile ...
%   torque_create_task torque_hold torque_job_list torque_job_status ...
%   torque_new_batch torque_print torque_rerun

if ~exist('run_mode','var') || isempty(run_mode)
  run_mode = 1;
end

for job_idx = 1:length(job_ids)
  job_id = job_ids(job_idx);
  
  task_in_fn = fullfile(ctrl.in_path_dir,sprintf('in_%d.mat',job_id));
  task_out_fn = fullfile(ctrl.out_path_dir,sprintf('out_%d.mat',job_id));
  
  if run_mode == 1
    fprintf('Loading input arguments and control parameters\n');
    task = load(task_in_fn);
    
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
    
    argsout = {};
    fprintf('Eval %s\n', eval_cmd);
    eval(eval_cmd);
    fprintf('Done eval\n');
    save(task_out_fn,'argsout');
    
  elseif run_mode == 2
    setenv('INPUT_PATH',ctrl.in_path_dir);
    setenv('OUTPUT_PATH',ctrl.out_path_dir);
    job_list_str = sprintf('%dd',job_id); job_list_str = job_list_str(1:end-1);
    setenv('JOB_LIST',job_list_str);
    setenv('CUSTOM_TORQUE','1');
    cur_dir = cd;
    cd(fullfile(fileparts(ctrl.sched.worker_fn)));
    worker_task;
    cd(cur_dir);
    
  elseif run_mode == 3
    setenv('INPUT_PATH',ctrl.in_path_dir);
    setenv('OUTPUT_PATH',ctrl.out_path_dir);
    job_list_str = sprintf('%dd',job_id); job_list_str = job_list_str(1:end-1);
    setenv('JOB_LIST',job_list_str);
    setenv('CUSTOM_TORQUE','1');
    system(ctrl.sched.worker_fn);
  end
end
return;
