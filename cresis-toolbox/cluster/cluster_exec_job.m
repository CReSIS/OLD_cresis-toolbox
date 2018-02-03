function cluster_exec_job(ctrl,task_ids,run_mode)
% cluster_exec_job(ctrl,task_ids,run_mode)
%
% Re-executes a specific job from the command line (runs it locally).
%
% Inputs:
% ctrl = ctrl structure returned from cluster_new_batch
%  .sched = scheduler structure
%   .worker_fn = path to worker
%  .in_fn_dir = input arguments directory
%  .out_fn_dir = output arguments directory
% task_ids = vector of task IDs, to run
% run_mode = optional scalar integer indicating how to run the job:
%   1: Run job as if you were in no scheduler mode [default]
%   2: Run job through uncompiled worker_task function
%   3: Run job through compiled worker_task function
%
% Author: John Paden
%
% See also: cluster_batch_list cluster_cleanup cluster_compile ...
%   cluster_create_task cluster_hold cluster_job_list cluster_job_status ...
%   cluster_new_batch cluster_print cluster_rerun

if ~exist('run_mode','var') || isempty(run_mode)
  run_mode = 1;
end

for task_idx = 1:length(task_ids)
  task_id = task_ids(task_idx);
  
  if run_mode == 1
    fprintf('  %s: Loading input arguments and control parameters\n', mfilename);

    % Create in/out filenames
    static_in_fn = fullfile(ctrl.in_fn_dir,'static.mat');
    dynamic_in_fn = fullfile(ctrl.in_fn_dir,'dynamic.mat');
    out_fn = fullfile(ctrl.out_fn_dir,sprintf('out_%d.mat',task_id));
    
    % Read input
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
    
    argsout = {};
    errorstruct = [];
    fprintf('  %s: Eval %s\n', mfilename, eval_cmd);
    try
      eval(eval_cmd);
    catch errorstruct
    end
    fprintf('  %s: Done eval\n', mfilename);
    save(out_fn,'argsout','errorstruct');
    
  elseif run_mode == 2
    setenv('INPUT_PATH',ctrl.in_fn_dir);
    setenv('OUTPUT_PATH',ctrl.out_fn_dir);
    job_list_str = sprintf('%dd',task_id); job_list_str = job_list_str(1:end-1);
    setenv('JOB_LIST',job_list_str);
    setenv('CUSTOM_CLUSTER','1');
    cluster_job;
    
  elseif run_mode == 3
    setenv('INPUT_PATH',ctrl.in_fn_dir);
    setenv('OUTPUT_PATH',ctrl.out_fn_dir);
    job_list_str = sprintf('%dd',task_id); job_list_str = job_list_str(1:end-1);
    setenv('JOB_LIST',job_list_str);
    setenv('CUSTOM_CLUSTER','1');
    system(ctrl.sched.cluster_job_fn);
  end
end
return;
