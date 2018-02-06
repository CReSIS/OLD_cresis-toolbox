function cluster_job(task_in_fn_dir,task_out_fn_dir,job_list)
% cluster_job(task_in_fn_dir,task_out_fn_dir,job_list)
%
% This M-file should be compiled with cluster_compile.
%
% Then when the distributed toolbox calls the cluster_job.sh command via
% the scheduler which calls the compiled version of this function.
%
% cluster_submit_batch.m shows how to use the scheduler from Matlab
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

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
  % Special merge of argsin cell array
  if isfield(sparam.static_param,'argsin')
    sparam_argsin_numel = numel(sparam.static_param.argsin);
  else
    sparam.static_param.argsin = {};
    sparam_argsin_numel = 0;
  end
  if isfield(dparam.(dparam_task_field),'argsin')
    dparam_argsin_numel = numel(dparam.(dparam_task_field).argsin);
  else
    dparam.(dparam_task_field).argsin = {};
    dparam_argsin_numel = 0;
  end
  for idx = 1:max(sparam_argsin_numel,dparam_argsin_numel)
    if idx <= sparam_argsin_numel
      if idx <= dparam_argsin_numel
        param.argsin{idx} = merge_structs(sparam.static_param.argsin{idx},dparam.(dparam_task_field).argsin{idx});
      else
        param.argsin{idx} = sparam.static_param.argsin{idx};
      end
    else
      param.argsin{idx} = dparam.(dparam_task_field).argsin{idx};
    end
  end
  
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
    fprintf('%s: %s\n', mfilename, param.notes);
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
