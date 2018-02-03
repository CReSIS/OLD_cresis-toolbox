function ctrl = cluster_task_update(ctrl,task_id)
% ctrl = cluster_task_update(ctrl,task_id)
%
% Updates a specific task's status. Support function for cluster_job_status.
%
% Inputs:
% ctrl: ctrl structure returned from cluster_new_batch
%  .job_id_list = Nx1 vector of cluster job IDs
%  .job_status = Nx1 vector of job status
%  .out_fn_dir = string containing the output directory
%  .error_mask = Nx1 vector of error status
%  .out = up-to-Nx1 vector of cells containing the outputs for each
%    job as they complete (read in from the output.mat files)
% job_id: the job id to check up on
%
% Outputs:
% ctrl = updated ctrl structure
%
% C -  Job is completed after having run.
% E -  Job is exiting after having run.
% H -  Job is held.
% Q -  job is queued, eligible to run or routed.
% R -  job is running.
% T -  job is being moved to new location.
% W -  job is waiting for its execution time
%      (-a option) to be reached.
% S -  (Unicos only) job is suspend.
%
% Author: John Paden
%
% See also: cluster_batch_list cluster_cleanup cluster_compile ...
%   cluster_create_task cluster_hold cluster_job_list cluster_job_status ...
%   cluster_new_batch cluster_print cluster_run

% Create in/out filenames
static_in_fn = fullfile(ctrl.in_fn_dir,'static.mat');
dynamic_in_fn = fullfile(ctrl.in_fn_dir,'dynamic.mat');
out_fn = fullfile(ctrl.out_fn_dir,sprintf('out_%d.mat',task_id));

% Read input
sparam = load(static_in_fn);
dparam_task_field = sprintf('dparam_%d',task_id);
dparam = load(dynamic_in_fn,dparam_task_field);
param = merge_structs(sparam.static_param,dparam.(dparam_task_field));

ctrl.notes{task_id} = param.notes;
ctrl.cpu_time(task_id) = param.cpu_time;
ctrl.mem(task_id) = param.mem;
ctrl.success{task_id} = param.success;

out_fn_exist_error = 1;
out_fn_load_error = 2;
argsout_exist_error = 4;
argsout_length_error = 8;
errorstruct_exist_error = 16;
errorstruct_contains_error = 32;
success_eval_error = 64;

error_mask = 0;

out_fn = fullfile(ctrl.out_fn_dir,sprintf('out_%d.mat',task_id));
if exist(out_fn,'file')
  success = robust_cmd('out = load(out_fn);',3);
  if ~success
    out = [];
    error_mask = error_mask + out_fn_load_error;
  end
else
  out = [];
  error_mask = error_mask + out_fn_exist_error;
end

if isfield(out,'argsout')
  if length(out.argsout) ~= param.num_args_out
    error_mask = error_mask + argsout_length_error;
  end
else
  error_mask = error_mask + argsout_exist_error;
end

if isfield(out,'errorstruct')
  if ~isempty(out.errorstruct)
    error_mask = error_mask + errorstruct_contains_error;
  end
else
  error_mask = error_mask + errorstruct_exist_error;
end

try
  eval(param.success);
catch ME
  error_mask = error_mask + success_eval_error;
end

ctrl.error_mask(task_id) = 0;
if ctrl.job_status(task_id) == 'T'
  % This job was not found in the cluster: if output exists, then assume it
  % is 'C' and update error_mask based on output files. Otherwise, assume
  % it has not been submitted yet.
  if ~bitand(error_mask,out_fn_exist_error)
    ctrl.job_status(task_id) = 'C';
    ctrl.error_mask(task_id) = error_mask;
  end
  
elseif ctrl.job_status(task_id) == 'C'
  % This job was found on the cluster and was completed: update error
  % flag based on output files. Make sure this job is not in the submission
  % queue.
  ctrl.submission_queue = ctrl.submission_queue(ctrl.submission_queue~=task_id);
  ctrl.error_mask(task_id) = error_mask;
  
else
  % This job was found in the cluster: if output exists, then assume it
  % is 'C' and update error_mask based on output files. Make sure this job
  % is not in the submission queue.
  ctrl.submission_queue = ctrl.submission_queue(ctrl.submission_queue~=task_id);
  if ~bitand(error_mask,out_fn_exist_error)
    ctrl.job_status(task_id) = 'C';
    ctrl.error_mask(task_id) = error_mask;
  end
  
end

if ctrl.error_mask(task_id)
  warning(' Job Error %d:%d/%d\n', ctrl.batch_id, task_id, ctrl.job_id_list(task_id));
  if bitand(ctrl.error_mask(task_id),out_fn_exist_error)
    fprintf('  Output file does not exist\n');
  end
  if bitand(ctrl.error_mask(task_id),out_fn_load_error)
    fprintf('  Output file exists, but failed to load\n');
  end
  if bitand(ctrl.error_mask(task_id),argsout_exist_error)
    fprintf('  argsout does not exist in output file\n');
  end
  if bitand(ctrl.error_mask(task_id),argsout_length_error)
    fprintf('  argsout is the wrong length in the output file\n');
  end
  if bitand(ctrl.error_mask(task_id),errorstruct_exist_error)
    fprintf('  errorstruct does not exist in output file\n');
  end
  if bitand(ctrl.error_mask(task_id),errorstruct_contains_error)
    fprintf('  errorstruct contains an error:\n');
    warning('%s',out.errorstruct.getReport);
  end
end

if ctrl.job_status(task_id) == 'C' && ctrl.error_mask(task_id)
  % Job is completed and has an error
  
  if ctrl.retries(task_id) < ctrl.cluster.max_retries
    % Update task to ctrl structure
    ctrl.submission_queue = cat(2,ctrl.submission_queue,task_id);
    new_job_status = 'T';
    new_job_id = -1;
    ctrl.job_id_list(task_id) = new_job_id;
    ctrl.job_status(task_id) = new_job_status;
    ctrl.error_mask(task_id) = 0;
    ctrl.retries(task_id) = ctrl.retries(task_id) + 1;
    
    % Update job IDs in job ID file
    fid = fopen(ctrl.job_id_fn,'r+');
    fseek(fid, 21*(task_id-1), -1);
    fprintf(fid,'%-20d\n', new_job_id);
    fclose(fid);
    
    % Print out retry message
    fprintf(' Retry %d Job %d:%d/%d %s (%s)\n', ctrl.retries(task_id), ctrl.batch_id, task_id, ctrl.job_id_list(task_id), param.notes, datestr(now));
  end

end
