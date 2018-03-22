function ctrl = cluster_task_update(ctrl,task_id,update_mode)
% ctrl = cluster_task_update(ctrl,task_id,update_mode)
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
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

if ~exist('update_mode','var') || isempty(update_mode)
  update_mode = 1;
end

% Create in/out filenames
static_in_fn = fullfile(ctrl.in_fn_dir,'static.mat');
dynamic_in_fn = fullfile(ctrl.in_fn_dir,'dynamic.mat');
out_fn = fullfile(ctrl.out_fn_dir,sprintf('out_%d.mat',task_id));

% Read input
sparam = load(static_in_fn);
if ~isfield(ctrl,'dparam') || numel(ctrl.dparam) < task_id || isempty(ctrl.dparam{task_id})
  tmp = load(dynamic_in_fn);
  ctrl.dparam = tmp.dparam;
end
param = merge_structs(sparam.static_param,ctrl.dparam{task_id});

ctrl.notes{task_id} = param.notes;
if ctrl.cpu_time(task_id) == 0
  ctrl.cpu_time(task_id) = param.cpu_time;
end
if ctrl.mem(task_id) == 0
  ctrl.mem(task_id) = param.mem;
end
ctrl.success{task_id} = param.success;

out_fn_exist_error = 1;
out_fn_load_error = 2;
argsout_exist_error = 4;
argsout_length_error = 8;
errorstruct_exist_error = 16;
errorstruct_contains_error = 32;
success_error = 64;
cluster_killed_error = 128;
walltime_exceeded_error = 256;
success_eval_error = 512;

error_mask = 0;

if update_mode && exist(ctrl.hold_fn,'file')
  fprintf('This batch has a hold. Run cluster_hold(ctrl) to remove. Either way, run dbcont to continue.\n');
  keyboard
end

% Get the job ID
job_id = ctrl.job_id_list(task_id);

% Get the task ID that the job ID is writing to
task_id_out = find(ctrl.job_id_list == job_id,1,'last');

if strcmpi(ctrl.cluster.type,'torque')
  % Look at stdout to determine the last task that this job started
  stdout_fn = fullfile(ctrl.stdout_fn_dir,sprintf('stdout_%d.txt',task_id_out));
  last_task_id = -1;
  if exist(stdout_fn,'file')
    if exist(stdout_fn,'file')
      try
        fid = fopen(stdout_fn);
        stdout_file_str = fread(fid,inf,'char=>char');
        stdout_file_str = stdout_file_str(:).';
        fclose(fid);
        % Find all task start messages in this job's stdout
        search_str = 'cluster_job: Load task ';
        idxs = regexp(stdout_file_str,search_str);
        if ~isempty(idxs)
          % Look at the task ID from the last start message
          last_task_id = sscanf(stdout_file_str(idxs(end)+length(search_str):end),'%d');
        end
      end
    end
  end
  if isempty(last_task_id)
    last_task_id = -1;
  end
  
  % If the last task started is the task that is getting updated, then
  % check to see if that task was killed during execution.
  error_fn = fullfile(ctrl.error_fn_dir,sprintf('error_%d.txt',task_id_out));
  if last_task_id == task_id && exist(error_fn,'file')
    try
      fid = fopen(error_fn);
      error_file_str = fread(fid,inf,'char=>char');
      error_file_str = error_file_str(:).';
      fclose(fid);
      if ~isempty(regexp(error_file_str,'PBS: job killed:'))
        error_mask = bitor(error_mask,cluster_killed_error);
      end
      if ~isempty(regexp(error_file_str,'PBS: job killed: walltime'))
        error_mask = bitor(error_mask,walltime_exceeded_error);
      end
    end
  end
end

if any(strcmpi(ctrl.cluster.type,{'torque','matlab','slurm'}))
  if ~bitand(error_mask,cluster_killed_error)
    out_fn = fullfile(ctrl.out_fn_dir,sprintf('out_%d.mat',task_id));
    % Sometimes the file system/matlab are slow in recognizing a file
    if update_mode && ctrl.job_status(task_id) == 'C' && ~exist(out_fn,'file')
      start_time = tic;
      while ~exist(out_fn,'file') && toc(start_time) < ctrl.cluster.file_check_pause;
        pause(5);
      end
      if ~exist(out_fn,'file')
        warning('Cluster batch %d task %d (%d) completed without producing out:\n  %s.', ctrl.batch_id, task_id, job_id, out_fn);
        %keyboard
      end
    end
  end
end

out_fn = fullfile(ctrl.out_fn_dir,sprintf('out_%d.mat',task_id));
if exist(out_fn,'file')
  success = robust_cmd('out = load(out_fn);',2);
  if ~success
    out = [];
    error_mask = bitor(error_mask,out_fn_load_error);
  end
else
  out = [];
  error_mask = bitor(error_mask,out_fn_exist_error);
end

if isfield(out,'argsout')
  if length(out.argsout) ~= param.num_args_out
    error_mask = bitor(error_mask,argsout_length_error);
  end
else
  error_mask = bitor(error_mask,argsout_exist_error);
end

if isfield(out,'errorstruct')
  if ~isempty(out.errorstruct)
    error_mask = bitor(error_mask,errorstruct_contains_error);
  end
else
  error_mask = bitor(error_mask,errorstruct_exist_error);
end

if isfield(out,'cpu_time_actual')
  ctrl.cpu_time_actual(task_id) = out.cpu_time_actual;
else
  ctrl.cpu_time_actual(task_id) = -1;
end

try
  eval(param.success); % Runs some form of "error_mask = bitor(error_mask,success_error);" on failure
catch success_eval_ME
  error_mask = bitor(error_mask,success_eval_error);
end

ctrl.error_mask(task_id) = 0;
if ctrl.job_status(task_id) == 'T'
  % This job was not found in the cluster: if output exists, then assume it
  % is 'C' and update error_mask based on output files. Otherwise, assume
  % it has not been submitted yet.
  if ~bitand(error_mask,out_fn_exist_error)
    ctrl.submission_queue = ctrl.submission_queue(ctrl.submission_queue~=task_id);
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

if update_mode && ctrl.error_mask(task_id)
  warning(' Job Error %d:%d/%d\n', ctrl.batch_id, task_id, job_id);
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
  if bitand(ctrl.error_mask(task_id),success_error)
    fprintf('  Task did not pass success criteria\n');
  end
  if bitand(ctrl.error_mask(task_id),success_eval_error)
    fprintf('  Task success condition failed to evaluate: %s\n', success_eval_ME.getReport);
  end
  if bitand(ctrl.error_mask(task_id),cluster_killed_error)
    fprintf('  Cluster killed this job\n');
  end
  if bitand(ctrl.error_mask(task_id),walltime_exceeded_error)
    fprintf('  Cluster killed this job due to wall time\n');
  end
end

if ctrl.cluster.cpu_time_mult*ctrl.cpu_time(task_id)*0.9 < ctrl.cpu_time_actual(task_id)
  warning(' %d:%d/%d: CPU time actual (%.0f sec) is more than 90%% of estimated time (%.0f sec). Consider revising estimates.', ...
    ctrl.batch_id, task_id, job_id, ctrl.cpu_time_actual(task_id), ctrl.cluster.cpu_time_mult*ctrl.cpu_time(task_id));
elseif ctrl.cpu_time_actual(task_id)>0 && ctrl.cluster.cpu_time_mult*ctrl.cpu_time(task_id)*0.3 > ctrl.cpu_time_actual(task_id)
  warning(' %d:%d/%d: CPU time actual (%.0f sec) is less than 30%% of estimated time (%.0f sec). Consider revising estimates.', ...
    ctrl.batch_id, task_id, job_id, ctrl.cpu_time_actual(task_id), ctrl.cluster.cpu_time_mult*ctrl.cpu_time(task_id));
end

if update_mode && ctrl.job_status(task_id) == 'C' && ctrl.error_mask(task_id)
  % Job is completed and has an error
  
  if ctrl.retries(task_id) < ctrl.cluster.max_retries
    if ~bitand(ctrl.error_mask(task_id),out_fn_exist_error)
      delete(out_fn);
    end

    % Move stdout and error files
    if any(strcmpi(ctrl.cluster.type,{'torque','slurm'}))
      stdout_fn = fullfile(ctrl.stdout_fn_dir,sprintf('stdout_%d.txt',task_id_out));
      error_fn = fullfile(ctrl.error_fn_dir,sprintf('error_%d.txt',task_id_out));
      attempt_stdout_fn = fullfile(ctrl.stdout_fn_dir,sprintf('stdout_%d_%d.txt',task_id_out, ctrl.retries(task_id)));
      attempt_error_fn = fullfile(ctrl.error_fn_dir,sprintf('error_%d_%d.txt',task_id_out, ctrl.retries(task_id)));
      if exist(stdout_fn,'file')
        movefile(stdout_fn,attempt_stdout_fn);
      end
      if exist(error_fn,'file')
        movefile(error_fn,attempt_error_fn);
      end
    end
    
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
