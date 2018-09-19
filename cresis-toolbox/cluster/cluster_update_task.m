function ctrl = cluster_update_task(ctrl,task_id,update_mode)
% ctrl = cluster_update_task(ctrl,task_id,update_mode)
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
  fprintf('This batch has a hold. Run <strong>cluster_hold(ctrl)</strong> to remove. Either way, run <strong>dbcont</strong> to continue.\n');
  keyboard
end

% Get the job ID
job_id = ctrl.job_id_list(task_id);

% Get the task ID that the job ID is writing to
task_id_out = find(ctrl.job_id_list == job_id,1,'last');

if any(strcmpi(ctrl.cluster.type,{'torque','slurm'}))
  % Extract information from stdout
  stdout_fn = fullfile(ctrl.stdout_fn_dir,sprintf('stdout_%d.txt',task_id_out));
  last_task_id = -1;
  if exist(stdout_fn,'file')
    if exist(stdout_fn,'file')
      stdout_file_str = '';
      try
        fid = fopen(stdout_fn);
        stdout_file_str = fread(fid,inf,'char=>char');
        stdout_file_str = stdout_file_str(:).';
        fclose(fid);
      end

      % Find the hostname of the execution node
      try
        idx = regexp(stdout_file_str,'hostname:');
        end_idx = idx+9 + find(stdout_file_str(idx+9:end)==' ',1)-1;
        hostname = stdout_file_str(idx+9:end_idx-1);
      catch
        hostname = '';
      end
      
      % Find the number of attempts to start the job
      try
        idx = regexp(stdout_file_str,'attempt:');
        attempt = sscanf(stdout_file_str(idx+8:end),'%d');
      catch
        attempt = -1;
      end
      if isempty(attempt)
        attempt = -1;
      end
      
      % Find the allowed maximum number of attempts to start the job
      try
        idx = regexp(stdout_file_str,'max_attempts:');
        max_attempts = sscanf(stdout_file_str(idx+13:end),'%d');
      catch
        max_attempts = -1;
      end
      if isempty(max_attempts)
        max_attempts = -1;
      end
      
      % Find the last task that this job started
      try
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

if update_mode && strcmpi(ctrl.cluster.type,'matlab') && ctrl.job_status(task_id) == 'C'
  stdout_fn = fullfile(ctrl.stdout_fn_dir,sprintf('stdout_%d.txt',task_id_out));
  error_fn = fullfile(ctrl.error_fn_dir,sprintf('error_%d.txt',task_id_out));
  
  if ~exist(error_fn,'file')
    fid = fopen(error_fn,'w');
    str = evalc(sprintf('ctrl.cluster.jm.Jobs(%d).Tasks(1)',job_id));
    fwrite(fid,str,'char');
    fclose(fid);
  end
  
  if ~exist(stdout_fn,'file')
    fid = fopen(stdout_fn,'w');
    fwrite(fid,ctrl.cluster.jm.Jobs(job_id).Tasks(1).Diary,'char');
    fclose(fid);
  end
end

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
  warning(' Job Error %d:%d/%d (lead task %d)\n', ctrl.batch_id, task_id, job_id, task_id_out);
  if any(strcmpi(ctrl.cluster.type,{'torque','slurm'}))
    fprintf('   hostname:%s attempt:%d max_attempts:%d\n', hostname, attempt, max_attempts);
  end
  if bitand(ctrl.error_mask(task_id),out_fn_exist_error)
    fprintf('  Output file does not exist: %s\n', out_fn);
  end
  if bitand(ctrl.error_mask(task_id),out_fn_load_error)
    fprintf('  Output file exists, but failed to load: %s\n', out_fn);
  end
  if bitand(ctrl.error_mask(task_id),argsout_exist_error)
    fprintf('  argsout does not exist in output file: %s\n', out_fn);
  end
  if bitand(ctrl.error_mask(task_id),argsout_length_error)
    fprintf('  argsout is the wrong length in the output file: %s\n', out_fn);
  end
  if bitand(ctrl.error_mask(task_id),errorstruct_exist_error)
    fprintf('  errorstruct does not exist in output file: %s\n', out_fn);
  end
  if bitand(ctrl.error_mask(task_id),errorstruct_contains_error)
    fprintf('  errorstruct contains an error:\n');
    warning('%s',out.errorstruct.getReport);
    if ctrl.cluster.stop_on_error
      keyboard
    end
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

if update_mode
  if ctrl.cluster.cpu_time_mult*ctrl.cpu_time(task_id)*0.9 < ctrl.cpu_time_actual(task_id)
    warning(' %d:%d/%d: CPU time actual (%.0f sec) is more than 90%% of estimated time (%.0f sec). Consider revising estimates.', ...
      ctrl.batch_id, task_id, job_id, ctrl.cpu_time_actual(task_id), ctrl.cluster.cpu_time_mult*ctrl.cpu_time(task_id));
  elseif ctrl.cpu_time_actual(task_id)>0 && ctrl.cluster.cpu_time_mult*ctrl.cpu_time(task_id)*0.3 > ctrl.cpu_time_actual(task_id)
    warning(' %d:%d/%d: CPU time actual (%.0f sec) is less than 30%% of estimated time (%.0f sec). Consider revising estimates.', ...
      ctrl.batch_id, task_id, job_id, ctrl.cpu_time_actual(task_id), ctrl.cluster.cpu_time_mult*ctrl.cpu_time(task_id));
  end
end

if update_mode && ctrl.job_status(task_id) == 'C' && ctrl.error_mask(task_id)
  % Job is completed and has an error
  
  % Copy stdout and error files
  if any(strcmpi(ctrl.cluster.type,{'torque','slurm'}))
    stdout_fn = fullfile(ctrl.stdout_fn_dir,sprintf('stdout_%d.txt',task_id_out));
    error_fn = fullfile(ctrl.error_fn_dir,sprintf('error_%d.txt',task_id_out));
    retry_stdout_fn = fullfile(ctrl.stdout_fn_dir,sprintf('stdout_%d_%d.txt',task_id, ctrl.retries(task_id)));
    retry_error_fn = fullfile(ctrl.error_fn_dir,sprintf('error_%d_%d.txt',task_id, ctrl.retries(task_id)));
    if exist(stdout_fn,'file')
      copyfile(stdout_fn,retry_stdout_fn);
    end
    if exist(error_fn,'file')
      copyfile(error_fn,retry_error_fn);
    end
  end
  
  if ctrl.retries(task_id) < ctrl.cluster.max_retries
    if ~bitand(ctrl.error_mask(task_id),out_fn_exist_error)
      delete(out_fn);
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
    
    % If there are no more tasks using the output files, then delete them
    % so that the new retry files start with a blank file.
    if any(strcmpi(ctrl.cluster.type,{'torque','slurm'}))
      if ~any(ctrl.job_id_list == job_id)
        if exist(stdout_fn,'file')
          delete(stdout_fn);
        end
        if exist(error_fn,'file')
          delete(error_fn);
        end
      end
    end
    
    % Print out retry message
    fprintf(' Retry %d Job %d:%d/%d %s (%s)\n', ctrl.retries(task_id), ctrl.batch_id, task_id, ctrl.job_id_list(task_id), param.notes, datestr(now));
  else
    fprintf(' Out of retries %d Job %d:%d/%d %s (%s)\n', ctrl.retries(task_id), ctrl.batch_id, task_id, ctrl.job_id_list(task_id), param.notes, datestr(now));
  end

end
