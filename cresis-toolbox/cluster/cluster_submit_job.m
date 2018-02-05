function [ctrl,new_job_id] = cluster_submit_job(ctrl,job_tasks,job_cpu_time,job_mem)
% [ctrl,new_job_id] = cluster_submit_job(ctrl,job_tasks,job_cpu_time,job_mem)
%
% Creates a new job. Calls to this function need to be proceeded by
% a single call to cluster_new_batch.m.
%
% Inputs:
% ctrl = ctrl structure returned from cluster_new_batch
%  .cluster = cluster parameters structure
%   .worker_fn = path to worker
%   .submit_arguments = submission arguments to add to qsub (-v currently
%     not supported)
%  .in_fn_dir = input arguments directory
%  .out_fn_dir = output arguments directory
%  .stdout_fn_dir = standard output directory
%  .error_fn_dir = error directory
%
% Outputs:
% ctrl = updated ctrl structure with new job
% job_id = ID of job (starts counting from one and never repeats)
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

if ~isfield(ctrl.cluster,'submit_pause')
  ctrl.cluster.submit_pause = 2;
end
if ~isfield(ctrl.cluster,'interactive')
  ctrl.cluster.interactive = 0;
end

%% Create the temporary file names
in_fn = ctrl.in_fn_dir;
out_fn = ctrl.out_fn_dir;
job_id = job_tasks(end); % Use last job ID for the stdout and error files
stdout_fn = fullfile(ctrl.stdout_fn_dir,sprintf('stdout_%d.txt',job_id));
error_fn = fullfile(ctrl.error_fn_dir,sprintf('error_%d.txt',job_id));

new_job_status = 'Q';

task_list_str = sprintf('%dd',job_tasks); task_list_str = task_list_str(1:end-1);

%% Submit the job
% =========================================================================
if strcmpi(ctrl.cluster.type,'torque')
  %% Run the qsub command
  
  worker = ctrl.cluster.worker_fn;
  [tmp worker_name] = fileparts(worker);
  
  submit_arguments = sprintf(ctrl.cluster.group_submit_arguments,length(submission_queue) * ctrl.cluster.group_walltime);
  
  % Add "qsub -m abe -M your@email.edu" to debug:
  if ctrl.cluster.interactive
    cmd = sprintf('qsub -I %s -e %s -o %s -v INPUT_PATH="%s",OUTPUT_PATH="%s",CUSTOM_TORQUE="1",JOB_LIST=''%s''', ...
      submit_arguments, error_fn, stdout_fn, in_fn, out_fn, job_list_str);
    fprintf('1. Run the command from the bash shell:\n  %s\n', cmd);
    fprintf('2. Once the interactive mode starts, run the command in the interactive shell:  %s\n', worker);
    fprintf('3. Once the job completes, exit the interactive shell which causes torque to realize the job is complete.\n');
    fprintf('4. In Matlab, set new_job_id to the torque job ID that you get from qsub. For example "2466505.m2" would need to have "new_job_id = 2466505" run.\n');
    fprintf('5. Once the job finishes, run "dbcont" in Matlab.\n');
    while ~exist('new_job_id','var')
      keyboard
    end
  else
    cmd = sprintf('qsub %s -e %s -o %s -v INPUT_PATH="%s",OUTPUT_PATH="%s",CUSTOM_TORQUE="1",JOB_LIST=''%s'' %s  </dev/null', ...
      submit_arguments, error_fn, stdout_fn, in_fn, out_fn, job_list_str, worker);
    [status,result] = robust_system(cmd);
    
    [job_id_str,result_tok] = strtok(result,'.');
    try
      new_job_id = str2double(job_id_str);
      if isnan(new_job_id)
        job_id_str
        warning('job_id_str expected numeric, but is not');
        keyboard;
      end
    catch
      job_id_str
      warning('job_id_str expected numeric, but is not');
      keyboard;
    end
  end
  ctrl.active_jobs = ctrl.active_jobs + 1;
  
elseif strcmpi(ctrl.cluster.type,'matlab')
  %% Create the job on the matlab cluster/job manager
  new_job = createJob(ctrl.cluster.jm);
  new_job_id = new_job.ID;
  task = createTask(new_job,@cluster_job,0,{in_fn,out_fn,task_list_str});
  ctrl.active_jobs = ctrl.active_jobs + 1;
  while ~strcmpi(task.State,'pending')
    warning('%s: pausing because task is not in pending state.', mfilename);
    pause(1);
  end
  submit(new_job);
  
elseif strcmpi(ctrl.cluster.type,'debug')
  %% Run the command now
  new_job_id = job_tasks(end);
  
  if ~isfield(ctrl.cluster,'run_mode')
    ctrl.cluster.run_mode = [];
  end
  cluster_exec_job(ctrl,job_tasks,ctrl.cluster.run_mode);
end

%% Update job IDs in job ID file
fid = fopen(ctrl.job_id_fn,'r+');
for task_idx = 1:length(job_tasks)
  task_id = job_tasks(task_idx);
  fseek(fid, 21*(task_id-1), -1);
  fprintf(fid,'%-20d\n', new_job_id);
  ctrl.job_id_list(task_id) = new_job_id;
  ctrl.job_status(task_id) = 'Q';
  ctrl.error_mask(task_id) = 0;
end
fclose(fid);

end
