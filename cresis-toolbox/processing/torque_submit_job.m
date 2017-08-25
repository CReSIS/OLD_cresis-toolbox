function [ctrl,new_job_id] = torque_submit_job(ctrl,submission_queue)
% [ctrl,new_job_id] = torque_submit_job(ctrl,submission_queue)
%
% Creates a new job. Calls to this function need to be proceeded by
% a single call to torque_new_batch.m.
%
% Inputs:
% ctrl = ctrl structure returned from torque_new_batch
%  .sched = scheduler structure
%   .worker_fn = path to worker
%   .submit_arguments = submission arguments to add to qsub (-v currently
%     not supported)
%  .in_path_dir = input arguments directory
%  .out_path_dir = output arguments directory
%  .stdout_path_dir = standard output directory
%  .error_path_dir = error directory
%
% Outputs:
% ctrl = updated ctrl structure with new job
% job_id = ID of job (starts counting from one and never repeats)
%
% Author: John Paden
%
% See also: torque_batch_list torque_cleanup torque_compile ...
%   torque_create_task torque_hold torque_job_list torque_job_status ...
%   torque_new_batch torque_print torque_rerun

if ~isfield(ctrl.sched,'test_mode')
  ctrl.sched.test_mode = 0;
end
if ~isfield(ctrl.sched,'group_submit_arguments')
  ctrl.sched.group_submit_arguments = ctrl.sched.submit_arguments;
end
if ~isfield(ctrl.sched,'group_walltime')
  ctrl.sched.group_walltime = 1;
end
if ~isfield(ctrl.sched,'submit_pause')
  ctrl.sched.submit_pause = 2;
end
if ~isfield(ctrl.sched,'interactive')
  ctrl.sched.interactive = 0;
end

%% Create the temporary file names
in_path = ctrl.in_path_dir;
out_path = ctrl.out_path_dir;
job_id = submission_queue(end); % Use last job ID for the stdout and error files
stdout_path = fullfile(ctrl.stdout_path_dir,sprintf('stdout_%d.txt',job_id));
error_path = fullfile(ctrl.error_path_dir,sprintf('error_%d.txt',job_id));

new_job_status = 'Q';

%% Create QSUB system command
worker = ctrl.sched.worker_fn;
[tmp worker_name] = fileparts(worker);

job_list_str = sprintf('%dd',submission_queue); job_list_str = job_list_str(1:end-1);
submit_arguments = sprintf(ctrl.sched.group_submit_arguments,length(submission_queue) * ctrl.sched.group_walltime);
% Add "qsub -m abe -M your@email.edu" to debug:
cmd = sprintf('qsub %s -e %s -o %s -v INPUT_PATH="%s",OUTPUT_PATH="%s",CUSTOM_TORQUE="1",JOB_LIST=''%s'' %s  </dev/null', ...
  submit_arguments, error_path, stdout_path, in_path, out_path, job_list_str, worker);

%% Run the qsub command
if ctrl.sched.test_mode
  % Create a fake job id since we are in test mode
  new_job_id = 1000000 + job_id;
  % Create all the output files for these jobs since we are in test mode
  for job_idx = 1:length(submission_queue)
    job_id = submission_queue(job_idx);
    out_fn = fullfile(out_path,sprintf('out_%d.mat',job_id));
    argsout = {1};
    save(out_fn,'argsout');
  end
  
elseif ctrl.sched.interactive
  cmd = sprintf('qsub -I %s -e %s -o %s -v INPUT_PATH="%s",OUTPUT_PATH="%s",CUSTOM_TORQUE="1",JOB_LIST=''%s''', ...
    submit_arguments, error_path, stdout_path, in_path, out_path, job_list_str);
  fprintf('1. Run the command from the bash shell:\n  %s\n', cmd);
  fprintf('2. Once the interactive mode starts, run the command in the interactive shell:  %s\n', worker);
  fprintf('3. Once the job completes, exit the interactive shell which causes torque to realize the job is complete.\n');
  fprintf('4. In Matlab, set new_job_id to the torque job ID that you get from qsub. For example "2466505.m2" would need to have "new_job_id = 2466505" run.\n');
  fprintf('5. Once the job finishes, run "dbcont" in Matlab.\n');
  keyboard
else
  status = -1;
  torque_attempts = 0;
  while status ~= 0
    try
      [status,result] = system(cmd);
      pause(ctrl.sched.submit_pause);
    catch
      cmd
      warning('system call failed');
      keyboard;
    end
    if status ~= 0
      warning('qsub failed %d %s', status, result);
      torque_attempts = torque_attempts + 1;
      delay_period = 3*2^(torque_attempts-1);
      fprintf('  Delaying %d seconds\n', delay_period)
      pause(delay_period);
      if torque_attempts > 10
        % There is potentially something wrong with the torque scheduler if
        % job submission fails this many times. Look into it before running
        % "dbcont" to keep trying to submit.
        keyboard;
      end
    else
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
  end
end

% Update job IDs in job ID file
fid = fopen(ctrl.job_id_file,'r+');
for job_idx = 1:length(submission_queue)
  job_id = submission_queue(job_idx);
  fseek(fid, 21*(job_id-1), -1);
  fprintf(fid,'%-20d\n', new_job_id);
  ctrl.job_id_list(job_id) = new_job_id;
  ctrl.job_status(job_id) = 'Q';
  ctrl.error_mask(job_id) = 0;
end
fclose(fid);

end
