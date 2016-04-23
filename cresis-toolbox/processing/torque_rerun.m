function ctrl = torque_rerun(ctrl,param)
% ctrl = torque_rerun(ctrl,param)
%
% Tries to rerun jobs that have failed up to the max specified.
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

if ~isfield(ctrl.sched,'submit_mode')
  ctrl.sched.submit_mode = 'group';
end

if ~isempty(ctrl.submission_queue)
  if strcmpi(ctrl.sched.submit_mode,'group')
    % Submit remaining jobs all at once
    [ctrl,new_job_id] = torque_submit_job(ctrl,ctrl.submission_queue);
    fprintf('Submitted these jobs in torque job %d:\n  %d', new_job_id, ctrl.submission_queue(1))
    if length(ctrl.submission_queue) > 1
      fprintf(', %d', ctrl.submission_queue(2:end));
    end
    fprintf('\n');
    ctrl.submission_queue = [];
  elseif strcmpi(ctrl.sched.submit_mode,'fill')
    submissions = ctrl.sched.num_submissions;
    jobs_per_submission = ceil(length(ctrl.submission_queue) / submissions);
    while ~isempty(ctrl.submission_queue)
      submission_queue = ctrl.submission_queue(1:min(jobs_per_submission,end));
      [ctrl,new_job_id] = torque_submit_job(ctrl,submission_queue);
      fprintf('Submitted these jobs in torque job %d/%d:\n  %d', ctrl.batch_id, new_job_id, submission_queue(1))
      if length(ctrl.submission_queue) > 1
        fprintf(', %d', submission_queue(2:end));
      end
      fprintf('\n');
      ctrl.submission_queue = ctrl.submission_queue(jobs_per_submission+1:end);
    end
  end
end

% Wait for jobs to complete
while ~all(ctrl.job_status == 'C' | ctrl.job_status == 'E')
  pause(3);
  if exist(fullfile(ctrl.batch_dir,'keyboard'), 'file')
    % Hold keyboard file exists
    keyboard
  end
  ctrl = torque_job_status(ctrl);
end
% Wait for jobs to complete phase 2 (no errors, but still exitting jobs)
while any(ctrl.job_status == 'E' & ctrl.error_mask == 0)
  warning(sprintf('Some jobs are still in the "Exiting" state. Consider waiting a few moments and then running:\nctrl = torque_job_status(ctrl);\nany(ctrl.job_status == ''E'' & ctrl.error_mask == 0) %% Should return zero when all jobs are done exiting\ndbcont\n'));
  pause(30);
  if exist(fullfile(ctrl.batch_dir,'keyboard'), 'file')
    % Hold keyboard file exists
    keyboard
  end
  ctrl = torque_job_status(ctrl);
end

retry = 0;
while retry < ctrl.sched.max_retries
  % Check errors for rerun
  if all(ctrl.error_mask == 0)
    return;
  end
  retry = retry + 1;
  
  if ctrl.sched.stop_on_fail
    warning('Errors occured, check them... if safe to proceed with retries then type dbcont');
    keyboard;
  end
  
  bad_job_ids = find(ctrl.error_mask);
  for job_id = reshape(bad_job_ids,[1 length(bad_job_ids)])
    ctrl = torque_resubmit_task(ctrl,job_id);
  end
  
  if ~isempty(ctrl.submission_queue)
    if strcmpi(ctrl.sched.submit_mode,'group')
      % Submit remaining jobs all at once
      [ctrl,new_job_id] = torque_submit_job(ctrl,ctrl.submission_queue);
      fprintf('Submitted these jobs in torque job %d:\n  %d', new_job_id, ctrl.submission_queue(1))
      if length(ctrl.submission_queue) > 1
        fprintf(', %d', ctrl.submission_queue(2:end));
      end
      fprintf('\n');
      ctrl.submission_queue = [];
    elseif strcmpi(ctrl.sched.submit_mode,'fill')
      submissions = ctrl.sched.num_submissions;
      jobs_per_submission = ceil(length(ctrl.submission_queue) / submissions);
      while ~isempty(ctrl.submission_queue)
        submission_queue = ctrl.submission_queue(1:min(jobs_per_submission,end));
        [ctrl,new_job_id] = torque_submit_job(ctrl,submission_queue);
        fprintf('Submitted these jobs in torque job %d:\n  %d', new_job_id, submission_queue(1))
        if length(ctrl.submission_queue) > 1
          fprintf(', %d', submission_queue(2:end));
        end
        fprintf('\n');
        ctrl.submission_queue = ctrl.submission_queue(jobs_per_submission+1:end);
      end
    end
  end
  
  % Wait for jobs to complete
  while ~all(ctrl.job_status == 'C' | ctrl.job_status == 'E')
    pause(3);
    if exist(fullfile(ctrl.batch_dir,'keyboard'), 'file')
      % Hold keyboard file exists
      keyboard
    end
    ctrl = torque_job_status(ctrl);
  end
  
end

return
