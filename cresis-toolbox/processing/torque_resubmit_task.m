function ctrl = torque_resubmit_task(ctrl,job_id)
% ctrl = torque_resubmit_task(ctrl,job_id)
%
% Resubmit a task.
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
% job_id = ID of job (starts counting from one and never repeats)
%
% Outputs:
% ctrl = updated ctrl structure with new job
%
% Author: John Paden
%
% See also: torque_batch_list torque_cleanup torque_compile ...
%   torque_create_task torque_hold torque_job_list torque_job_status ...
%   torque_new_batch torque_print torque_rerun

if ~isfield(ctrl.sched,'submit_mode')
  ctrl.sched.submit_mode = 'group';
end
if ~isfield(ctrl.sched,'group_size')
  ctrl.sched.group_size = 1;
end

%% Check number of processes in queue
if ~isfield(ctrl,'num_jobs_in_queue')
  ctrl = torque_job_status(ctrl);
  ctrl.num_jobs_in_queue = sum(ctrl.job_status == 'Q');
end
while ctrl.num_jobs_in_queue >= ctrl.sched.max_in_queue
  pause(1);
  ctrl = torque_job_status(ctrl);
  ctrl.num_jobs_in_queue = sum(ctrl.job_status == 'Q');
  if exist(fullfile(ctrl.batch_dir,'keyboard'), 'file')
    keyboard
  end
end
ctrl.num_jobs_in_queue = ctrl.num_jobs_in_queue + 1;

in_path = fullfile(ctrl.in_path_dir,sprintf('in_%d.mat',job_id));

%% Not actually submitting the job right now... putting it into the submission queue
ctrl.submission_queue = cat(2,ctrl.submission_queue,job_id);
new_job_status = 'T';
new_job_id = -1;
%% Get the torque job id for this job and store that
ctrl.job_id_list(job_id) = new_job_id;
ctrl.job_status(job_id) = new_job_status;
ctrl.error_mask(job_id) = 0;
% Write new job ID to job_id file
fid = fopen(ctrl.job_id_file,'r+');
fseek(fid, 21*(job_id-1), -1);
fprintf(fid,'%-20d\n', new_job_id);
fclose(fid);

%% Submit this job to torque if submission queue is full enough
if strcmpi(ctrl.sched.submit_mode,'group') && length(ctrl.submission_queue) >= ctrl.sched.group_size
  ctrl = torque_submit_job(ctrl,ctrl.submission_queue);
  ctrl.submission_queue = [];
end

load(in_path,'param');
if ~isempty(param.notes)
  fprintf('  rerun %d:%d/%d: %s (%s)\n', ctrl.batch_id, job_id, ctrl.job_id_list(job_id), param.notes, datestr(now));
end

end
