function [ctrl,job_id] = torque_create_task(ctrl,taskfunction,num_args_out,argsin,param)
% [ctrl,job_id] = torque_create_task(ctrl,taskfunction,num_args_out,argsin,param)
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
% task_function = function handle of job, this function handle tells
%   worker_task.m what to run
% num_args_out = number of output arguments to expect
% argsin = cell vector of input arguments
% param = structure of task creation parameters
%  .conforming: the function to be called is a custom torque conforming
%    function and its first return parameter returns a scalar value
%    indicating success when it is a 1 and failure otherwise
%    default is false
%  .notes: optional note to print after successful submission of job
%    default is empty (nothing is written out)
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

if ~exist('param','var') || isempty(param)
  param.conforming  = false;
  param.notes = '';
end
if ~isfield(param,'conforming')
  param.conforming = false;
end
if ~isfield(param,'notes')
  param.notes = '';
end
if ~isfield(ctrl.sched,'submit_mode')
  ctrl.sched.submit_mode = 'group';
end
if ~isfield(ctrl.sched,'group_size')
  ctrl.sched.group_size = 1;
end

%% Check for hold on this batch
if exist(fullfile(ctrl.batch_dir,'keyboard'), 'file')
  % Hold keyboard file exists
  keyboard
end

%% Get the new job_id
if ~isfield(ctrl,'job_id')
  [fid,msg] = fopen(ctrl.job_id_file,'r');
  if fid < 1
    warning ('Could not open job id list file %s\n', ctrl.job_id_file);
    ctrl.job_id_list = [];
    ctrl.error_mask = [];
    ctrl.job_status = [];
    return;
  end
  ctrl.job_id_list = textscan(fid,'%f');
  fclose(fid);
  ctrl.job_id_list = ctrl.job_id_list{1};
  ctrl.job_id = length(ctrl.job_id_list);
end
ctrl.job_id = ctrl.job_id + 1;
job_id = ctrl.job_id;

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

%% Save input arguments to the input arguments directory
% save(in_path,'taskfunction','argsin','num_args_out','param'); <-- % REPLACED THIS LINE TO FIX IU FILESERVER ISSUE
attempts = 0;
in_path = fullfile(ctrl.in_path_dir,sprintf('in_%d.mat',job_id));
while attempts < 3
  try
    save(in_path,'taskfunction','argsin','num_args_out','param','-v7.3');
    attempts = -1;
    break;
  catch ME
    ME
    attempts = attempts + 1;
  end
end
if attempts >= 0
  keyboard
end

%% Not actually submitting the job right now... putting it into the submission queue
ctrl.submission_queue = cat(2,ctrl.submission_queue,job_id);
new_job_status = 'T';
new_job_id = -1;
%% Get the torque job id for this job and store that
ctrl.job_id_list(end+1) = new_job_id;
ctrl.job_status(end+1) = new_job_status;
ctrl.error_mask(end+1) = 0;
% Write new job ID to job_id file
fid = fopen(ctrl.job_id_file,'a');
fprintf(fid,'%-20d\n', ctrl.job_id_list(end));
fclose(fid);

%% Submit this job to torque if submission queue is full enough
if strcmpi(ctrl.sched.submit_mode,'group') && length(ctrl.submission_queue) >= ctrl.sched.group_size
  ctrl = torque_submit_job(ctrl,ctrl.submission_queue);
  ctrl.submission_queue = [];
end

if ~isempty(param.notes)
  fprintf('  sub %d:%d/%d: %s (%s)\n', ctrl.batch_id, job_id, ctrl.job_id_list(end), param.notes, datestr(now));
end

end
