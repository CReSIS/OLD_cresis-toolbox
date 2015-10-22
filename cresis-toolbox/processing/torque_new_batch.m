function ctrl = torque_new_batch(param)
% ctrl = torque_new_batch(param)
%
% Creates a new batch
% 1. Creates a temporary batch directory in param.data_location
% 2. Creates stdout, error, in, and out directories
% 3. Initializes the ctrl structure
%
% Input:
%  param = structure (like gRadar)
%   .sched structure (like gRadar.sched)
%    .data_location = location to store temporary job files
%
% Output:
% ctrl = batch control structure
%  .sched = equal to param.sched structure
%  .batch_dir = path to directory containing batch information
%  .job_id_file = path to file containing torque job ids
%  .in_path_dir = input arguments directory
%  .out_path_dir = output arguments directory
%  .stdout_path_dir = stdout directory
%  .error_path_dir = error directory
%  .user_file = path of file containing the username who owns batch
%  .user = string containing username who owns batch
%  .job_id = number of jobs created, N (which is zero since we just created
%    this batch)
%  .job_id_list = Nx1 vector of torque job IDs
%  .job_status = Nx1 vector of job status
%  .error_mask = Nx1 vector of error status
%
% Author: John Paden
%
% See also: torque_batch_list torque_cleanup torque_compile ...
%   torque_create_task torque_hold torque_job_list torque_job_status ...
%   torque_new_batch torque_print torque_rerun

if ~exist('param','var') || isempty(param)
  global gRadar;
  param = gRadar;
end

ctrl.sched = param.sched;

%% Create directory to store temporary files
% Find the first unique and unused batch_id
% Assign batch_id, batch_dir
ctrls = torque_batch_list(param);

ctrl.batch_id = 1;
done = 0;
while ~done
  done = 1;
  for batch_idx = 1:length(ctrls)
    if ctrls{batch_idx}.batch_id == ctrl.batch_id
      done = 0;
      ctrl.batch_id = ctrl.batch_id + 1;
      break;
    end
  end
end

% Each batch of jobs creates a unique directory
[tmp tmp_name] = fileparts(tempname);
ctrl.batch_dir = fullfile(param.sched.data_location,sprintf('batch_%i_%s', ctrl.batch_id, tmp_name));

%% Initialize ctrl structure and create batch directory and 4 subdirectories
ctrl.job_id = 0;
ctrl.job_id_list = [];
ctrl.job_status = [];
ctrl.error_mask = [];
ctrl.submission_queue = [];

ctrl.in_path_dir = fullfile(ctrl.batch_dir,'in');
mkdir(ctrl.in_path_dir)
ctrl.out_path_dir = fullfile(ctrl.batch_dir,'out');
mkdir(ctrl.out_path_dir)
ctrl.stdout_path_dir = fullfile(ctrl.batch_dir,'stdout');
mkdir(ctrl.stdout_path_dir)
ctrl.error_path_dir = fullfile(ctrl.batch_dir,'error');
mkdir(ctrl.error_path_dir)

% Create user file
cmd = 'id -un';
try
  [status,result] = system(cmd);
catch
  cmd
  warning('system call failed');
  %   keyboard;
end
ctrl.user = result(1:end-1);
ctrl.user_file = fullfile(ctrl.batch_dir,'user_file');
% Write username to user file
fid = fopen(ctrl.user_file,'w');
fprintf(fid,'%s\n', ctrl.user);
fclose(fid);

% Create job list file
ctrl.job_id_file = fullfile(ctrl.batch_dir,'job_id_file');
% Create empty job list file
fid = fopen(ctrl.job_id_file,'w');
fclose(fid);

return;
