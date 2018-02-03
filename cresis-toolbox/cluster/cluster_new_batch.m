function ctrl = cluster_new_batch(param)
% ctrl = cluster_new_batch(param)
%
% Creates a new batch
% 1. Creates a temporary batch directory in param.data_location
% 2. Creates stdout, error, in, and out directories
% 3. Initializes the ctrl structure
%
% Input:
%  param: structure with cluster information
%   .cluster: structure (default is gRadar.cluster)
%    .data_location: location to store temporary job files
%
% Output:
% ctrl = batch control structure
%  .cluster = set equal to param.cluster structure
%  .batch_dir = path to directory containing batch information
%  .job_id_fn = path to file containing cluster job ids associated with
%    each task
%  .in_fn_dir = input arguments directory
%  .out_fn_dir = output arguments directory
%  .stdout_fn_dir = stdout directory
%  .error_fn_dir = error directory
%  .user_file = path of file containing the username who owns batch
%  .user = string containing username who owns batch
%  .task_id = number of tasks created, N (which is zero since we just created
%    this batch)
%  .job_id_list = Nx1 vector of cluster job IDs associated with each task
%    (this should match the contents of the job_id_fn)
%  .job_status = Nx1 vector of job status
%  .error_mask = Nx1 vector of error status
%
% Example:
%  cluster_new_batch; % Usually called with no arguments
%
% Author: John Paden
%
% See also: cluster_batch_list cluster_cleanup cluster_compile ...
%   cluster_create_task cluster_hold cluster_job_list cluster_job_status ...
%   cluster_new_batch cluster_print cluster_rerun

%% Input arguments check
if ~exist('param','var') || isempty(param)
  global gRadar;
  param = gRadar;
end

ctrl.cluster = param.cluster;

if strcmpi(ctrl.cluster.type,'none')
  error('%s should not be called with ctrl.cluster.type="%s"', mfilename,ctrl.cluster.type);
end

if ~isfield(ctrl.cluster,'max_jobs_active') || isempty(ctrl.cluster.max_jobs_active)
  ctrl.cluster.max_jobs_active = 1;
end

if ~isfield(ctrl.cluster,'max_time_per_job') || isempty(ctrl.cluster.max_time_per_job)
  ctrl.cluster.max_time_per_job = 0;
end

if ~isfield(ctrl.cluster,'max_retries') || isempty(ctrl.cluster.max_retries)
  ctrl.cluster.max_retries = 1;
end

%% Create directory to store temporary files
% Find the first unique and unused batch_id
% Assign batch_id, batch_dir
ctrls = cluster_get_batch_list(param);

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
ctrl.batch_dir = fullfile(param.cluster.data_location,sprintf('batch_%i_%s', ctrl.batch_id, tmp_name));

%% Initialize ctrl structure and create batch directory and 4 subdirectories
ctrl.task_id = 0;
ctrl.job_id_list = [];
ctrl.job_status = '';
ctrl.error_mask = [];
ctrl.submission_queue = [];
ctrl.active_jobs = 0;

ctrl.in_fn_dir = fullfile(ctrl.batch_dir,'in');
mkdir(ctrl.in_fn_dir)
ctrl.out_fn_dir = fullfile(ctrl.batch_dir,'out');
mkdir(ctrl.out_fn_dir)
ctrl.stdout_fn_dir = fullfile(ctrl.batch_dir,'stdout');
mkdir(ctrl.stdout_fn_dir)
ctrl.error_fn_dir = fullfile(ctrl.batch_dir,'error');
mkdir(ctrl.error_fn_dir)

%% Get the job manager for the matlab cluster interface
if strcmpi(ctrl.cluster.type,'matlab')
  ctrl.cluster.jm = parcluster();
end

%% Create empty job id list file (one 20 character line per task)
ctrl.job_id_fn = fullfile(ctrl.batch_dir,'job_id_file');
fid = fopen(ctrl.job_id_fn,'w');
fclose(fid);

return;
