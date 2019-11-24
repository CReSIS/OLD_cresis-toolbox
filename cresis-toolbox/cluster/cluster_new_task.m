function [ctrl,task_id] = cluster_new_task(ctrl,sparam,dparam,varargin)
% [ctrl,task_id] = cluster_new_task(ctrl,sparam,dparam,varargin)
%
% Creates a new job. Calls to this function need to be proceeded by
% a single call to cluster_new_batch.m.
%
% Inputs:
% ctrl: ctrl structure returned from cluster_new_batch
%  .sched: scheduler structure
%  .in_fn_dir: input arguments directory
%  .out_fn_dir: output arguments directory
%  .stdout_fn_dir: standard output directory
%  .error_fn_dir: error directory
% sparam: Struct containing static arguments to task
%  .task_function: function handle of job, this function handle tells
%    cluster_job.m what to run
%  .argsin: cell vector of input arguments (default is {})
%  .num_args_out: number of output arguments to expect (default is 0)
%  .notes: optional note to print after successful submission of job
%    default is '' (nothing is written out)
%  .cpu_time: maximum cpu time of this task in seconds
%  .mem: maximum memory usage of this task in bytes (default is 0)
%  .file_success: cell array of files that must exist to determine the
%    task a success. If the file is a .mat file, then it must have the
%    file_version variable and not be marked for deletion.
% dparam: Dynamic arguments to task (will be merged with sparam when the
%   task runs). This structure has all the same fields as sparam. After
%   merging, sparam and dparam must have all the required fields. It does
%   not matter where the field is set. For example, "task_function" could
%   be set from dparam and "mem" could be set from sparam. Note that if a
%   field is set in both structures, then dparam will overrule sparam.
%
% Outputs:
% ctrl: updated ctrl structure with new job
% task_id: ID of task (integer, starts counting from one and never repeats)
%
% EXAMPLES:
%   See cluster_submit_batch.m
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

%% Input checking
if ~exist('sparam','var')
  error('sparam must be specified');
end
if ~exist('dparam','var')
  error('dparam must be specified');
end
if ~isfield(sparam,'task_function') || isempty(sparam.task_function)
  error('sparam.task_function must be specified');
end
if isempty(which(sparam.task_function))
  error('sparam.task_function (%s) is not in the path', sparam.task_function);
end
if isfield(dparam,'task_function') && isempty(which(dparam.task_function))
  error('dparam.task_function is not in the path');
end
if ~isfield(sparam,'argsin') || isempty(sparam.argsin)
  sparam.argsin = {};
end
if ~isfield(sparam,'num_args_out')
  sparam.num_args_out = 0;
end
if ~isfield(sparam,'notes')
  sparam.notes = '';
end
if ~isfield(sparam,'cpu_time') || isempty(sparam.cpu_time)
  sparam.cpu_time = 0;
end
if ~isfield(sparam,'mem') || isempty(sparam.mem)
  sparam.mem = 0;
end
if ~isfield(sparam,'success') || isempty(sparam.success)
  sparam.success = '';
end
if ~isfield(sparam,'file_success') || isempty(sparam.file_success)
  sparam.file_success = {};
end
sparam.file_version = ctrl.cluster.file_version;
  
dparam_save = 1;
for param_idx = 1:2:length(varargin)
  eval(sprintf('%s = varargin{param_idx+1};', varargin{param_idx}));
end

%% Check for hold on this batch
if exist(fullfile(ctrl.batch_dir,'hold'), 'file')
  % Hold file exists
  keyboard
end

%% Get the new task id
if ~isfield(ctrl,'task_id')
  [fid,msg] = fopen(ctrl.task_id_file,'r');
  if fid < 1
    warning ('Could not open job id list file %s\n', ctrl.task_id_file);
    ctrl.job_id_list = [];
    ctrl.error_mask = [];
    ctrl.job_status = [];
    return;
  end
  ctrl.job_id_list = textscan(fid,'%f');
  fclose(fid);
  ctrl.job_id_list = ctrl.job_id_list{1};
  ctrl.task_id = length(ctrl.job_id_list);
end
ctrl.task_id = ctrl.task_id + 1;
task_id = ctrl.task_id;


%% Create task input file
if task_id == 1
  static_in_fn = fullfile(ctrl.in_fn_dir,'static.mat');
  static_param = sparam;
  robust_save(static_in_fn,ctrl.cluster.file_version,'static_param');
  if dparam_save
    dynamic_in_fn = fullfile(ctrl.in_fn_dir,'dynamic.mat');
    dynamic_param = [];
    dynamic_param.dparam{task_id} = dparam;
    robust_save(dynamic_in_fn,ctrl.cluster.file_version,'-struct','dynamic_param','dparam');
  else
    ctrl.dparam{task_id} = dparam;
  end
else
  % Add the dynamic parameters
  if dparam_save
    dynamic_in_fn = fullfile(ctrl.in_fn_dir,'dynamic.mat');
    dynamic_param = load(dynamic_in_fn);
    dynamic_param.dparam{task_id} = dparam;
    robust_save(dynamic_in_fn,ctrl.cluster.file_version,'-struct','dynamic_param','dparam');
  else
    ctrl.dparam{task_id} = dparam;
  end
end
param = merge_structs(sparam,dparam);

%% Add task to ctrl structure
ctrl.submission_queue = cat(2,ctrl.submission_queue,task_id);
new_job_status = 'T';
new_job_id = -1;
ctrl.job_id_list(end+1) = new_job_id;
ctrl.job_status(end+1) = new_job_status;
ctrl.error_mask(end+1) = 0;
ctrl.retries(end+1) = 0;
ctrl.notes{end+1} = param.notes;
ctrl.cpu_time(end+1) = param.cpu_time;
ctrl.mem(end+1) = param.mem;
ctrl.success{end+1} = param.success;
ctrl.file_success{end+1} = param.file_success;
ctrl.cpu_time_actual(end+1) = -1;

%% Write new task ID to task_id file
fid = fopen(ctrl.job_id_fn,'a');
fprintf(fid,'%-20d\n', ctrl.job_id_list(end));
fclose(fid);

%% Print out debug information to stdout
if ~isempty(param.notes)
  fprintf('  task %d:%d: %s (%s)\n', ctrl.batch_id, task_id, param.notes, datestr(now));
end

end
