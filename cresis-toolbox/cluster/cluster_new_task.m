function [ctrl,task_id] = cluster_new_task(ctrl,sparam,dparam)
% [ctrl,task_id] = cluster_new_task(ctrl,sparam,dparam)
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
% sparam:
%  .task_function: function handle of job, this function handle tells
%    cluster_job.m what to run
%  .argsin: cell vector of input arguments (default is {})
%  .num_args_out: number of output arguments to expect (default is 0)
%  .notes: optional note to print after successful submission of job
%    default is '' (nothing is written out)
%  .cpu_time: maximum cpu time of this task in seconds
%  .mem: maximum memory usage of this task in bytes (default is 0)
%  .success: string to be evaluated by "eval" to determine task success
%    (default is 'success=1;')
% dparam:
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
  error('sparam.cpu_time must be specified');
end
if ~isfield(sparam,'mem') || isempty(sparam.mem)
  sparam.mem = 0;
end
if ~isfield(sparam,'success') || isempty(sparam.success)
  sparam.success = '';
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
  robust_save(static_in_fn,'static_param');
  dynamic_in_fn = fullfile(ctrl.in_fn_dir,'dynamic.mat');
  dparam_task_field = sprintf('dparam_%d',task_id);
  dynamic_param = struct();
  dynamic_param.(dparam_task_field) = dparam;
  robust_save(dynamic_in_fn,'-v7.3','-struct','dynamic_param',dparam_task_field);
else
  % Add the dynamic parameters
  dynamic_in_fn = fullfile(ctrl.in_fn_dir,'dynamic.mat');
  dparam_task_field = sprintf('dparam_%d',task_id);
  dynamic_param = struct();
  dynamic_param.(dparam_task_field) = dparam;
  robust_save(dynamic_in_fn,'-append','-struct','dynamic_param',dparam_task_field);
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

%% Write new task ID to task_id file
fid = fopen(ctrl.job_id_fn,'a');
fprintf(fid,'%-20d\n', ctrl.job_id_list(end));
fclose(fid);

%% Print out debug information to stdout
if ~isempty(param.notes)
  fprintf('  task %d:%d: %s (%s)\n', ctrl.batch_id, task_id, param.notes, datestr(now));
end

end
