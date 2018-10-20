function ctrl = cluster_get_batch(ctrl,force_check,update_mode)
% ctrl = cluster_get_batch(ctrl,force_check,update_mode)
%
% Updates task status information from the cluster. May be called from any
% matlab process since it rebuilds the ctrl structure.
%
% Inputs:
% ctrl = Must be specified. Identifies the batch.
% force_check: default to true, forces cluster_update_task to be run on all
%   tasks
% update_mode: Scalar integer indicating mode of operation:
%   0: Only getting task information, but do not updating tasks.
%   1: Get task information AND update/print tasks as required (default)
%
% Outputs:
% ctrl = ctrl structure for the specified batch_id
%
% C -  Job is completed after having run/
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

if ~exist('force_check','var') || isempty(force_check)
  force_check = true;
end
if ~exist('update_mode','var') || isempty(update_mode)
  update_mode = 1;
end

if isnumeric(ctrl)
  ctrl_is_struct = false;
  batch_id = ctrl;
  global gRadar;
  ctrl = [];
  ctrl.cluster = gRadar.cluster;
else
  ctrl_is_struct = true;
end

if isfield(ctrl,'batch_id')
  batch_id = ctrl.batch_id;
end

if isempty(batch_id)
  error('%s called with no batch_id. Either ctrl.batch_id or batch_id must be specified.',mfilename);
end

ctrls = cluster_get_batch_list(ctrl);

match_idxs = [];
for batch_idx = 1:length(ctrls)
  if ctrls{batch_idx}.batch_id == batch_id
    match_idxs(end+1) = batch_idx;
  end
end

if isempty(match_idxs)
  error('Batch %d not found.', batch_id);
end
if ~ctrl_is_struct || ~isfield(ctrl,'batch_dir')
  if length(match_idxs) > 1
    
    fprintf('<strong>There are %d matches for this batch_id %d.</strong>\n', length(match_idxs), batch_id);
    fprintf('Select a batch from the list below:\n');
    % This for-loop from cluster_get_batch_list.m:
    for match_idxs_idx = 1:length(match_idxs)
      batch_idx = match_idxs(match_idxs_idx);
      % Create input filenames
      static_in_fn = fullfile(ctrls{batch_idx}.in_fn_dir,'static.mat');
      dynamic_in_fn = fullfile(ctrls{batch_idx}.in_fn_dir,'dynamic.mat');
      % Load input filenames
      if exist(static_in_fn,'file')
        sparam = load(static_in_fn);
      else
        warning('Missing %s', static_in_fn);
        sparam = [];
        sparam.static_param = [];
      end
      if exist(dynamic_in_fn,'file')
        tmp = load(dynamic_in_fn);
        ctrls{batch_idx}.dparam = tmp.dparam;
        if isempty(ctrls{batch_idx}.dparam)
          ctrls{batch_idx}.dparam = {[]};
        end
      else
        ctrls{batch_idx}.dparam = {[]};
      end
      
      fprintf('<strong>%d</strong>: Batch %d %s\n', match_idxs_idx, ctrls{batch_idx}.batch_id, ctrls{batch_idx}.batch_dir);
      
      param = merge_structs(sparam.static_param,ctrls{batch_idx}.dparam{1});
      if ~isfield(param,'task_function')
        param.task_function = '';
      end
      if ~isfield(param,'notes')
        param.notes = '';
      end
      fprintf('    task_function: %s\n', param.task_function);
      fprintf('    notes: %s\n', param.notes);
    end
    
    % Check with user
    uinput = [];
    while length(uinput)~=1 || ~isnumeric(uinput) || uinput < 1 || uinput > length(match_idxs)
      uinput = input('? ');
    end
    match_idxs = match_idxs(uinput);
  end
  
  ctrl = ctrls{match_idxs};
  
else
  found = false;
  for match_idxs_idx = 1:length(match_idxs)
    if strcmpi(ctrls{match_idxs(match_idxs_idx)}.batch_dir,ctrl.batch_dir)
      found = true;
      ctrl = merge_structs(ctrls{match_idxs(match_idxs_idx)}, ctrl);
      break;
    end
  end
  if ~found
    error('Batch %d found, but original batch seems to be deleted since ctrl.batch_dir does not match an existing batch.', batch_id);
  end
end

[fid,msg] = fopen(ctrl.job_id_fn,'r');
if fid < 1
  warning ('Could not open job id list file %s\n', ctrl.job_id_fn);
  ctrl.job_id_list = [];
  ctrl.error_mask = [];
  ctrl.job_status = [];
  return;
end
ctrl.job_id_list = textscan(fid,'%f');
fclose(fid);
ctrl.job_id_list = ctrl.job_id_list{1};
ctrl.task_id = length(ctrl.job_id_list);

% Start by assuming that all jobs are in the to be submitted queue
ctrl.submission_queue = 1:numel(ctrl.job_id_list);
ctrl.job_status = char('T'*ones(1,numel(ctrl.job_id_list)));

% Start by assuming that there are no errors in any tasks
ctrl.error_mask = zeros(size(ctrl.job_status));

% Start by assuming no retries yet
ctrl.retries = zeros(size(ctrl.job_status));

% Start by assuming zero cpu time and memory requirements
ctrl.cpu_time = zeros(size(ctrl.job_status));
ctrl.mem = zeros(size(ctrl.job_status));

job_status_found = zeros(size(ctrl.job_status));

%% Update task status for each task using cluster interface
ctrl.active_jobs = 0;
if any(strcmpi(ctrl.cluster.type,{'torque','matlab','slurm'}))
  if strcmpi(ctrl.cluster.type,'torque')
    % Runs qstat command
    % -----------------------------------------------------------------------
    [system,user_name] = robust_system('whoami </dev/null');
    user_name = user_name(1:end-1);
    cmd = sprintf('qstat -u %s </dev/null', user_name);
    [status,result] = robust_system(cmd);
    
  elseif strcmpi(ctrl.cluster.type,'matlab')
    status = 0;
    result = 'NA';
    
  elseif strcmpi(ctrl.cluster.type,'slurm')
    % Runs squeue command
    % -----------------------------------------------------------------------
    [system,user_name] = robust_system('whoami </dev/null');
    user_name = user_name(1:end-1);
    cmd = sprintf('squeue --users=%s </dev/null', user_name);
    [status,result] = robust_system(cmd);
  end
  
  % Parse qstat command results
  % -----------------------------------------------------------------------
  if ~isempty(result)
    if strcmpi(ctrl.cluster.type,'torque')
      qstat_res = textscan(result,'%s %s %s %s %s %s %s %s %s %s %s %s','Headerlines',5,'Delimiter',sprintf(' \t'),'MultipleDelimsAsOne',1);
      for idx = 1:size(qstat_res{1},1)
        qstat_res{1}{idx} = str2double(strtok(qstat_res{1}{idx},'.'));
        if qstat_res{10}{idx} ~= 'C'
          ctrl.active_jobs = ctrl.active_jobs + 1;
        end
      end
      qstat_res{7} = cell2mat(qstat_res{1});
      qstat_res{5} = qstat_res{10};
      
    elseif strcmpi(ctrl.cluster.type,'matlab')
      qstat_res{5} = {};
      qstat_res{7} = [];
      if ~isfield(ctrl.cluster,'jm')
        ctrl.cluster.jm = parcluster;
      end
      IDs = cell2mat({ctrl.cluster.jm.Jobs.ID});
      States = {ctrl.cluster.jm.Jobs.State};
      for job_idx = 1:length(IDs)
        qstat_res{7}(job_idx,1) = IDs(job_idx);
        if any(strcmpi(States(job_idx),{'finished','failed'}))
          qstat_res{5}{job_idx,1} = 'C';
        elseif any(strcmpi(States(job_idx),{'running'}))
          ctrl.active_jobs = ctrl.active_jobs + 1;
          qstat_res{5}{job_idx,1} = 'R';
        else
          ctrl.active_jobs = ctrl.active_jobs + 1;
          qstat_res{5}{job_idx,1} = 'Q';
        end
      end
      
    elseif strcmpi(ctrl.cluster.type,'slurm')
      qstat_res = textscan(result,'%s %s %s %s %s %s %s %s','Headerlines',1,'Delimiter',sprintf(' \t'),'MultipleDelimsAsOne',1);
      for idx = 1:size(qstat_res{1},1)
        qstat_res{1}{idx} = str2double(qstat_res{1}{idx});
        qstat_res{5}{idx} = qstat_res{5}{idx}(1);
        if qstat_res{5}{idx} ~= 'C'
          ctrl.active_jobs = ctrl.active_jobs + 1;
        end
      end
      qstat_res{7} = cell2mat(qstat_res{1});
      
    end
    % Loop through all the jobs that qstat returned
    % qstat_res{7}(idx): JOB ID
    % qstat_res{5}{idx}: JOB STATUS
    for idx = 1:size(qstat_res{5},1)
      task_ids = find(qstat_res{7}(idx)==ctrl.job_id_list);
      % Qstat returns all jobs, just look at jobs in this batch
      while ~isempty(task_ids)
        task_id = task_ids(1);
        task_ids = task_ids(2:end);
        job_status_found(task_id) = 1;
        if qstat_res{5}{idx} ~= ctrl.job_status(task_id)
          if ctrl.job_status(task_id) ~= 'C'
            % Only update job if it is not complete
            new_job_status = qstat_res{5}{idx};
            % Debug print
            if update_mode
              fprintf(' QJob %d:%d/%d status changed to %s (%s)\n', ctrl.batch_id, task_id, ctrl.job_id_list(task_id), new_job_status, datestr(now))
            end
            ctrl.job_status(task_id) = new_job_status;
          end
        end
      end
    end
  end
end

if ~force_check
  static_in_fn = fullfile(ctrl.in_fn_dir,'static.mat');
  dynamic_in_fn = fullfile(ctrl.in_fn_dir,'dynamic.mat');
  if isempty(ctrl.job_status)
    sparam = [];
  else
    sparam = load(static_in_fn);
  end
  for task_id = 1:length(ctrl.job_status)
    if ~isfield(ctrl,'dparam') || numel(ctrl.dparam) < task_id || isempty(ctrl.dparam{task_id})
      tmp = load(dynamic_in_fn);
      ctrl.dparam = tmp.dparam;
    end
    param = merge_structs(sparam.static_param,ctrl.dparam{task_id});
    ctrl.cpu_time(task_id) = param.cpu_time;
    ctrl.mem(task_id) = param.mem;
  end
  
  return
end

for task_id = 1:length(ctrl.job_status)
  ctrl = cluster_update_task(ctrl,task_id,update_mode);
end
