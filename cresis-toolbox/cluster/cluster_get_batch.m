function ctrl = cluster_get_batch(ctrl,batch_id,force_check)
% ctrl = cluster_get_batch(ctrl,batch_id,force_check)
%
% Updates task status information from the cluster. May be called from any
% matlab process since it rebuilds the ctrl structure.
%
% Inputs:
% ctrl = Default is to leave empty, []. However, if ctrl is provided:
%   ctrl.cluster.data_location is used
%   ctrl.batch_id is used
% batch_id: integer for which batch to get jobs info for
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
% See also: cluster_batch_list cluster_cleanup cluster_compile ...
%   cluster_create_task cluster_hold cluster_job_list cluster_job_status ...
%   cluster_new_batch cluster_print cluster_rerun

if ~exist('force_check','var') || isempty(force_check)
  force_check = true;
end

if isempty(ctrl)
  global gRadar;
  ctrl.cluster = gRadar.cluster;
end

if isfield(ctrl,'batch_id')
  batch_id = ctrl.batch_id;
end

ctrls = cluster_get_batch_list(ctrl);

found = 0;
for batch_idx = 1:length(ctrls)
  if ctrls{batch_idx}.batch_id == batch_id
    found = 1;
    ctrl = ctrls{batch_idx};
    break;
  end
end

if found == 0
  fprintf('Batch %d not found\n', batch_id);
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
ctrl.job_id = length(ctrl.job_id_list);

% Start by assuming that all jobs are in the to be submitted queue
ctrl.submission_queue = 1:numel(ctrl.job_id_list);
ctrl.job_status = char('T'*ones(1,numel(ctrl.job_id_list)));

% Start by assuming that there are no errors in any tasks
ctrl.error_mask = zeros(size(ctrl.job_status));

% Start by assuming no retries yet
ctrl.retries = zeros(size(ctrl.job_status));

job_status_found = zeros(size(ctrl.job_status));

if ~force_check
  return
end

%% Update task status for each task using cluster interface
ctrl.active_jobs = 0;
if any(strcmpi(ctrl.cluster.type,{'torque','matlab','slurm'}))
  if strcmpi(ctrl.cluster.type,'torque')
    % Runs qstat command
    % -----------------------------------------------------------------------
    [system,user_name] = robust_system('whoami');
    user_name = user_name(1:end-1);
    cmd = sprintf('qstat -u %s </dev/null', user_name);
    [status,result] = robust_system(cmd);
    
  elseif strcmpi(ctrl.cluster.type,'matlab')
    status = 0;
    result = 'NA';
    
  elseif strcmpi(ctrl.cluster.type,'slurm')
    [system,user_name] = robust_system('whoami');
    user_name = user_name(1:end-1);
    status = 0;
    result = 'NA';
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
      for job_idx = 1:length(ctrl.cluster.jm.Jobs)
        qstat_res{7}(job_idx,1) = ctrl.cluster.jm.Jobs(job_idx).ID;
        if strcmpi(ctrl.cluster.jm.Jobs(job_idx).State,'finished')
          qstat_res{5}{job_idx,1} = 'C';
        else
          ctrl.active_jobs = ctrl.active_jobs + 1;
          qstat_res{5}{job_idx,1} = 'Q';
        end
      end
      
    elseif strcmpi(ctrl.cluster.type,'slurm')
      qstat_res{5} = {};
      qstat_res{7} = [];
      
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
            fprintf(' QJob %d:%d/%d status changed to %s (%s)\n', ctrl.batch_id, task_id, ctrl.job_id_list(task_id), new_job_status, datestr(now))
            ctrl.job_status(task_id) = new_job_status;
          end
        end
      end
    end
  end
end

ctrl.active_jobs = 0;
for task_id = 1:length(ctrl.job_status)
  ctrl = cluster_update_task(ctrl,task_id);
end
