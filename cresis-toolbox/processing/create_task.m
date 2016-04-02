function [ctrl,job_id,task_id] = create_task(ctrl,fh,num_out_args,in_args)
%
% Inputs:
%  ctrl = control structure controlling the batch behavior
%   Fields assigned by user:
%    .sched = structure with scheduler information
%      .max_in_queue = 
%      .max_tasks_per_jobs = 
%    .fd = file dependencies for job creation
%   Fields assigned by this function:
%    .jobs = cell vector of structures
%      .job = job handle
%      .status = status of job wrt this function (queuing, submitted, done)
%      .tasks = list of task handles
%      .args_out = output arguments from tasks
%      .error_mask = a bit mask vector the same size as the list of
%         task handles indicating up to three events:
%         1: an error message has been printed for this task and you should
%            be able to scroll up and find it
%            (task failed due to uncaught exception from call to "error")
%         2: return variable success was empty (tasks probably never ran)
%         4: return variable success was zero (tasks failed, but returned normally)
%      .print_mask = a boolean mask vector the same size as the list of
%         task handles indicating if an error message has been printed
%         for this task
%  fh = function handle to be executed
%  num_out_args = number of output arguments from fh 
%  in_args = input arguments for fh
%  
%
% Output:
%  ctrl = control structure which is modified by this function
%  task = returned by Matlab's createTask command
%
%
% Used with the torque scheduler for two reasons:
%  1. When the number of jobs is large, the creation of tasks is really
%     slow, so this function allows tasks to be submitted in smaller
%     batches.
%  2. The torque scheduler has a maximum number of jobs in the queue
%     that is far less than the number of jobs being submitted.  So this
%     function does them in batches.
%
% Used with Matlab scheduler for two reasons:
%  1. The matlab scheduler stores all the job information in RAM so that
%     if there are too many tasks created at once, there may be a memory
%     allocation failure or memory management inefficiencies.
%  2. Nice to use the same code for IU and KU
%
% Other requirements:
%  1. All tasks must return a boolean success flag (1 = success)

if ~isfield(ctrl.sched,'ver')
  % Two versions supported. Original version of Matlab interface
  % (e.g. Matlab 2011b) is "1" and new version (e.g. Matlab 2014b) is "2"
  ctrl.sched.ver = 1;
end

if strcmpi(ctrl.cmd,'task')
  if length(ctrl.jobs) == 0 || (~isempty(ctrl.sched.max_tasks_per_jobs) && length(ctrl.jobs{end}.tasks) == ctrl.sched.max_tasks_per_jobs)
    % ====================================================================
    % If this is the first job or the last job is full, create a new job
    % ====================================================================
    if strcmp(ctrl.sched.type,'jobmanager')
      if isfield(ctrl.sched,'name') && ~isempty(ctrl.sched.name)
        jm = findResource('scheduler','type',ctrl.sched.type,'LookupUrl', ...
          ctrl.sched.url,'name',ctrl.sched.name);
      else
        jm = findResource('scheduler','type',ctrl.sched.type,'LookupUrl', ...
          ctrl.sched.url);
      end
      ctrl.jobs{end+1}.job = createJob(jm);
      set(ctrl.jobs{end}.job,'MaximumNumberOfWorkers',ctrl.sched.cluster_size);
      set(ctrl.jobs{end}.job,'RestartWorker',true);
    elseif strcmp(ctrl.sched.type,'torque')
      if isempty(ctrl.sched.url)
        jm = findResource('scheduler','type',ctrl.sched.type);
      else
        jm = findResource('scheduler','type',ctrl.sched.type,'LookupUrl',ctrl.sched.url);
      end
      set(jm,'DataLocation',ctrl.sched.data_location);
      set(jm,'SubmitArguments',ctrl.sched.submit_arguments);
      ctrl.jobs{end+1}.job = createJob(jm);
    elseif strcmpi(ctrl.sched.type,'local')
      if ctrl.sched.ver == 1
        jm = findResource('scheduler','type',ctrl.sched.type);
        if ctrl.sched.cluster_size < get(jm,'ClusterSize')
          set(jm,'ClusterSize',ctrl.sched.cluster_size);
        end
        if isfield(ctrl.sched,'data_location') && ~isempty(ctrl.sched.data_location)
          set(jm,'DataLocation',ctrl.sched.data_location);
        end
      else
        jm = parcluster();
        if ctrl.sched.cluster_size < get(jm,'NumWorkers')
          set(jm,'NumWorkers',ctrl.sched.cluster_size);
        end
        if isfield(ctrl.sched,'data_location') && ~isempty(ctrl.sched.data_location)
          set(jm,'JobStorageLocation',ctrl.sched.data_location);
        end
      end
      ctrl.jobs{end+1}.job = createJob(jm);
    end
      if ctrl.sched.ver == 1
        set(ctrl.jobs{end}.job,'FileDependencies',ctrl.fd);
      else
        set(ctrl.jobs{end}.job,'AdditionalPaths',ctrl.sched.AdditionalPaths);
      end
    ctrl.jobs{end}.status = 'queuing';
  end
  
  if ctrl.num_tasks_in_queue == ctrl.sched.max_in_queue
    % ====================================================================
    % Queue is full, need to wait until a job/task completes
    % ====================================================================
    
    % ------------------------------------------------------------------
    % Wait for a task to complete, by using a time out
    % and going through jobs in a round robin fashion.
    while ctrl.num_tasks_in_queue == ctrl.sched.max_in_queue
      ctrl = create_task_status(ctrl);
      if ctrl.num_tasks_in_queue == ctrl.sched.max_in_queue
        % Pause to avoid constantly polling
        if strcmp(ctrl.sched.type,'torque')
          % qstat operation is expensive on torque, avoid calling too often
          pause(30);
        else
          pause(2);
        end
      end
    end
    
    % ------------------------------------------------------------------
    % Destroy jobs for which all tasks have completed
    ctrl = create_task_cleanup(ctrl, false);
  end
  
  % ------------------------------------------------------------------
  % Submit task
  if ~isfield(ctrl.jobs{end},'tasks')
    ctrl.jobs{end}.tasks(1) = createTask(ctrl.jobs{end}.job,fh,num_out_args,in_args);
  else
    ctrl.jobs{end}.tasks(end+1) = createTask(ctrl.jobs{end}.job,fh,num_out_args,in_args);
  end
  job_id = length(ctrl.jobs);
  task_id = length(ctrl.jobs{end}.tasks);
  ctrl.num_tasks_in_queue = ctrl.num_tasks_in_queue + 1;
  
  if length(ctrl.jobs{end}.tasks) == ctrl.sched.max_tasks_per_jobs
    % --------------------------------------------------------
    % Submit job since full of tasks
    fprintf('Submitting %d tasks in job %s %d (%s)\n', ...
      length(ctrl.jobs{end}.tasks), get(ctrl.jobs{end}.job,'Name'), get(ctrl.jobs{end}.job,'ID'), datestr(now));
    submit(ctrl.jobs{end}.job);
    ctrl.jobs{end}.error_mask = zeros(size(ctrl.jobs{end}.tasks));
    ctrl.jobs{end}.print_mask = zeros(size(ctrl.jobs{end}.tasks));
    ctrl.jobs{end}.error_idxs = [];
    ctrl.jobs{end}.status = 'submitted';
  end
  
elseif strcmpi(ctrl.cmd,'init')
  % ====================================================================
  % Initialize submission variables
  % ====================================================================
  ctrl.jobs = {};
  ctrl.num_tasks_in_queue = 0;
  ctrl.error_mask = 0;
  
elseif strcmpi(ctrl.cmd,'done')
  % ====================================================================
  % Submit remaining tasks and wait for all jobs to complete and clean up
  % ====================================================================
  if length(ctrl.jobs) == 0
    return;
  end
  
  % --------------------------------------------------------------------
  % Submit remaining tasks
  if strcmpi(ctrl.jobs{end}.status, 'queuing')
    % --------------------------------------------------------
    % Submit job since full of tasks
    fprintf('Submitting %d tasks in job %s %d (%s)\n', ...
      length(ctrl.jobs{end}.tasks), get(ctrl.jobs{end}.job,'Name'), get(ctrl.jobs{end}.job,'ID'), datestr(now));
    submit(ctrl.jobs{end}.job);
    ctrl.jobs{end}.error_mask = zeros(size(ctrl.jobs{end}.tasks));
    ctrl.jobs{end}.print_mask = zeros(size(ctrl.jobs{end}.tasks));
    ctrl.jobs{end}.error_idxs = [];
    ctrl.jobs{end}.status = 'submitted';
  end

  % --------------------------------------------------------------------
  % Wait for all tasks to complete from each job, by using a timeout
  % and going through jobs in a round robin fashion.
  while ctrl.num_tasks_in_queue > 0
    ctrl = create_task_status(ctrl);
    if ctrl.num_tasks_in_queue > 0
      % Pause to avoid constantly polling
      if strcmp(ctrl.sched.type,'torque')
        % qstat operation is expensive on torque, avoid calling too often
        pause(30);
      else
        pause(2);
      end
    end
  end
  
  % --------------------------------------------------------------------
  % Destroy jobs for which all tasks have completed
  ctrl = create_task_cleanup(ctrl, true);
  
end

return;


% ====================================================================
% ====================================================================
% ctrl = create_task_status(ctrl)
%
% Function for getting the status of each task
% ====================================================================
% ====================================================================
function ctrl = create_task_status(ctrl)

for job_idx = 1:length(ctrl.jobs)
  if ~strcmpi(ctrl.jobs{job_idx}.status,'submitted')
    continue;
  end
  if ctrl.sched.ver == 1
    get(ctrl.jobs{job_idx}.job,'state');
  else
    get(ctrl.jobs{job_idx}.job,'State');
  end
  for task_idx = 1:length(ctrl.jobs{job_idx}.tasks)
    % Check for errors (ME = matlab exceptions)
    ME = get(ctrl.jobs{job_idx}.tasks(task_idx),'Error');
    
    %     if isempty(ME)
    %       keyboard;
    %     end
    if (isstruct(ME) || strcmpi(class(ME),'MException')) ...
        && ~isempty(ME) && (~isempty(ME.identifier) || ~isempty(ME.message)) ...
        && ~ctrl.jobs{job_idx}.error_mask(task_idx)
      ctrl.jobs{job_idx}.error_mask(task_idx) = 1;
      ctrl.jobs{job_idx}.error_idxs(end+1) = task_idx;
      warning('Job %d Task %d: %s\n  %s: line %d', job_idx, task_idx, ...
        ME.message, ME.stack(1).name, ME.stack(1).line);
    end;
    % Check for state
    if ctrl.sched.ver == 1
      job_done = (strcmpi(get(ctrl.jobs{job_idx}.tasks(task_idx),'state'),'finished') ...
          || strcmpi(get(ctrl.jobs{job_idx}.tasks(task_idx),'state'),'failed'));
    else
      job_done = (strcmpi(get(ctrl.jobs{job_idx}.tasks(task_idx),'State'),'finished') ...
          || strcmpi(get(ctrl.jobs{job_idx}.tasks(task_idx),'State'),'failed'));
    end
    if ~ctrl.jobs{job_idx}.print_mask(task_idx) && job_done
      ctrl.jobs{job_idx}.print_mask(task_idx) = 1;
      ctrl.num_tasks_in_queue = ctrl.num_tasks_in_queue - 1;
      fprintf('  %d: %d of %d tasks finished (%s)\n', ...
        job_idx, sum(ctrl.jobs{job_idx}.print_mask), ...
        numel(ctrl.jobs{job_idx}.tasks), datestr(now));
    end
  end
end

return;


% ====================================================================
% ====================================================================
% ctrl = create_task_cleanup(ctrl, wait_for_completion)
%
% Function for cleaning up tasks that have completed. If
% wait_for_completion is true, then all jobs will be cleaned up
% and the function will wait until they are done before doing so.
% ====================================================================
% ====================================================================
function ctrl = create_task_cleanup(ctrl, wait_for_completion)

for job_idx = 1:length(ctrl.jobs)
  if strcmpi(ctrl.jobs{job_idx}.status,'done')
    continue;
  end
  
  if ctrl.sched.ver == 1
    job_done = (strcmpi(get(ctrl.jobs{job_idx}.job,'state'),'finished') ...
      || strcmpi(get(ctrl.jobs{job_idx}.job,'state'),'failed'));
  else
    job_done = (strcmpi(get(ctrl.jobs{job_idx}.job,'State'),'finished') ...
      || strcmpi(get(ctrl.jobs{job_idx}.job,'State'),'failed'));
  end
  
  if wait_for_completion
    while ~(strcmpi(ctrl.jobs{job_idx}.status,'submitted') && job_done)
      fprintf('Waiting for Matlab to change job state to finished/failed for job %d\n', job_idx);
      % Pause to avoid constantly polling
      if strcmp(ctrl.sched.type,'torque')
        % qstat operation is expensive on torque, avoid calling too often
        pause(30);
      else
        pause(2);
      end
    end
  else
    if ~(strcmpi(ctrl.jobs{job_idx}.status,'submitted') && job_done)
      continue;
    end
  end
  
  % -----------------------------------------------------------------
  % Print out any unprinted messages
  for task_idx = 1:length(ctrl.jobs{job_idx}.tasks)
    % Check for errors (ME = matlab exceptions)
    ME = get(ctrl.jobs{job_idx}.tasks(task_idx),'Error');
    if (isstruct(ME) || strcmpi(class(ME),'MException')) ...
        && ~isempty(ME) && (~isempty(ME.identifier) || ~isempty(ME.message)) ...
        && ~ctrl.jobs{job_idx}.error_mask(task_idx)
      ctrl.jobs{job_idx}.error_mask(task_idx) = 1;
      ctrl.jobs{job_idx}.error_idxs(end+1) = task_idx;
      warning('Job %d Task %d: %s\n  %s: line %d', job_idx, task_idx, ...
        ME.message, ME.stack(1).name, ME.stack(1).line);
    end;
    % Check for state
    if ctrl.sched.ver == 1
      job_done = (strcmpi(get(ctrl.jobs{job_idx}.tasks(task_idx),'state'),'finished') ...
          || strcmpi(get(ctrl.jobs{job_idx}.tasks(task_idx),'state'),'failed'));
    else
      job_done = (strcmpi(get(ctrl.jobs{job_idx}.tasks(task_idx),'State'),'finished') ...
          || strcmpi(get(ctrl.jobs{job_idx}.tasks(task_idx),'State'),'failed'));
    end
    if ~ctrl.jobs{job_idx}.print_mask(task_idx) && job_done
      ctrl.jobs{job_idx}.print_mask(task_idx) = 1;
      ctrl.num_tasks_in_queue = ctrl.num_tasks_in_queue - 1;
      fprintf('  %d: %d of %d tasks finished (%s)\n', ...
        job_idx, sum(ctrl.jobs{job_idx}.print_mask), ...
        numel(ctrl.jobs{job_idx}.tasks), datestr(now));
    end
  end
  
  % -----------------------------------------------------------------
  % Get output arguments
  try
    ctrl.jobs{job_idx}.args_out = getAllOutputArguments(ctrl.jobs{job_idx}.job);
  catch ME
    warning(ME.getReport);
    keyboard
  end
  
  % -----------------------------------------------------------------
  % Check success flag and add to error queue any tasks that did
  % not set the success flag unless they are already in the queue.
  for task_idx = 1:length(ctrl.jobs{job_idx}.tasks)
    if size(ctrl.jobs{job_idx}.args_out,2) < 1 ...
        || isempty(ctrl.jobs{job_idx}.args_out{task_idx,1})
      ctrl.jobs{job_idx}.error_mask(task_idx) = ctrl.jobs{job_idx}.error_mask(task_idx) + 2;
      if all(ctrl.jobs{job_idx}.error_idxs ~= task_idx)
        ctrl.jobs{job_idx}.error_idxs(end+1) = task_idx;
      end
    elseif ctrl.jobs{job_idx}.args_out{task_idx,1} ~= 1
      ctrl.jobs{job_idx}.error_mask(task_idx) = ctrl.jobs{job_idx}.error_mask(task_idx) + 4;
      ctrl.jobs{job_idx}.error_idxs(end+1) = task_idx;
      if all(ctrl.jobs{job_idx}.error_idxs ~= task_idx)
        ctrl.jobs{job_idx}.error_idxs(end+1) = task_idx;
      end
    end
  end
  
  % -----------------------------------------------------------------
  % Check for existence of each type of error and set error mask
  % to this type of error if it has not already been set
  if any(mod(floor(ctrl.jobs{job_idx}.error_mask/1),2)) ...
      && ~mod(floor(ctrl.error_mask/1),2)
    ctrl.error_mask = ctrl.error_mask + 1;
  end
  if any(mod(floor(ctrl.jobs{job_idx}.error_mask/2),2)) ...
      && ~mod(floor(ctrl.error_mask/2),2)
    ctrl.error_mask = ctrl.error_mask + 2;
  end
  if any(mod(floor(ctrl.jobs{job_idx}.error_mask/4),2)) ...
      && ~mod(floor(ctrl.error_mask/4),2)
    ctrl.error_mask = ctrl.error_mask + 4;
  end
  
  % -----------------------------------------------------------------
  % Check for errors
  if ~isempty(ctrl.jobs{job_idx}.error_idxs)
    warning('Errors found in job idx %d. Task idxs are printed below.', job_idx);
    sort(ctrl.jobs{job_idx}.error_idxs)
    ctrl.done = false;
    if ctrl.sched.stop_on_fail && ctrl.error_mask ~= 2
      fprintf('After debugging the error, type dbcont\n');
      keyboard;
    else
      fprintf('Skipping debugging since this might be a recoverable error.\n');
      fprintf('If not, uncomment the keyboard command here.\n');
      % keyboard
    end
  end
  
  % -----------------------------------------------------------------
  % Clean up job
  destroy(ctrl.jobs{job_idx}.job);
  
  ctrl.jobs{job_idx}.status = 'done';
end

