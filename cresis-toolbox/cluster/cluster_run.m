function ctrl_chain = cluster_run(ctrl_chain,block)
% ctrl_chain = cluster_run(ctrl_chain,block)
%
% Submits jobs in a list of batch chains. Each chain in the list runs in
% parallel. Batches within a chain are run in series.
%
% Inputs:
% ctrl_chain: cell array of chains that can be run in parallel
%  ctrl_chain{chain}: cell array of batches that must be run in series (stages)
%   ctrl_chain{chain}{stage}: control structure for a batch
% block: logical (if true, this function will block until all jobs are
%   completed). If false, this function just passes through all chains one
%   time (so that cluster status can be polled).
%
% Outputs:
% ctrl_chain: updated list of batch chains that was passed in
%
% Example:
% % If there is just one control structure to run called ctrl
% ctrl_chain = {{ctrl}};
% ctrl_chain = cluster_run(ctrl_chain);
%
% % Let ctrl1 be in chain 1, let ctrl2a and ctrl2b be in chain 2
% ctrl_chain = {{ctrl1},{ctrl2a,ctrl2b}};
% ctrl_chain = cluster_run(ctrl_chain);
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

if iscell(ctrl_chain)
  %% Input checking
  if ~exist('block','var') || isempty(block)
    block = true;
  end
  
  %% Traverse chain list
  active_stage = ones(numel(ctrl_chain),1);
  first_run = ones(numel(ctrl_chain),1);
  while any(isfinite(active_stage))
    for chain = 1:numel(ctrl_chain)
      if isfinite(active_stage(chain))
        % 1. There is at least one batch left to run in this chain
        ctrl = ctrl_chain{chain}{active_stage(chain)};
        
        % 2. If this is the first loop of cluster_run, force a complete
        %   update of the job status information.
        if first_run(chain)
          ctrl = cluster_get_batch(ctrl);
          first_run(chain) = false;
        else
          ctrl = cluster_update_batch(ctrl);
          pause(ctrl.cluster.stat_pause);
        end
        
        % 3. Submit jobs from the active stage for each parallel control structure
        %   ctrl.max_active_jobs.
        ctrl = cluster_run(ctrl);
        
        % 4. Update ctrl_chain
        ctrl_chain{chain}{active_stage(chain)} = ctrl;
        
        % 5. If all jobs completed in a batch and:
        %    If no errors, move to the next stage
        %    If errors and out of retries, stop chain
        if all(ctrl.job_status=='C')
          if ~any(ctrl.error_mask)
            active_stage(chain) = active_stage(chain) + 1;
            first_run(chain) = true;
            if active_stage(chain) > numel(ctrl_chain{chain})
              % Chain is complete
              active_stage(chain) = inf;
            end
          elseif all(ctrl.retries >= ctrl.cluster.max_retries | ~ctrl.error_mask)
            % Stop chain
            active_stage(chain) = -inf;
          end
        end
        
        % 6. Check to see if a hold has been placed on this batch
        if exist(ctrl.hold_fn,'file')
          fprintf('This batch has a hold. Run cluster_hold(ctrl) to remove. Either way, run dbcont to continue.\n');
          keyboard
        end
      end
    end
    if ~block
      break;
    end
  end

  failed_chains = find(active_stage == -inf);
  for chain=1:length(failed_chains)
    fprintf('Chain %d failed (%s)\n', chain, datestr(now));
    for stage=1:numel(ctrl_chain{chain})
      ctrl = ctrl_chain{chain}{stage};
      if ~any(ctrl.error_mask)
        fprintf('Stage %d succeeded\n', stage);
      else
        fprintf('  Stage %d failed (%d of %d tasks failed)\n', stage, sum(ctrl.error_mask~=0), length(ctrl.error_mask));
      end
    end
  end
  
elseif isstruct(ctrl_chain)
  ctrl = ctrl_chain;
  
  if strcmpi(ctrl.cluster.type,'none')
    return;
  end

  job_tasks = [];
  job_cpu_time = 0;
  job_mem = 0;
  while ~isempty(ctrl.submission_queue) && ctrl.active_jobs < ctrl.cluster.max_jobs_active
    % Get task from queue
    task_id = ctrl.submission_queue(1);

    if isempty(job_tasks) ...
        && ctrl.cluster.max_time_per_job > 0 ...
        && ctrl.cluster.max_time_per_job < job_cpu_time + ctrl.cpu_time(task_id)
      error('ctrl.cluster.max_time_per_job is less than task %d time %.0f', task_id, ctrl.cpu_time(task_id));
    end
    if ctrl.cluster.max_time_per_job < job_cpu_time + ctrl.cpu_time(task_id)
      [ctrl,new_job_id] = cluster_submit_job(ctrl,job_tasks,job_cpu_time,job_mem);
      fprintf('Submitted these tasks in cluster job %d/%d:\n  %d', ctrl.batch_id, new_job_id, job_tasks(1))
      if length(job_tasks) > 1
        fprintf(', %d', job_tasks(2:end));
      end
      fprintf('\n');
      job_tasks = [];
      job_cpu_time = 0;
      job_mem = 0;
      pause(ctrl.cluster.submit_pause);
    end
    job_tasks(end+1) = task_id;
    job_cpu_time = job_cpu_time + ctrl.cpu_time(task_id);
    job_mem = max(job_mem, ctrl.mem(task_id));
    ctrl.submission_queue = ctrl.submission_queue(2:end);
  end
  
  if ctrl.active_jobs < ctrl.cluster.max_jobs_active && ~isempty(job_tasks)
    [ctrl,new_job_id] = cluster_submit_job(ctrl,job_tasks,job_cpu_time,job_mem);
    fprintf('Submitted these tasks in cluster job %d/%d:\n  %d', ctrl.batch_id, new_job_id, job_tasks(1))
    if length(job_tasks) > 1
      fprintf(', %d', job_tasks(2:end));
    end
    fprintf('\n');
    pause(ctrl.cluster.submit_pause);
    
  else
    % Put jobs back in the queue because they can't be run yet
    ctrl.submission_queue = cat(2,job_tasks,ctrl.submission_queue);
  end
  
  % Return the updated ctrl
  ctrl_chain = ctrl;
end

return
