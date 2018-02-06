function cluster_cleanup(ctrl)
% cluster_cleanup(ctrl)
%
% Delete batches (jobs from queue and temporary files)
%
% Inputs:
% ctrl = Several options which specify which batches to act on
%   1. Pass in a cluster batch ctrl structure (only needs "batch_dir"
%      defined)
%   2. A vector of batch ids to apply hold to
%
% Examples:
%   cluster_cleanup(ctrl)
%   cluster_cleanup([1 3])
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

if nargin == 0
  answer = input('Are you sure you want to cleanup all cluster jobs? [y/N] ','s');
  if isempty(regexpi(answer,'y'))
    return
  end
end

%% Handle case where batch IDs have been passed in
if nargin == 0 || ~isstruct(ctrl)
  ctrls = cluster_get_batch_list;
  for batch_idx = 1:length(ctrls)
    if nargin == 0 || any(ctrls{batch_idx}.batch_id == ctrl)
      fprintf('  Deleting jobs in batch %d\n', ctrls{batch_idx}.batch_id);
      cluster_cleanup(ctrls{batch_idx});
    else
      fprintf('  Skipping %d\n', ctrls{batch_idx}.batch_id);
    end
  end
  return;
end

%% Delete the jobs on the cluster
if any(strcmpi(ctrl.cluster.type,{'torque','matlab','slurm'}))
  ctrl = cluster_get_batch(ctrl,ctrl.batch_id,false);
  
  % For each job in the batch, delete the job
  for job_id = 1:length(ctrl.job_id_list)
    if ctrl.job_status(job_id) ~= 'C'
      % Only delete jobs that have not been completed (completed jobs
      % are effectively deleted already)
      if strcmpi(ctrl.cluster.type,'torque')
        cmd = sprintf('qdel -W 60 -a %i', ctrl.job_id_list(job_id));
        try; [status,result] = system(cmd); end
      
      elseif strcmpi(ctrl.cluster.type,'matlab')
        for job_idx = length(ctrl.cluster.jm.Jobs):-1:1
          if ~isempty(ctrl.cluster.jm.Jobs(job_idx).ID == ctrl.job_id_list)
            try; delete(ctrl.cluster.jm.Jobs(job_idx)); end;
          end
        end
        
      elseif strcmpi(ctrl.cluster.type,'slurm')
        cmd = sprintf('scancel %i', ctrl.job_id_list(job_id));
        try; [status,result] = system(cmd); end
        
      end
    end
  end
end

%% Finally, remove the batch directory containing all the batch information
if exist(ctrl.batch_dir,'dir')
  fprintf('  %s: Removing %s\n', mfilename, ctrl.batch_dir);
  robust_rmdir(ctrl.batch_dir);
end

end
