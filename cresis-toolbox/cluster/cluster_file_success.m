function error_mask = cluster_file_success(fns_list)
% error_mask = cluster_file_success(fns_list)
%
% Check a list of files to see if they are successfully created. Support
% function for using the cluster.
%
% Inputs:
% ctrl: ctrl structure returned from cluster_new_batch
%  .job_id_list = Nx1 vector of cluster job IDs
%  .job_status = Nx1 vector of job status
%  .out_fn_dir = string containing the output directory
%  .error_mask = Nx1 vector of error status
%  .out = up-to-Nx1 vector of cells containing the outputs for each
%    job as they complete (read in from the output.mat files)
% job_id: the job id to check up on
%
% Outputs:
% ctrl = updated ctrl structure
%
% C -  Job is completed after having run.
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

file_success_error = 1024;
file_success_corrupt_error = 2048;

error_mask = 0;
for file_idx=1:length(fns_list)
  fn = fns_list{file_idx};
  [~,~,fn_ext] = fileparts(fn);
  if ~exist(fn,'file')
    error_mask = bitor(error_mask,file_success_error);
  else
    if strcmpi(fn_ext,'.mat') && ~ct_file_lock_check(fn,4)
      error_mask = bitor(error_mask,file_success_corrupt_error);
    end
  end
end
