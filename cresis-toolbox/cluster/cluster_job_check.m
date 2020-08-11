function is_cluster_job = cluster_job_check()
% is_cluster_job = cluster_job_check()
%
% Call to determine if the current process is running on the cluster. Uses
% the param.cluster.is_cluster_job field.
%
% Inputs:
% global variable gRadar (uses gRadar.cluster.is_cluster_job)
%
% Outputs:
% is_cluster_job: logical scalar. True if this current process is a cluster
% job called from cluster_job.m (param.cluster.is_cluster_job is true).
% Returns false otherwise.
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

global gRadar;
if isfield(gRadar,'cluster') && isfield(gRadar.cluster,'is_cluster_job') && gRadar.cluster.is_cluster_job
  is_cluster_job = true;
else
  is_cluster_job = false;
end
