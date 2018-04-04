function [ctrl_chain,stats] = cluster_print_chain(ctrl_chain)
% [ctrl_chain,stats] = cluster_print_chain(ctrl_chain)
%
% Prints information about a list of batch chains. Also returns the
% information.
%
% Inputs:
% ctrl_chain: cell array of chains that can be run in parallel
%  ctrl_chain{chain}: cell array of batches that must be run in series (stages)
%   ctrl_chain{chain}{stage}: control structure for a batch
%
% Outputs:
% ctrl_chain: updated list of batch chains that was passed in
% stats: statistics about all the batch chains
%
% Example:
%   [ctrl_chain,stats] = cluster_print_chain(ctrl_chain);
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

if iscell(ctrl_chain)
  %% Traverse chain list
  stats.cpu_time = [];
  stats.mem = [];
  stats.error_mask = [];
  stats.job_status = [];
  stats.retries = [];
  for chain = 1:numel(ctrl_chain)
    fprintf('Chain %d\n', chain);
    for stage=1:numel(ctrl_chain{chain})
      fprintf('  Stage %d (Batch %d)\n', stage, ctrl_chain{chain}{stage}.batch_id);
      [ctrl_chain{chain}{stage},ctrl_stats] = cluster_print_chain(ctrl_chain{chain}{stage});
      stats.cpu_time = cat(2,stats.cpu_time,ctrl_stats.cpu_time);
      stats.mem = cat(2,stats.mem,ctrl_stats.mem);
      stats.error_mask = cat(2,stats.error_mask,ctrl_stats.error_mask);
      stats.retries = cat(2,stats.retries,ctrl_stats.retries);
      stats.job_status = cat(2,stats.job_status,ctrl_stats.job_status);
    end
  end
  fprintf('====================================================\n');
  fprintf('Number of tasks: %.0f, %.0f/%.0f/%.0f/%.0f C/R/Q/T, %.0f error, %.0f retries\n', ...
    numel(stats.cpu_time), sum(stats.job_status=='C'), sum(stats.job_status=='R'), sum(stats.job_status=='Q'), sum(stats.job_status=='T'), sum(stats.error_mask~=0), sum(stats.retries));
  fprintf('Max CPU time: %.0f min\n', max(stats.cpu_time)/60);
  fprintf('Max mem: %.0f MB\n', max(stats.mem)/1e6);
  fprintf('Mean CPU time: %.0f min\n', mean(stats.cpu_time)/60);
  fprintf('Mean mem: %.0f MB\n', mean(stats.mem)/1e6);  
  
elseif isstruct(ctrl_chain)
  ctrl = ctrl_chain;
  
  ctrl = cluster_get_batch(ctrl,[],[],0);
  
  stats.cpu_time = ctrl.cpu_time;
  stats.mem = ctrl.mem;
  stats.error_mask = ctrl.error_mask;
  stats.retries = ctrl.retries;
  stats.job_status = ctrl.job_status;
  
  % Create in/out filenames
  static_in_fn = fullfile(ctrl.in_fn_dir,'static.mat');
  
  % Read input
  sparam = load(static_in_fn);
  if ~isempty(ctrl.dparam)
    param = merge_structs(sparam.static_param,ctrl.dparam{1});
    fprintf('    task_function: %s\n', param.task_function);
    fprintf('    notes: %s\n', param.notes);
  end
  
  fprintf('    Number of tasks: %.0f, %.0f/%.0f/%.0f/%.0f C/R/Q/T, %.0f error, %.0f retries\n', ...
    numel(stats.cpu_time), sum(stats.job_status=='C'), sum(stats.job_status=='R'), sum(stats.job_status=='Q'), sum(stats.job_status=='T'), sum(stats.error_mask~=0), sum(stats.retries));
  fprintf('    Max CPU time: %.0f min\n', max(stats.cpu_time)/60);
  fprintf('    Max mem: %.0f MB\n', max(stats.mem)/1e6);
  fprintf('    Mean CPU time: %.0f min\n', mean(stats.cpu_time)/60);
  fprintf('    Mean mem: %.0f MB\n', mean(stats.mem)/1e6);
  
  if max(stats.cpu_time) > ctrl.cluster.max_time_per_job
    warning('Max task cpu time exceeds ctrl.cluster.max_time_per_job setting.');
  end
  
  % Update output
  ctrl_chain = ctrl;
  
elseif isnumeric(ctrl_chain)
  ctrl_chain = cluster_load_chain([],ctrl_chain);
  ctrl_chain = cluster_print_chain(ctrl_chain);
  
end

return
