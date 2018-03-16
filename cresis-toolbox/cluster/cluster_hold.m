function cluster_hold(ctrl,hold_state)
% cluster_hold(ctrl,hold_state)
%
% Places or removes hold on specified batches.
%
% Inputs:
% ctrl = Several options which specify which batches to act on
%   1. Pass in a cluster batch ctrl structure (only needs "batch_dir"
%      defined)
%   2. A vector of batch ids to apply hold to
% hold_state = mode must be one of the following
%   0: removes hold
%   1: applies hold
%   []: if empty or undefined, the hold_state is toggled
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

if exist('ctrl','var') && isstruct(ctrl)
  %% This section actually does the placement/removal of holds
  if ~exist('hold_state','var') || isempty(hold_state)
    % When no hold state passed in, then toggle the hold state
    if exist(ctrl.hold_fn,'file')
      hold_state = 0;
    else
      hold_state = 1;
    end
  end
  
  if hold_state == 1
    fid = fopen(ctrl.hold_fn,'w');
    fclose(fid);
    
  elseif hold_state == 0
    if exist(ctrl.hold_fn,'file')
      delete(ctrl.hold_fn);
    end
    return
  end
  
else
  %% Handle case where a non-structure method of identifying the batch was used
  ctrls = cluster_get_batch_list;
  for batch_idx = 1:length(ctrls)
    if ~exist('ctrl','var') || any(ctrls{batch_idx}.batch_id == ctrl)
      if ~exist('hold_state','var') || isempty(hold_state)
        % When no hold state passed in, then toggle the hold state
        if exist(ctrls{batch_idx}.hold_fn,'file')
          hold_state = 0;
        else
          hold_state = 1;
        end
      end
      if hold_state == 1
        fprintf(' Placing hold on batch %d\n', ctrls{batch_idx}.batch_id);
      elseif hold_state == 0
        fprintf(' Removing hold on batch %d\n', ctrls{batch_idx}.batch_id);
      end
      cluster_hold(ctrls{batch_idx},hold_state);
    else
      fprintf(' Skipping %d\n', ctrls{batch_idx}.batch_id);
    end
  end
  return;
end
