function active_stage = cluster_chain_stage(ctrl_chain)
% active_stage = cluster_chain_stage(ctrl_chain)
%
% Checks to see if a cluster chain is done running.
%
% INPUTS:
% ctrl_chain: cluster chain of control/batch structures
%
% OUTPUTS:
% active_stage: vector the same size as ctrl_chain representing the stage of
%   each chain.
%   Values from 1 to numel(chain) represent the active stage.
%   Value of inf means the chain is complete.
%   Value of -inf means chain stopped (usually because the chain is
%   out of retries due to errors).
%
% EXAMPLES:
% active_stage = cluster_chain_stage(ctrl_chain)
%
% Authors: John Paden
%
% See also: cluster_batch_list cluster_cleanup cluster_compile ...
%   cluster_create_task cluster_hold cluster_job_list cluster_job_status ...
%   cluster_new_batch cluster_print cluster_rerun

%% Traverse chain list
active_stage = ones(size(ctrl_chain));
for chain = 1:numel(ctrl_chain)
  for stage = 1:numel(ctrl_chain{chain})
    % 1. There is at least one batch left to run in this chain
    ctrl = ctrl_chain{chain}{stage};
    
    if all(ctrl.job_status=='C')
      if ~any(ctrl.error_mask)
        active_stage(chain) = active_stage(chain) + 1;
        if active_stage(chain) > numel(ctrl_chain{chain})
          % Chain is complete
          active_stage(chain) = inf;
        end
      elseif all(ctrl.retries >= ctrl.cluster.max_retries)
        % Chain is stopped
        active_stage(chain) = -inf;
        break;
      end
    else
      % Chain is working on this stage
      break;
    end
    
  end
end
