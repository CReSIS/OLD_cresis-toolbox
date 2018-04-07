function [ctrl_chain,chain_fn] = cluster_load_chain(chain_fn)
% [ctrl_chain,chain_fn] = cluster_load_chain(chain_fn)
%
% Load a control chain from a file. Usually saved by cluster_save_chain.m.
%
% INPUTS
% chain_fn: A control chain ID (positive integer) or filename string. If a number,
%   gRadar.cluster.data_location determines the directory to load cluster chains and batches
%
% OUTPUTS
% ctrl_chain: control chain loaded from the control chain file specified by
%   the input
%
% EXAMPLES
% [ctrl_chain,chain_fn] = cluster_load_chain(1);
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list,
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

if ischar(chain_fn)
  load(chain_fn);
  return;
  
elseif isnumeric(chain_fn)
  load_chain_id = chain_fn(1);
  
  global gRadar;
  
  chain_fns = get_filenames(gRadar.cluster.data_location,'chain_','','',struct('type','f'));
  
  for chain_idx = 1:length(chain_fns)
    chain_fn = chain_fns{chain_idx};
    [~,chain_fn_name] = fileparts(chain_fn);
    chain_id = strtok(chain_fn_name(7:end),'_');
    chain_id = str2double(chain_id);
    if chain_id == load_chain_id
      load(chain_fn);
      return;
    end
  end
  
  error('Chain ID %d not found.', load_chain_id);
end
