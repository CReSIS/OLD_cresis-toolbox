function cluster_hold(ctrl_chain,hold_state,mode)
% cluster_hold(ctrl_chain,hold_state,mode)
%
% Places or removes hold on specified batches.
%
% Inputs:
% ctrl_chain = Several options to specify
%   1. Vector of chain IDs or batch IDs, if batch ID, then mode must be set to 'batch'
%   2. A ctrl structure identifying a batch
%   3. A chain (cell array of ctrl)
%   4. A list of chains (cell array of chains)
% hold_state = mode must be one of the following
%   0: removes hold
%   1: applies hold
%   []: if empty or undefined, the hold_state is toggled
% mode: Only used if ctrl_chain is an integer array, default is 'chain'.
%   Integer array will be treated as chains if mode is 'chain' and batches
%   if mode is 'batch'.
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

%% Input check
if nargin == 0 || isempty(ctrl_chain)
  answer = input('Are you sure you want to hold all cluster jobs? [y/N] ','s');
  if isempty(regexpi(answer,'y'))
    return
  end
  
  ctrl_chain = cluster_get_batch_list;
end

if ~exist('hold_state','var')
  hold_state = [];
end

%% Get a list of all batches
ctrls = cluster_get_batch_list;
ctrls_mask = logical(zeros(size(ctrls)));

%% Determine which batches
if isstruct(ctrl_chain)
  % This is an array of batch control structures
  for batch_idx = 1:length(ctrl_chain)
    batch_id = ctrl_chain(batch_idx).batch_id;
    for batch_idx = 1:length(ctrls)
      if ctrls{batch_idx}.batch_id == batch_id
        ctrls_mask(batch_idx) = true;
      end
    end
  end
elseif iscell(ctrl_chain)
  for chain_idx = 1:length(ctrl_chain)
    if ~iscell(ctrl_chain{chain_idx})
      % This is a control chain
      chain = {ctrl_chain{chain_idx}};
    else
      % This is a list of control chains
      chain = ctrl_chain{chain_idx};
    end
    for batch_idx = 1:length(chain)
      batch_id = chain{batch_idx}.batch_id;
      for batch_idx = 1:length(ctrls)
        if ctrls{batch_idx}.batch_id == batch_id
          ctrls_mask(batch_idx) = true;
        end
      end
    end
  end
  
elseif isnumeric(ctrl_chain)
  if nargin >= 3 && ~isempty(mode) && strcmpi(mode,'batch')
    % This is a list of batch IDs
    for idx = 1:length(ctrl_chain)
      batch_id = ctrl_chain(idx);
      for batch_idx = 1:length(ctrls)
        if ctrls{batch_idx}.batch_id == batch_id
          ctrls_mask(batch_idx) = true;
        end
      end
    end
    
  else
    % This is a list of chain IDs
    for idx = 1:length(ctrl_chain)
      try
        ctrl_chain = cluster_load_chain(ctrl_chain(idx));
      catch
        continue
      end
      cluster_hold(ctrl_chain,hold_state);
    end
    return;
  end
end
ctrls = ctrls(ctrls_mask);

%% Place hold on each batch
for ctrl_idx = 1:length(ctrls)
  ctrl = ctrls{ctrl_idx};
  
  if isempty(hold_state)
    % When no hold state passed in, then toggle the hold state
    if exist(ctrl.hold_fn,'file')
      hold_state = 0;
    else
      hold_state = 1;
    end
  end
  if hold_state == 1
    fprintf(' Placing hold on batch %d\n', ctrl.batch_id);
  elseif hold_state == 0
    fprintf(' Removing hold on batch %d\n', ctrl.batch_id);
  end
  
  if hold_state == 1
    fid = fopen(ctrl.hold_fn,'w');
    fclose(fid);
    
  elseif hold_state == 0
    if exist(ctrl.hold_fn,'file')
      delete(ctrl.hold_fn);
    end
  end
end

