function ctrl_chain = cluster_set_dparam(ctrl_chain,varargin)
% ctrl_chain = cluster_set_dparam(ctrl_chain,varargin)
%
% cluster_set_dparam(ctrl_chain,'argsin{1}.cpu_time',100);
% cluster_set_dparam(ctrl,'argsin{1}',40);
% cluster_set_dparam(ctrl,'cpu_time',100);
%
% Author: John Paden

if iscell(ctrl_chain)
  %% Traverse chain lists
  for chain = 1:numel(ctrl_chain)
    fprintf('Chain %d\n', chain);
    for stage=1:numel(ctrl_chain{chain})
      fprintf('  Stage %d\n', stage);
      [ctrl_chain{chain}{stage}] = cluster_set_dparam(ctrl_chain{chain}{stage},varargin{:});
    end
  end
  
elseif isstruct(ctrl_chain)
  %% Set batch
  ctrl = ctrl_chain;
  
  % Load the dparam file
  dynamic_in_fn = fullfile(ctrl.in_fn_dir,'dynamic.mat');
  dparam = load(dynamic_in_fn);
  
  for arg_idx = 1:2:numel(varargin)
    % Create and run the set commands
    if 1
      cmd = sprintf('dparam.dparam{task_id}.%s = varargin{arg_idx+1};', varargin{arg_idx});
    else
      field_name_list = regexp(varargin{arg_idx}, '\.+', 'split');
      cmd = 'dparam.dparam{task_id}';
      for idx=1:length(field_name_list)
        cmd = cat(2,cmd,sprintf('.(field_name_list{%d})',idx));
      end
      cmd = cat(2,cmd,sprintf(' = varargin{%d};',arg_idx+1));
    end
    try
      for task_id = 1:ctrl.task_id
        eval(cmd)
      end
    catch ME
      warning(ME.getReport);
    end
  end
  
  % Save the updated dparam file
  save(dynamic_in_fn,'-append','-struct','dparam','dparam');
end
