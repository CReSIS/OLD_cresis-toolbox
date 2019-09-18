function ctrl_chain = cluster_set_chain(ctrl_chain,varargin)
% ctrl_chain = cluster_set_chain(ctrl_chain,varargin)
%
% ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.cpu_time_mult',2,'cluster.mem_mult',1);
% ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.type','debug');
% ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.desired_time_per_job',24*60*60);

if isnumeric(ctrl_chain)
  ctrl_chain = cluster_load_chain(ctrl_chain);
end

if iscell(ctrl_chain)
  %% Traverse chain lists
  for chain = 1:numel(ctrl_chain)
    fprintf('Chain %d\n', chain);
    for stage=1:numel(ctrl_chain{chain})
      fprintf('  Stage %d\n', stage);
      [ctrl_chain{chain}{stage}] = cluster_set_chain(ctrl_chain{chain}{stage},varargin{:});
    end
  end
  
elseif isstruct(ctrl_chain)
  %% Set batch
  ctrl = ctrl_chain;
  
  for arg_idx = 1:2:numel(varargin)
    % Create and run the set commands
    if 1
      cmd = sprintf('ctrl.%s = varargin{arg_idx+1};', varargin{arg_idx});
    else
      field_name_list = regexp(varargin{arg_idx}, '\.+', 'split');
      cmd = 'ctrl';
      for idx=1:length(field_name_list)
        cmd = cat(2,cmd,sprintf('.(field_name_list{%d})',idx));
      end
      cmd = cat(2,cmd,sprintf(' = varargin{%d};',arg_idx+1));
    end
    try
      eval(cmd)
    catch ME
      warning(ME.getReport);
    end
  end
  
  % Update output
  ctrl_chain = ctrl;
end
