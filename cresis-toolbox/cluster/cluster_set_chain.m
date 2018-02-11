function ctrl_chain = cluster_set_chain(ctrl_chain,varargin)
% ctrl_chain = cluster_set_chain(ctrl_chain,varargin)
%
% ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.cpu_time_mult',2,'cluster.mem_mult',1);

if iscell(ctrl_chain)
  %% Traverse chain list
  for chain = 1:numel(ctrl_chain)
    fprintf('Chain %d\n', chain);
    for stage=1:numel(ctrl_chain{chain})
      fprintf('  Stage %d\n', stage);
      [ctrl_chain{chain}{stage}] = cluster_set_chain(ctrl_chain{chain}{stage},varargin{:});
    end
  end
  
elseif isstruct(ctrl_chain)
  ctrl = ctrl_chain;
  
  for arg_idx = 1:2:numel(varargin)
    field_name_list = regexp(varargin{arg_idx}, '\.+', 'split');
    cmd = 'ctrl';
    for idx=1:length(field_name_list)
      cmd = cat(2,cmd,sprintf('.(field_name_list{%d})',idx));
    end
    cmd = cat(2,cmd,sprintf(' = varargin{%d};',arg_idx+1));
    try
      eval(cmd)
    catch ME
      warning(ME.getReport);
    end
  end
  
  % Update output
  ctrl_chain = ctrl;
end

return


end
