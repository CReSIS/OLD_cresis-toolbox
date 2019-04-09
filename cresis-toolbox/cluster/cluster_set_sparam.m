function ctrl_chain = cluster_set_sparam(ctrl_chain,varargin)
% ctrl_chain = cluster_set_sparam(ctrl_chain,varargin)
%
% cluster_set_sparam(ctrl_chain,'argsin{1}.records.file.base_dir','/archive/');
% cluster_set_sparam(ctrl,'argsin{1}.records.file.base_dir','/archive/');
% cluster_set_sparam(ctrl,'argsin{1}',20);
%
% Author: John Paden

if iscell(ctrl_chain)
  %% Traverse chain lists
  for chain = 1:numel(ctrl_chain)
    fprintf('Chain %d\n', chain);
    for stage=1:numel(ctrl_chain{chain})
      fprintf('  Stage %d\n', stage);
      [ctrl_chain{chain}{stage}] = cluster_set_sparam(ctrl_chain{chain}{stage},varargin{:});
    end
  end
  
elseif isstruct(ctrl_chain)
  %% Set batch
  ctrl = ctrl_chain;
  
  % Load the sparam file
  static_in_fn = fullfile(ctrl.in_fn_dir,'static.mat');
  sparam = load(static_in_fn);
  
  for arg_idx = 1:2:numel(varargin)
    % Create and run the set commands
    if 1
      cmd = sprintf('sparam.static_param.%s = varargin{arg_idx+1};', varargin{arg_idx});
    else
      field_name_list = regexp(varargin{arg_idx}, '\.+', 'split');
      cmd = 'sparam.static_param';
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
  
  % Save the updated sparam
  save(static_in_fn,'-append','-struct','sparam','static_param');
end
