function ctrl = cluster_save_dparam(ctrl)

if isfield(ctrl,'dparam')
  dynamic_in_fn = fullfile(ctrl.in_fn_dir,'dynamic.mat');
  if exist(dynamic_in_fn,'file')
    % Dynamic input arguments file already exists
    dynamic_param = load(dynamic_in_fn,'dparam');
    if ~isfield(dynamic_param,'dparam')
      dynamic_param.dparam = {};
    end
    % Merge file inputs with ctrl inputs
    for task_id = 1:length(ctrl.dparam)
      if ~isempty(ctrl.dparam{task_id})
        dynamic_param.dparam{task_id} = ctrl.dparam{task_id};
      end
    end
    save(dynamic_in_fn,ctrl.cluster.file_version,'-struct','dynamic_param','dparam');
  else
    % Dynamic input arguments file does not exist
    save(dynamic_in_fn,ctrl.cluster.file_version,'-struct','ctrl','dparam');
  end
end

