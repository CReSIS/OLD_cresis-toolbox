function success = layer_tracker_combine_task(param)
% success = layer_tracker_combine_task(param)
%
% Combines temporary layer data from layer_tracker_task.m and stores the
% output with opsCopyLayers. See run_layer_tracker_2D.m.

method_idx = 1;
track = param.layer_tracker.track(method_idx);

%% Create output directory path
if strcmp(param.layer_tracker.layer_params.source,'ops')
  tmp_out_fn_dir_dir = ct_filename_out(param,'ops','layer_tracker_tmp');
else
  tmp_out_fn_dir_dir = ct_filename_out(param,param.layer_tracker.layer_params.layerdata_source,'layer_tracker_tmp');
end

%% Load in all temporary files
gps_time = [];
twtt = [];
for frm = param.cmd.frms
  tmp_out_fn_name = sprintf('m%03d_%s.mat', method_idx, track.method);
  tmp_out_fn = fullfile(tmp_out_fn_dir_dir,sprintf('layer_tracker_%03d', frm),tmp_out_fn_name);
  fprintf('Loading %s (%s)\n', tmp_out_fn, datestr(now));
  tmp = load(tmp_out_fn);
  
  Nx = length(tmp.gps_time);
  gps_time(1,end+(1:Nx)) = tmp.gps_time;
  twtt(:,end+(1:Nx)) = tmp.twtt;
end

%% opsCopyLayer copy struct
param.layer_tracker.copy_param.gaps_fill.method = 'preserve_gaps';
copy_param = param.layer_tracker.copy_param;
copy_param.layer_source.existence_check = false;
copy_param.layer_source.source = 'custom';
copy_param.layer_source.gps_time = gps_time;
copy_param.layer_dest = param.layer_tracker.layer_params;
copy_param.layer_dest.existence_check = false;
if ~isfield(copy_param.layer_dest,'name') || isempty(copy_param.layer_dest.name)
  automated_name_en = true;
else
  automated_name_en = false;
end
%% Copy each layer
for layer_idx = 1:size(twtt,1)
  copy_param.layer_source.twtt = twtt(layer_idx,:);

  if automated_name_en
    % Create an automatic name
    switch track.method
      case 'lsm'
        if layer_idx <= length(track.lsm.storeIter)
          copy_param.layer_dest.name = sprintf('%s_%s_surface_%03d', ...
            track.name,track.method,layer_idx);
        else
          copy_param.layer_dest.name = sprintf('%s_%s_bottom_%03d', ...
            track.name,track.method,layer_idx-length(track.lsm.storeIter));
        end
      case {'mcmc','stereo'}
        if layer_idx == 1
          copy_param.layer_dest.name = sprintf('%s_%s_surface',track.name,track.method);
        else
          copy_param.layer_dest.name = sprintf('%s_%s_bottom',track.name,track.method);
        end
      case 'viterbi'
        copy_param.layer_dest.name = sprintf('%s_%s_bottom',track.name,track.method);
      otherwise
        copy_param.layer_dest.name = sprintf('%s_%s_surface',track.name,track.method);
    end
  end
  fprintf('opsCopyLayers %s %s (%s)\n', param.day_seg, copy_param.layer_dest.name, datestr(now));
  opsCopyLayers(param,copy_param);
end

fprintf('Done (%s)\n', datestr(now));
success = true;
