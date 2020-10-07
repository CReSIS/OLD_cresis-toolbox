function success = layer_tracker_combine_task(param)
% success = layer_tracker_combine_task(param)
%
% Combines temporary layer data from layer_tracker_task.m and stores the
% output with opsCopyLayers. See run_layer_tracker.m.
%
% Authors: Anjali Pare, John Paden
%
% See also: layer_tracker.m, layer_tracker_combine_task.m,
% layer_tracker_task.m, run_layer_tracker.m, run_layer_tracker_tune.m

%% Create output directory path
if strcmp(param.layer_tracker.layer_params.source,'ops')
  tmp_out_fn_dir_dir = ct_filename_out(param,'ops','layer_tracker_tmp');
else
  tmp_out_fn_dir_dir = ct_filename_out(param,param.layer_tracker.layer_params.layerdata_source,'layer_tracker_tmp');
end

layer_dest.name = [];
layer_source.twtt = [];
layer_source.gps_time = [];

%% Track loop
% Combine and Copy each tracking result
for track_idx = 1:length(param.layer_tracker.track)
  track = param.layer_tracker.track{track_idx};
  
  %% Track: Load in all temporary files
  gps_time = [];
  twtt = [];
  for frm = param.cmd.frms
    tmp_out_fn_name = sprintf('t%03d_%s.mat', track_idx, track.method);
    tmp_out_fn = fullfile(tmp_out_fn_dir_dir,sprintf('layer_tracker_%03d', frm),tmp_out_fn_name);
    fprintf('Loading %s (%s)\n', tmp_out_fn, datestr(now));
    tmp = load(tmp_out_fn);
    
    Nx = length(tmp.gps_time);
    gps_time(1,end+(1:Nx)) = tmp.gps_time;
    twtt(:,end+(1:Nx)) = tmp.twtt;
  end
  
  %% Track: opsCopyLayer copy struct
  param.layer_tracker.copy_param.gaps_fill.method = 'preserve_gaps';
  copy_param = param.layer_tracker.copy_param;
  copy_param.layer_source.existence_check = false;
  copy_param.layer_source.source = 'custom';
  copy_param.layer_source.gps_time = gps_time;
  copy_param.layer_dest = param.layer_tracker.layer_params;
  copy_param.layer_dest.existence_check = false;
  if ~isfield(param.layer_tracker.track{track_idx},'layer_names') || isempty(param.layer_tracker.track{track_idx}.layer_names)
    automated_name_en = true;
  else
    automated_name_en = false;
    layer_dest.name = param.layer_tracker.track{track_idx}.layer_names;
  end
  
  %% Track: Copy each layer
  for layer_idx = 1:size(twtt,1)
    if automated_name_en
      switch track.method
        case 'lsm'
          if layer_idx <= length(track.lsm.storeIter)
            layer_dest.name{end+1} = sprintf('%s_%s_surface_%03d', ...
              track.name,track.method,layer_idx);
          else
            layer_dest.name{end+1} = sprintf('%s_%s_bottom_%03d', ...
              track.name,track.method,layer_idx-length(track.lsm.storeIter));
          end
        case {'mcmc','stereo'}
          if layer_idx == 1
            layer_dest.name{end+1} = sprintf('%s_%s_surface',track.name,track.method);
          else
            layer_dest.name{end+1} = sprintf('%s_%s_bottom',track.name,track.method);
          end
        case 'viterbi'
          layer_dest.name{end+1} = sprintf('%s_%s_bottom',track.name,track.method);
        otherwise
          layer_dest.name{end+1} = sprintf('%s_%s_surface',track.name,track.method);
      end
    end
    
    layer_source.twtt{end+1} = twtt(layer_idx,:);
    layer_source.gps_time{end+1} = gps_time;
  end
  
end

copy_param.layer_source.gps_time = layer_source.gps_time;
copy_param.layer_source.twtt = layer_source.twtt;
copy_param.layer_dest.name = layer_dest.name;
fprintf('opsCopyLayers %s (%s)\n', param.day_seg, datestr(now));
opsCopyLayers(param,copy_param);

fprintf('Done (%s)\n', datestr(now));
success = true;
