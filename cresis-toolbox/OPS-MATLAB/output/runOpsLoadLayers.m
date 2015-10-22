% script runOpsLoadLayers.m
%
% Runs opsLoadLayers.m

if 1
  %% User Settings
  params = read_param_xls(ct_filename_param('rds_param_2008_Greenland_TO.xls'),'20080715_05','post');
  
  layer_params = [];
  
  idx = 0;
  ref_idx = 2;
  
  idx = idx + 1;
  layer_params(idx).name = 'atm';
  layer_params(idx).source = 'ops';
  layer_params(idx).echogram_source = 'qlook';
  layer_params(idx).layerdata_source = 'layerData';
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'ops';
  layer_params(idx).echogram_source = 'qlook';
  layer_params(idx).layerdata_source = 'layerData';
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'records';
  layer_params(idx).echogram_source = 'qlook';
  layer_params(idx).layerdata_source = 'layerData';
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'lidar';
  layer_params(idx).echogram_source = 'qlook';
  layer_params(idx).layerdata_source = 'layerData';
  
  %% Automated Section
  
  %% Load each of the day segments
  % =====================================================================
  layers = {};
  day_seg = {};
  for param_idx = 1:length(params)
    param = params(param_idx);
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
      continue;
    end
    param = merge_structs(param,gRadar);
    fprintf('opsLoadLayers %s\n', param.day_seg);
    layers{end+1} = opsLoadLayers(param,layer_params);
    day_seg{end+1} = param.day_seg;
  end
  
  for seg_idx = 1:length(layers)
    % Interpolate onto reference
    lay_idxs = [1:ref_idx-1 ref_idx+1:length(layers{seg_idx})];
    
    layers{seg_idx}(ref_idx).twtt_ref = layers{seg_idx}(ref_idx).twtt;
    
    master = [];
    master.GPS_time = layers{seg_idx}(ref_idx).gps_time;
    master.Latitude = layers{seg_idx}(ref_idx).lat;
    master.Longitude = layers{seg_idx}(ref_idx).lon;
    master.Elevation = layers{seg_idx}(ref_idx).elev;
    for lay_idx = lay_idxs
      ops_layer = [];
      ops_layer{1}.gps_time = layers{seg_idx}(lay_idx).gps_time;
      ops_layer{1}.type = layers{seg_idx}(lay_idx).type;
      ops_layer{1}.quality = layers{seg_idx}(lay_idx).quality;
      ops_layer{1}.twtt = layers{seg_idx}(lay_idx).twtt;
      ops_layer{1}.type(isnan(ops_layer{1}.type)) = 2;
      ops_layer{1}.quality(isnan(ops_layer{1}.quality)) = 1;
      lay = opsInterpLayersToMasterGPSTime(master,ops_layer,[300 60]);
      layers{seg_idx}(lay_idx).twtt_ref = lay.layerData{1}.value{2}.data;
    end
  end

end

status = {};
for seg_idx = 1:length(layers)
  
  figure(1); clf;
  plot(layers{seg_idx}(2).gps_time, layers{seg_idx}(2).twtt, '.')
  hold on;
  plot(layers{seg_idx}(1).gps_time, layers{seg_idx}(1).twtt, 'r.')
  plot(layers{seg_idx}(3).gps_time, layers{seg_idx}(3).twtt, 'g.')
  plot(layers{seg_idx}(4).gps_time, layers{seg_idx}(4).twtt, 'k.')
  grid on
  
  figure(2); clf;
  plot(layers{seg_idx}(ref_idx).gps_time, layers{seg_idx}(ref_idx).twtt, '.')
  for lay_idx = lay_idxs
    hold on;
    plot(layers{seg_idx}(ref_idx).gps_time, layers{seg_idx}(lay_idx).twtt_ref, 'r.')
    grid on
    hold off;
  end
    
  fprintf('%s', day_seg{seg_idx});
  figure(3); clf;
  plot_modes = {'r.' 'b.' 'g.' 'k.'};
  for lay_idx = lay_idxs
    plot(layers{seg_idx}(ref_idx).gps_time, layers{seg_idx}(lay_idx).twtt_ref - layers{seg_idx}(ref_idx).twtt, plot_modes{lay_idx})
    mean_offset = nanmean(layers{seg_idx}(lay_idx).twtt_ref - layers{seg_idx}(ref_idx).twtt);
    fprintf('\t%.3f\t%.3f\t%.3f', ...
      1e6*mean_offset, ...
      1e6*nanmedian(layers{seg_idx}(lay_idx).twtt_ref - layers{seg_idx}(ref_idx).twtt), ...
      1e6*nanstd(layers{seg_idx}(lay_idx).twtt_ref - layers{seg_idx}(ref_idx).twtt), ...
      1e6*nanmax(abs(layers{seg_idx}(lay_idx).twtt_ref - layers{seg_idx}(ref_idx).twtt - mean_offset)));
    hold on;
    grid on
  end
  hold off;
  fprintf('\n');
  
  %status{seg_idx} = input(sprintf('%s: ', day_seg{seg_idx}), 's');

  pause;
  
end

layers_fn = ct_filename_tmp(rmfield(param,'day_seg'),'','surf_layers','surf_layers.mat')
layers_fn_dir = fileparts(layers_fn);
if ~exist(layers_fn_dir,'dir')
  mkdir(layers_fn_dir);
end
save(layers_fn,'-v7.3','day_seg','layers')



