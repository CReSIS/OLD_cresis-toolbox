% script runOpsLoadLayers.m
%
% Example script for running opsLoadLayers.m. Demonstrates a few of the
% most common operations to be performed with opsLoadLayers and supports
% interpolation of various layer sources for better comparison.
%
% Authors: John Paden

% =====================================================================
%% User Settings
% =====================================================================

params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'));
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20120330_04');
params = ct_set_params(params,'cmd.frms',[2]);

layers = {};
layer_params = [];
idx = 0;

%% Example Inputs (just choose one)
if 0
  %% Load layers from layerData file version 1  
  ref_idx = 1;
  idx = idx + 1;
  layer_params(idx).layers2load = []; % leave empty for all;
  layer_params(idx).name = 'layers'; % surface, bottom or layers for all.
  layer_params(idx).source = 'layerData_ver1';
  layer_params(idx).layerdata_source = 'layerData_ver1';
  
elseif 0
  %% Compare echogram surface and a MacGregor's layer in ops  
  ref_idx = 1;
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'echogram';
  layer_params(idx).echogram_source = 'CSARP_post/qlook';
  idx = idx + 1;
  layer_params(idx).name = 'layer_1950-001988';
  layer_params(idx).source = 'ops';

elseif 0
  %% Load a single layer from the echogram
  ref_idx = 1;
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'echogram';
  layer_params(idx).echogram_source = 'qlook';

elseif 0
  %% Load a single layer from the layerData file
  ref_idx = 1;
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'layerData';
  layer_params(idx).layerdata_source = 'layerData';
elseif 1
  %% Load multiple snow layers from the layerData file (koenig's snow layers)
  %% Currently only work for a single frame
  ref_idx = 1;
  idx = idx + 1;
  layer_params(idx).layers2load = []; % leave empty for all;
  layer_params(idx).name = 'layers';  % surface, bottom or layers for all.
  layer_params(idx).source = 'layerData_koenig';
  layer_params(idx).layerdata_source = 'layerData_koenig';
elseif 0
  %% Compare echogram, layerData, and records
  ref_idx = 1;
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'echogram';
  layer_params(idx).echogram_source = 'qlook';
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'layerData';
  layer_params(idx).layerdata_source = 'layerData';
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'records';
 
elseif 0
  %% Compare echogram and custom layers in layerData
  ref_idx = 1;
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'lidar';
  layer_params(idx).lidar_source = 'awi';
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'layerData';
  layer_params(idx).layerdata_source = 'layerData';
  layer_params(idx).twtt_offset = 0;
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'echogram';
  layer_params(idx).echogram_source = 'standard';
  layer_params(idx).twtt_offset = 0;% 1.18603e-07;
  idx = idx + 1;
  layer_params(idx).name = 'GIMP';
  layer_params(idx).source = 'layerData';
  layer_params(idx).layerdata_source = 'layerData';
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'records';
  layer_params(idx).twtt_offset = 0;
  
elseif 0
  %% Compare OPS surface (records and OPS) to ATM data (raw and OPS)
  ref_idx = 2;
  idx = idx + 1;
  layer_params(idx).name = 'atm';
  layer_params(idx).source = 'ops';
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'ops';
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'records';
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'lidar';
end

% =====================================================================
%% Automated Section
% =====================================================================

global gRadar;

%% Load each of the day segments
layers = {};
day_seg = {};
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) ...
      || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  param = merge_structs(param,gRadar);
  
  fprintf('opsLoadLayers %s\n', param.day_seg);
  layers{end+1} = opsLoadLayers(param,layer_params);
  day_seg{end+1} = param.day_seg;
end

% =====================================================================
%% Example Section
% =====================================================================

if 1
  %% Plot all layers versus relative time using layerData of version 1
  lay_idxs = 1:length(layers{ref_idx}.twtt);
  figure(1);clf
  for lyr_idx = lay_idxs
    figure(1);hold on;plot(layers{ref_idx}.gps_time-layers{ref_idx}.gps_time(1),layers{ref_idx}.twtt{lyr_idx}*1e6);
  end
  xlabel('relative time (s)');
  ylabel('twtt(\mus)')
  title(day_seg{ref_idx},'Interpreter','none');
end

if 0
  %% Interpolate all layers onto a common reference (ref_idx)
  for seg_idx = 1:length(layers)
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

if 0
  %% Simple plot of all layers versus time
  % ref_color: Reference plot color
  ref_color = 'k.';
  % Layer indexes for comparison layers
  lay_idxs = [1:ref_idx-1 ref_idx+1:length(layers{seg_idx})];
  % plot_modes: Different colors used to plot
  plot_modes = {'r.','g.','c.','b.','m.'};
  
  for seg_idx = 1:length(layers)
    if ~isfield(layer_params(ref_idx),'twtt_offset') ...
        || isempty(layer_params(ref_idx).twtt_offset)
      layer_params(ref_idx).twtt_offset = 0;
    end
    figure(1); clf;
    h_axes = axes('Parent',1);
    h_plot = plot(layers{seg_idx}(ref_idx).gps_time, ...
      layers{seg_idx}(ref_idx).twtt + layer_params(ref_idx).twtt_offset, ...
      'k.','Parent',h_axes);
    title(sprintf('%s', day_seg{seg_idx}),'Interpreter','none','Parent',h_axes);
    legend_strs = {sprintf('Ref %d', ref_idx)};
    hold(h_axes,'on');
    grid(h_axes,'on');
    xlabel('GPS time (sec)','Parent',h_axes);
    ylabel('TWTT (us)','Parent',h_axes);
    
    figure(2); clf(2);
    h_axes_comp = axes('Parent',2);
    h_plot_comp = plot(layers{seg_idx}(ref_idx).gps_time, ...
      layers{seg_idx}(ref_idx).twtt_ref + layer_params(ref_idx).twtt_offset ...
      - (layers{seg_idx}(ref_idx).twtt_ref + layer_params(ref_idx).twtt_offset), ...
      'k.','Parent',h_axes_comp);
    title(sprintf('%s', day_seg{seg_idx}),'Interpreter','none','Parent',h_axes_comp);
    hold(h_axes_comp,'on');
    grid(h_axes_comp,'on');
    xlabel('GPS time (sec)','Parent',h_axes_comp);
    ylabel('TWTT Difference (us)','Parent',h_axes_comp);
    
    fprintf('Layer Index\tMedian Offset\tMean Offset\n');
    
    for lay_idx = lay_idxs
      if ~isfield(layer_params(lay_idx),'twtt_offset') ...
          || isempty(layer_params(lay_idx).twtt_offset)
        layer_params(lay_idx).twtt_offset = 0;
      end
      
      h_plot(end+1) = plot(layers{seg_idx}(lay_idx).gps_time, ...
        layers{seg_idx}(lay_idx).twtt + layer_params(lay_idx).twtt_offset, ...
        plot_modes{mod(lay_idx-1,length(plot_modes))+1},'Parent',h_axes);
      legend_strs{end+1} = sprintf('Layer %d', lay_idx);

      h_plot_comp(end+1) = plot(layers{seg_idx}(ref_idx).gps_time, ...
        layers{seg_idx}(lay_idx).twtt_ref + layer_params(lay_idx).twtt_offset ...
        - (layers{seg_idx}(ref_idx).twtt_ref + layer_params(ref_idx).twtt_offset), ...
        plot_modes{mod(lay_idx-1,length(plot_modes))+1},'Parent',h_axes_comp);
      
      fprintf('%11d\t%13g\t%11g\n', lay_idx, ...
        nanmedian(layers{seg_idx}(lay_idx).twtt_ref + layer_params(lay_idx).twtt_offset ...
        - (layers{seg_idx}(ref_idx).twtt_ref + layer_params(ref_idx).twtt_offset)), ...
        nanmean(layers{seg_idx}(lay_idx).twtt_ref + layer_params(lay_idx).twtt_offset ...
        - (layers{seg_idx}(ref_idx).twtt_ref + layer_params(ref_idx).twtt_offset)));
    end
    legend(h_plot, legend_strs);
    legend(h_plot_comp, legend_strs);
    linkaxes([h_axes_comp h_axes],'x');
    
    if seg_idx ~= length(layers)
      pause;
    end
  end
end

if 0
  %% Compare N layers
  % ref_color: Reference plot color
  ref_color = 'k.';
  % Layer indexes for comparison layers
  lay_idxs = [1:ref_idx-1 ref_idx+1:length(layers{seg_idx})];
  % plot_modes: Different colors used to plot
  plot_modes = {'r.','g.','c.','b.','m.'};
  
  status = {};
  for seg_idx = 1:length(layers)

    figure(1); clf;
    plot(layers{seg_idx}(2).gps_time, layers{seg_idx}(ref_idx).twtt, 'k.');
    title(sprintf('%s', day_seg{seg_idx}),'Interpreter','none');
    hold on;
    for lay_idx = lay_idxs
      plot(layers{seg_idx}(lay_idx).gps_time, layers{seg_idx}(lay_idx).twtt, ...
        plot_modes{mod(lay_idx-1,length(plot_modes))+1});
    end
    h_axis = gca;
    grid on
    xlabel('GPS time (sec)');
    ylabel('TWTT ({\mieu}s)');
    
    figure(2); clf;
    plot(layers{seg_idx}(ref_idx).gps_time, layers{seg_idx}(ref_idx).twtt, '.');
    title(sprintf('%s', day_seg{seg_idx}),'Interpreter','none');
    for lay_idx = lay_idxs
      hold on;
      plot(layers{seg_idx}(ref_idx).gps_time, layers{seg_idx}(lay_idx).twtt_ref, 'r.')
      grid on
      hold off;
    end
    h_axis(2) = gca;
    xlabel('GPS time (sec)');
    ylabel('TWTT ({\mieu}s)');
    
    figure(3); clf;
    for lay_idx = lay_idxs
      plot(layers{seg_idx}(ref_idx).gps_time, layers{seg_idx}(lay_idx).twtt_ref ...
        - layers{seg_idx}(ref_idx).twtt, plot_modes{lay_idx});
      title(sprintf('%s', day_seg{seg_idx}),'Interpreter','none');
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
    h_axis(3) = gca;
    xlabel('GPS time (sec)');
    ylabel('TWTT ({\mieu}s)');
    fprintf('\n');

    linkaxes(h_axis,'x');
    
    if seg_idx ~= length(layers)
      pause;
    end
    
  end
end

if 0
  %% Example for saving output
  layers_fn = ct_filename_tmp(rmfield(param,'day_seg'),'','surf_layers','surf_layers.mat')
  layers_fn_dir = fileparts(layers_fn);
  if ~exist(layers_fn_dir,'dir')
    mkdir(layers_fn_dir);
  end
  save(layers_fn,'-v7.3','day_seg','layers')
end


