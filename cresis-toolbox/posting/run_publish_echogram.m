function run_publish_echogram
%
% Example function which sets up parameters and calls publish_echogram
% and publish_map for several frames and annotates the maps
% for each frame.
%
% Author: John Paden
%
% See also: publish_echogram

% ========================================================================
% User Settings
% ========================================================================
fns = {};
fns{end+1} = 'Y:/ct_data/rds/2009_Greenland_TO/CSARP_mvdr/20090409_01/Data_20090409_01_005.mat';
lay_fns = {};
lay_fns{end+1} = 'Y:/ct_data/rds/2009_Greenland_TO/CSARP_layerData/20090409_01/Data_20090409_01_005.mat';

create_map = false;
map_param.type = 'contour';
map_param.location = 'greenland';
map_param.fig_hand = 1;
map_param.map_title = '';
map_param.decimate_seg = false;

echo_param.fig_hand = 2;
echo_param.num_x_tics = 4;
echo_param.depth = '[-500 2000]'; % depth range to plot (-500 to 3500 typical, deep antarctica -500 to 4000)
echo_param.elev_comp = true;
echo_param.plot_quality = 1;

create_outputs = false;
out_fn_dir = '~/';

% ========================================================================
% Automated Section
% ========================================================================

fprintf('==============================================================\n');
fprintf('run_publish_echogram\n');
fprintf('==============================================================\n');

% Setup map figure
if create_map
  map_info = publish_map('setup',map_param);
end

% Load data files and plot echograms, collect map info from each file
mdata = {};
lay = {};
X = [];
Y = [];
for fn_idx = 1:length(fns)
  mdata{fn_idx} = load(fns{fn_idx});
  if create_map
    [frame_X{fn_idx},frame_Y{fn_idx}] = projfwd(map_info.proj, ...
      mdata{fn_idx}.Latitude,mdata{fn_idx}.Longitude);
    X = cat(2,X,frame_X{fn_idx});
    Y = cat(2,Y,frame_Y{fn_idx});
  end
  lay{fn_idx} = load(lay_fns{fn_idx});
  [fn_dir fn_name] = fileparts(fns{fn_idx});
  echo_param.frm_id = fn_name(end-14:end);
  echo_info = publish_echogram(echo_param,mdata{fn_idx},lay{fn_idx});
  echo_fn = fullfile(out_fn_dir,sprintf('Data_%s.fig',echo_param.frm_id));
  %saveas(echo_param.fig_hand, echo_fn);
end

if create_map
  % Plot the map for each frame
  map_param.day_seg_x = X;
  map_param.day_seg_y = Y;
  map_param.frame_X = frame_X;
  map_param.frame_Y = frame_Y;
  map_info = publish_map('plot',map_param,map_info);
  legend([map_info.h_frm{1} map_info.h_start{1}], ...
    'Frame', 'Start', 'Location', 'Southeast')
  
  % Add text annotations to regional map for each segment
  for fn_idx = 1:length(fns)
    axes(map_info.ah_region);
    h = text(frame_X{fn_idx}(1)/1000,frame_Y{fn_idx}(1)/1000, ...
      sprintf('%d', fn_idx), 'FontSize', 14);
    [fn_dir fn_name] = fileparts(fns{fn_idx});
    % Print the frame ID out
    fprintf('%d: %s\n', fn_idx, fn_name(end-14:end));
  end
  
  if create_outputs
    % Save outputs
    map_fn = fullfile(out_fn_dir,sprintf('Map_%s.fig',echo_param.frm_id));
    saveas(map_param.fig_hand, map_fn);
  end
end

if create_outputs
  % Save outputs
  echo_fn = fullfile(out_fn_dir,sprintf('Echo_%s.fig',echo_param.frm_id));
  saveas(echo_param.fig_hand, echo_fn);
end

return;
