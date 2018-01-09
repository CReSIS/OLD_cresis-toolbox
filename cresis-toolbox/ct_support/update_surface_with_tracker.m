% script update_surface_with_tracker
%
% Updates the Surface variable in an echogram using one of three surface
% detections methods.
%
% surf: Parameters used to control the tracker. Only general parameters
%   used within update_surface_with_tracker are included here. For more
%   details about tracker specific parameters, look at the help for:
%   tracker_snake_simple, tracker_snake_manual_gui, tracker_threshold,
%   tracker_max, tracker_snake_simple
% .en: must be true for update_surface_with_tracker to run on a segment
% .medfilt: scalar which specifies the length of a median filter to apply
%   to the tracked data. Similar to medfilt1, but it handles edge
%   conditions and allows for a threshold parameter to be set.
% .medfilt_threshold: scalar used with medfilt. Median filter will only
%   update in the difference of the median filter is more than this
%   threshold. Default is 0 which means every pixel is updated.
% .feedthru: all values in the radar image which exceed the power mask set
%   by the time and power_dB vectors will be set to zero power. To get
%   these values a typical procedure is: 
%     plot(mdata.Time, lp(mean(mdata.Data,2)))
%     [surf.feedthru.time,surf.feedthru.power_dB] = ginput; % Press enter to end the input
%   .time: N length vector of two way travel times
%   .power_dB: N length vector of power in dB
% .max_diff: used by tracker routines (optional, default is inf)
% .min_bin: used by tracker routines
% .max_bin: used by tracker routines (optional, default is not used)
% .manual: used by tracker routines (optional, default is false)
% .init.method: used by tracker routines (optional, default is not used)
% .init.lidar_source: used by tracker routines (optional, default is not used)
%
% For more details about tracker specific parameters:
% tracker_snake_simple, tracker_snake_manual_gui, tracker_threshold,
% tracker_max, tracker_snake_simple
%
% Example:
%   See run_update_surface_with_tracker.m to run.
%
% Author: John Paden

%% Automated Section
% ----------------------------------------------------------------------
physical_constants;
global gRadar;

if ~exist('debug_time_guard','var') || isempty(debug_time_guard)
  debug_time_guard = 50e-9;
end

% gimp_geoid_loaded = false;

for param_idx = 1:length(params)
  param = params(param_idx);
  param = merge_structs(gRadar,param);
  
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  orig_surf = merge_structs(param.get_heights.surf,surf_override);
  
  if ~orig_surf.en
    continue;
  end
  
  if ~isfield(orig_surf,'max_diff')
    orig_surf.max_diff = inf;
  end
  
  if ~isfield(orig_surf,'medfilt_threshold')
    orig_surf.medfilt_threshold = 0;
  end
  
  fprintf('Updating surface %s (%s)\n', param.day_seg, datestr(now,'HH:MM:SS'));
  
  %% Load in GIMP and Geoid
  load_surface_land_dems = false;
  if isfield(orig_surf,'init') && strcmpi(orig_surf.init.method,'dem') ...
      && (~exist('load_surface_land_dems_finished','var') || ~load_surface_land_dems_finished)
    load_surface_land_dems = true;
  end
  
  if load_surface_land_dems
    sea_surface.fn = ct_filename_gis([],fullfile('world','dtu_meansealevel','DTU10MSS_1min.nc'));
    sea_surface.lat = ncread(sea_surface.fn,'lat');
    sea_surface.lon = ncread(sea_surface.fn,'lon');
    sea_surface.elev = ncread(sea_surface.fn,'mss').';
    if 0
      sea_surface.fn = ct_filename_gis([],fullfile('world','egm96_geoid','WW15MGH.DAC'));
      [sea_surface.lat,sea_surface.lon,sea_surface.elev] = egm96_loader(sea_surface.fn);
    end
    
    if strcmpi(params(end).post.ops.location,'arctic')
      land_surface.fn = ct_filename_gis([],'greenland/DEM/GIMP/gimpdem_90m.tif');
    elseif strcmpi(params(end).post.ops.location,'antarctic')
      land_surface.fn = ct_filename_gis([],'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_surface.tif');
    end
    [land_surface.dem, land_surface.R, tmp] = geotiffread(land_surface.fn);
    land_surface.dem = double(land_surface.dem);
    land_surface.dem(land_surface.dem == 32767) = NaN;
    land_surface.proj = geotiffinfo(land_surface.fn);
    
    load_surface_land_dems_finished = true;
    
    if 0
      % Debug Plot
      figure(1); clf;
      land_surface.x = land_surface.R(3,1) + land_surface.R(2,1)*(0:size(land_surface.dem,2)-1);
      land_surface.y = land_surface.R(3,2) + land_surface.R(1,2)*(0:size(land_surface.dem,1)-1);
      imagesc(land_surface.x,land_surface.y,land_surface.dem)
      set(gca,'YDir','normal');
    end
  end

  %% Load "frames" for this segment
  load(ct_filename_support(param,'','frames'));
  
  %% Determine valid frames to process
  if isempty(param.cmd.frms)
    param.cmd.frms = 1:length(frames.frame_idxs);
  end
  % Remove frames that do not exist from param.cmd.frms list
  [valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
  if length(valid_frms) ~= length(param.cmd.frms)
    bad_mask = ones(size(param.cmd.frms));
    bad_mask(keep_idxs) = 0;
    warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
      param.cmd.frms(find(bad_mask,1)));
    param.cmd.frms = valid_frms;
  end
  
  %% Get all the frames for this segment
  if any(strcmpi(save_sources,'ops'))
    opsAuthenticate(param,false);
    sys = ct_output_dir(param.radar_name);
    ops_param = struct('properties',[]);
    ops_param.properties.season = param.season_name;
    ops_param.properties.segment = param.day_seg;
    [status,ops_seg_data] = opsGetSegmentInfo(sys,ops_param);
  end
  
  %% Load LIDAR data if exists
  if isfield(orig_surf,'init') && isfield(orig_surf.init,'lidar_source')
    layer_params = [];
    lay_idx = 1;
    layer_params(lay_idx).name = 'surface';
    layer_params(lay_idx).source = 'lidar';
    layer_params(lay_idx).lidar_source = orig_surf.init.lidar_source;
    
    layers = opsLoadLayers(param,layer_params);
    
    % Ensure that layer gps times are monotonically increasing
    lay_idx = 1;
    layers_fieldnames = fieldnames(layers(lay_idx));
    [~,unique_idxs] = unique(layers(lay_idx).gps_time);
    for field_idx = 1:length(layers_fieldnames)-1
      if ~isempty(layers(lay_idx).(layers_fieldnames{field_idx}))
        layers(lay_idx).(layers_fieldnames{field_idx}) = layers(lay_idx).(layers_fieldnames{field_idx})(unique_idxs);
      end
    end
  end

  %% Track each of the frames
  for frm_idx = 1:length(param.cmd.frms)
    frm = param.cmd.frms(frm_idx);
    
    if ~exist('echogram_img','var')
      echogram_img = 0;
    end
    if echogram_img == 0
      data_fn = fullfile(ct_filename_out(param,echogram_source,''), ...
        sprintf('Data_%s_%03d.mat', param.day_seg, frm));
    else
      data_fn = fullfile(ct_filename_out(param,echogram_source,''), ...
        sprintf('Data_img_%02d_%s_%03d.mat', echogram_img, param.day_seg, frm));
    end
    fprintf('  Updating %s (%s)\n', data_fn, datestr(now,'HH:MM:SS'));
    
    if ~exist(data_fn)
      fprintf('    Missing file\n');
      continue;
    end
    
    %% Load echogram data
    mdata = load(data_fn);
    
    %% Update surf structure for this particular echogram
    %   - min_bin field is specified as two way travel time, but needs to
    %   be converted to range bins
    surf = orig_surf;
    surf.min_bin = find(mdata.Time > orig_surf.min_bin, 1);
    if isfield(orig_surf,'max_bin') && ~isempty(orig_surf.max_bin)
      surf.max_bin = find(mdata.Time > orig_surf.max_bin, 1);
    else
      surf.max_bin = inf;
    end
    dt = mdata.Time(2) - mdata.Time(1);

    if ~isfield(orig_surf,'manual')
      orig_surf.manual = false;
    end
    
    if isfield(orig_surf,'feedthru')
      %% Optional feed through removal
      
      % Interpolate feed through power levels on to data time axis
      feedthru_threshold = interp1(orig_surf.feedthru.time,orig_surf.feedthru.power_dB,mdata.Time);
      feedthru_threshold = interp_finite(feedthru_threshold,-inf);
      
      % Set all data to zero that does not exceed the feed through
      % threshold power
      for rline=1:size(mdata.Data,2)
        mdata.Data(:,rline) = mdata.Data(:,rline) .* (lp(mdata.Data(:,rline)) > feedthru_threshold);
      end
    end
    
    %% Initialize surf.dem for dem and lidar methods
    surf.dem = nan(size(mdata.GPS_time));

    %% Interpolate GIMP and Geoid
    if isfield(orig_surf,'init') && strcmpi(orig_surf.init.method,'dem')
      mdata.sea_dem = interp2(sea_surface.lon,sea_surface.lat,sea_surface.elev,mod(mdata.Longitude,360),mdata.Latitude);
      [mdata.x,mdata.y] = projfwd(land_surface.proj,mdata.Latitude,mdata.Longitude);
      mdata.land_dem = interp2(land_surface.dem,(mdata.x-land_surface.R(3,1))/land_surface.R(2,1)+1, ...
        (mdata.y-land_surface.R(3,2))/land_surface.R(1,2)+1);
      
      % Merge GIMP and Geoid
      surf.dem = mdata.land_dem;
      surf.dem(isnan(surf.dem)) = mdata.sea_dem(isnan(surf.dem));
      surf.dem = (mdata.Elevation - surf.dem) / (c/2);
      surf.dem = interp1(mdata.Time,1:length(mdata.Time),surf.dem + surf.init.dem_offset);
      surf.dem = interp_finite(surf.dem,length(mdata.Time)/2);
    end
    
    %% Load LIDAR data if exists
    if isfield(orig_surf,'init') && isfield(orig_surf.init,'lidar_source')
      % Interpolate LIDAR onto RADAR time
      lidar_interp_gaps_dist = [150 75];
      ops_layer = [];
      ops_layer{1}.gps_time = layers(lay_idx).gps_time;
      ops_layer{1}.type = layers(lay_idx).type;
      ops_layer{1}.quality = layers(lay_idx).quality;
      ops_layer{1}.twtt = layers(lay_idx).twtt;
      ops_layer{1}.type(isnan(ops_layer{1}.type)) = 2;
      ops_layer{1}.quality(isnan(ops_layer{1}.quality)) = 1;
      lay = opsInterpLayersToMasterGPSTime(mdata,ops_layer,lidar_interp_gaps_dist);
      atm_layer = lay.layerData{1}.value{2}.data;
      atm_layer = interp1(mdata.Time,1:length(mdata.Time),atm_layer);
      surf.dem = merge_vectors(atm_layer, surf.dem);
    end
  
    surf.max_diff = orig_surf.max_diff/dt;
    
    %% Track the surface
    if surf.manual
      [new_surface,pnt] = tracker_snake_simple(mdata.Data,surf);
      fprintf('  Press F1 for help\n');
      layer = tracker_snake_manual_gui(lp(mdata.Data),pnt);

    elseif strcmpi(surf.method,'threshold')
      new_surface = tracker_threshold(mdata.Data,surf);
    elseif strcmpi(surf.method,'max')
      new_surface = tracker_max(mdata.Data,surf);
    elseif strcmpi(surf.method,'snake')
      new_surface = tracker_snake_simple(mdata.Data,surf);
    else
      error('Not a supported surface tracking method.');
    end

    %% Run median filtering on tracked surface
    if isfield(surf,'medfilt') && ~isempty(surf.medfilt)
      % OLD METHOD: new_surface = medfilt1(new_surface,surf.medfilt);
      for rline=1:length(new_surface)
        rlines = rline + (-surf.medfilt:surf.medfilt);
        rlines = rlines(rlines>=1 & rlines<=length(new_surface));
        if abs(new_surface(rline) - nanmedian(new_surface(rlines))) > surf.medfilt_threshold
          new_surface(rline) = nanmedian(new_surface(rlines));
        end
      end
    end
    
    % Convert from range bins to two way travel time
    Surface = interp1(1:length(mdata.Time), mdata.Time, new_surface);

    if debug_level > 0
      % Debug result
      figure(1); clf;
      imagesc([],mdata.Time,lp(mdata.Data));
      colormap(1-gray(256));
      hold on;
      plot(Surface,'m');
      if isfield(orig_surf,'init') && any(strcmpi(orig_surf.init.method,{'dem','lidar'}))
        plot(interp1(1:length(mdata.Time),mdata.Time,surf.dem),'g')
        plot(interp1(1:length(mdata.Time),mdata.Time,surf.dem-surf.max_diff),'r')
        plot(interp1(1:length(mdata.Time),mdata.Time,surf.dem+surf.max_diff),'b')
      end
      hold off;
      ylim([min(Surface)-debug_time_guard max(Surface)+debug_time_guard])
      keyboard
    end
    
    %% Save result
    if any(strcmpi(save_sources,'echogram'))
      save(data_fn,'-append','Surface');
    end
    
    if any(strcmpi(save_sources,'layerdata'))
      layer_fn = fullfile(ct_filename_out(param,layerdata_source,''), ...
        sprintf('Data_%s_%03d.mat', param.day_seg, frm));
      if exist(layer_fn,'file')
        % Load the layerData file
        lay = load(layer_fn);
        % Update the surface auto picks
        lay.layerData{1}.value{2}.data = interp1(mdata.GPS_time,Surface,lay.GPS_time);
        lay.layerData{1}.value{2}.data = interp_finite(lay.layerData{1}.value{2}.data);
        % Append the new results back to the layerData file
        save(layer_fn,'-append','-struct','lay','layerData');
      end
    end
    
    if any(strcmpi(save_sources,'ops'))
      % OPS query to get the point path ID's
      ops_param = struct('properties',[]);
      ops_param.properties.location = param.post.ops.location;
      ops_param.properties.season = param.season_name;
      ops_param.properties.start_gps_time = ops_seg_data.properties.start_gps_time(frm);
      ops_param.properties.stop_gps_time = ops_seg_data.properties.stop_gps_time(frm);
      
      sys = ct_output_dir(param.radar_name);
      [status,data] = opsGetPath(sys,ops_param);
      
      % Write the new layer information to these point path ID's
      ops_param = struct('properties',[]);
      ops_param.properties.point_path_id = data.properties.id;
      ops_param.properties.twtt = interp_finite(interp1(mdata.GPS_time,Surface,data.properties.gps_time));
      ops_param.properties.type = 2*ones(size(ops_param.properties.twtt));
      ops_param.properties.quality = 1*ones(size(ops_param.properties.twtt));
      ops_param.properties.lyr_name = save_ops_layer;
      
      opsCreateLayerPoints(sys,ops_param);
    end
  end
  
end






