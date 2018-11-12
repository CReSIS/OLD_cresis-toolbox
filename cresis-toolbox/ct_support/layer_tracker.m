function layer_tracker(param,param_override)
% layer_tracker(param,param_override)
%
% Tracks a layer using one of the layer tracking methods.
%
% track: Parameters used to control the tracker. Only general parameters
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
%     [track.feedthru.time,track.feedthru.power_dB] = ginput; % Press enter to end the input
%   .time: N length vector of two way travel times
%   .power_dB: N length vector of power in dB
% .method: string containing the method to use for tracking as long as
%   manual mode is not enabled. Options:
%    'threshold': runs tracker_threshold (generally the best for surface altimetry)
%    'snake': runs tracker_snake_simple
%    'max': tracker_max
% .manual: Optional, default is false. If true, the manual version of
%   tracker_snake_simple is run with tracker_snake_manual_gui.
% .max_diff: used by tracker routines (optional, default is inf)
% .min_bin: used by tracker routines
% .max_bin: used by tracker routines (optional, default is not used)
% .init.method: used by tracker routines (optional, default is not used)
% .init.dem_layer: dem layer struct to be used with opsLoadLayers
%   (usually lidar surface)
%
% For more details about tracker specific parameters:
% tracker_snake_simple, tracker_snake_manual_gui, tracker_threshold,
% tracker_max
%
% Example:
%   See run_layer_tracker.m to run.
%
% Author: John Paden


%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks
% =====================================================================

physical_constants;

% debug_level: set to 1 to enable a debug plot
if ~isfield(param.layer_tracker,'debug_level') || isempty(param.layer_tracker.debug_level)
  param.layer_tracker.debug_level = 0;
end
debug_level = param.layer_tracker.debug_level;

% debug_time_guard: Vertical band around surface in debug plot
if ~isfield(param.layer_tracker,'debug_time_guard') || isempty(param.layer_tracker.debug_time_guard)
  param.layer_tracker.debug_time_guard = 50e-9;
end
debug_time_guard = param.layer_tracker.debug_time_guard;

% echogram_img: To choose an image besides the base (0) image
if ~isfield(param.layer_tracker,'echogram_img') || isempty(param.layer_tracker.echogram_img)
  param.layer_tracker.echogram_img = 0;
end
echogram_img = param.layer_tracker.echogram_img;

% echogram_source: location of echogram data used for tracking
if ~isfield(param.layer_tracker,'echogram_source') || isempty(param.layer_tracker.echogram_source)
  error('An echogram_source must be specified.');
end
echogram_source = param.layer_tracker.echogram_source;

layer_params = param_override.layer_tracker.layer_params;

track = merge_structs(param.qlook.surf,param_override.layer_tracker.track);

if ~track.en
  return;
end

if ~isfield(track,'detrend') || isempty(track.detrend)
  track.detrend = 0;
end
if ischar(track.detrend)
  % Load detrend file generated by detrend
  % e.g. ct_tmp/echogram_stats/snow/2011_Greenland_P3/stats_20110329_02.mat
  detrend = load(track.detrend,'dt','bins','min_means');
  detrend.time = detrend.dt*detrend.bins;
end

if ~isfield(track,'init') || isempty(track.init)
  track.init = [];
end
if ~isfield(track.init,'method') || isempty(track.init.method)
  track.init.method = '';
end
if ~any(strcmpi(track.init.method,{'snake','medfilt','dem','lidar'}))
  error('Unsupported surface init method %s. Disable by setting to an empty string (default setting).', track.init.method);
end
if ~isfield(track.init,'max_diff') || isempty(track.init.max_diff)
  track.init.max_diff = inf;
end
if ~isfield(track.init,'threshold_dB') || isempty(track.init.threshold_dB)
  track.init.threshold_dB = 15;
end
if ~isfield(track.init,'threshold_noise_rng') || isempty(track.init.threshold_noise_rng)
  track.init.threshold_noise_rng = [0 -inf inf];
end

if ~isfield(track,'filter') || isempty(track.filter)
  track.filter = [1 1];
end
if length(track.filter) == 1
  warning('Deprecated surf.filter format. Should specify 2 element vector that specifies the multilooks in [cross-track along-track].');
  track.filter = [1 track.filter(1)];
end
if any(mod(track.filter,2) == 0)
  error('Surface filter lengths must be odd. layer_tracker.track.filter = [%d %d].', layer_tracker.track.filter);
end

if ~isfield(track,'filter_trim') || isempty(track.filter_trim)
  track.filter_trim = [0 0];
end

if ~isfield(track,'fixed_value') || isempty(track.fixed_value)
  track.fixed_value = 0;
end

if ~isfield(track,'min_bin') || isempty(track.min_bin)
  track.min_bin = 0;
end

if ~isfield(track,'max_bin') || isempty(track.max_bin)
  track.max_bin = inf;
end

if ~isfield(track,'medfilt') || isempty(track.medfilt)
  track.medfilt = 0;
end
if ~isfield(track,'medfilt_threshold') || isempty(track.medfilt_threshold)
  track.medfilt_threshold = inf;
end

if ~isfield(track,'prefilter_trim') || isempty(track.prefilter_trim)
  track.prefilter_trim = [0 0];
end

if ~isfield(track,'max_rng') || isempty(track.max_rng)
  track.max_rng = [0 0];
end

if ~isfield(track,'sidelobe_rows') || isempty(track.sidelobe_rows) || ~isfield(track,'sidelobe_dB') || isempty(track.sidelobe_dB)
  track.sidelobe_rows = [];
  track.sidelobe_dB = [];
end


%% Load in ocean mask, land DEM, and sea surface DEM
global load_surface_land_dems_finished;
global load_surface_land_dems_day_seg;
global ocean_shp_all;
global ocean_shp_bb;
global ocean_shp_day_seg;
global ocean_shp_bb_day_seg;
global land_surface;
global sea_surface;
load_surface_land_dems = false;
if isfield(track,'init') && strcmpi(track.init.method,'dem')
  if isempty(load_surface_land_dems_finished) ...
      || ~load_surface_land_dems_finished
    load_surface_land_dems = true;
    load_surface_land_dems_day_seg = '';
  end
end

if load_surface_land_dems
  % Load ocean mask shape file (-180 to +180 lon)
  ocean_mask_fn = ct_filename_gis([],fullfile('world','land_mask','Land_Mask_IDL_jharbeck','GSHHS_f_L1.shp'));
  warning off;
  %ocean_shp_all = shaperead(ocean_mask_fn, 'BoundingBox', [min_lon min_lat; max_lon max_lat]);
  ocean_shp_all = shaperead(ocean_mask_fn);
  warning on;
  % All bounding boxes of every shape
  ocean_shp_bb = [ocean_shp_all(:).BoundingBox];
  
  % Load land DEM
  if strcmpi(param.post.ops.location,'arctic')
    if 1
      land_surface.fn = ct_filename_gis([],'greenland/DEM/GIMP/gimpdem_90m.tif');
      land_surface.bad_value = 32767;
    else
      % PADEN: Load Arctic DEM corresponding to this segment
    end
  elseif strcmpi(param.post.ops.location,'antarctic')
    land_surface.fn = ct_filename_gis([],'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_surface.tif');
    land_surface.bad_value = 32767;
  end
  [land_surface.dem_all, land_surface.R, tmp] = geotiffread(land_surface.fn);
  land_surface.x_all = land_surface.R(3,1) + land_surface.R(2,1)*(0:size(land_surface.dem_all,2)-1);
  land_surface.y_all = land_surface.R(3,2) + land_surface.R(1,2)*(0:size(land_surface.dem_all,1)-1);
  land_surface.proj = geotiffinfo(land_surface.fn);
  
  load_surface_land_dems_finished = true;
end

if ~strcmpi(param.day_seg,load_surface_land_dems_day_seg) ...
    && isfield(track,'init') && strcmpi(track.init.method,'dem')
  % Load records file
  records_fn = ct_filename_support(param,'','records');
  records = load(records_fn);
  min_lat = min(records.lat);
  max_lat = max(records.lat);
  % Handle longitude in a special way because it wraps around.
  mean_lon = angle(mean(exp(1i*records.lon/180*pi)))*180/pi;
  max_lon = mean_lon + max(angle(exp(1i*(records.lon-mean_lon)/180*pi)))*180/pi;
  min_lon = mean_lon + min(angle(exp(1i*(records.lon-mean_lon)/180*pi)))*180/pi;
  
  [records.x,records.y] = projfwd(land_surface.proj,records.lat,records.lon);
  min_x = min(records.x);
  max_x = max(records.x);
  min_y = min(records.y);
  max_y = max(records.y);
  
  % Get all ocean shapes within the data segment bounding box. Shape is
  % not included if any of the following holds:
  %  - Bottom of the shape is above the top of the segment, >max_lat
  %  - Top of the shape is below the bottom of the segment, <min_lat
  %  - Left side of the shape is to the right of the segment, >max_lon
  %  - Right side of the shape is to the left of the segment, <min_lon
  % Handle longitude in a special way because it wraps around.
  if isempty(ocean_shp_bb)
    ocean_shp_day_seg = [];
  else
    rel_min_lon = mean_lon + angle(exp(1i*(ocean_shp_bb(1,1:2:end) - mean_lon)/180*pi))*180/pi;
    rel_max_lon = mean_lon + angle(exp(1i*(ocean_shp_bb(2,1:2:end) - mean_lon)/180*pi))*180/pi;
    bb_good_mask = ~(ocean_shp_bb(1,2:2:end)>max_lat | ocean_shp_bb(2,2:2:end)<min_lat ...
      | rel_min_lon>max_lon | rel_max_lon<min_lon);
    ocean_shp_day_seg = ocean_shp_all(bb_good_mask);
    % All bounding boxes of every shape
    ocean_shp_bb_day_seg = [ocean_shp_day_seg(:).BoundingBox];
  end
  
  if 0
    % Debug code to check bounding box code
    figure(1); clf;
    for idx=1:length(ocean_shp_day_seg)
      if length(ocean_shp_day_seg(idx).X) > 2000
        plot(ocean_shp_day_seg(idx).X(1:5:end),ocean_shp_day_seg(idx).Y(1:5:end))
      else
        plot(ocean_shp_day_seg(idx).X,ocean_shp_day_seg(idx).Y)
      end
      hold on;
    end
    plot(records.lon(1:100:end),records.lat(1:100:end),'k.')
  end
  
  % Load sea level data (0 to 360 lon data)
  if 0
    % EGM-96
    sea_surface.fn = ct_filename_gis([],'world\egm96_geoid\WW15MGH.DAC');
    points = [];
    [sea_surface.lat,sea_surface.lon,sea_surface.elev] = egm96_loader(sea_surface.fn);
    points.lon = [points.lon 360];
    sea_surface.elev = [sea_surface.elev sea_surface.elev(:,1)];
    [sea_surface.lon,sea_surface.lat] = meshgrid(sea_surface.lon,sea_surface.lat);
  else
    % Load DTU mean sea level
    sea_surface.fn = ct_filename_gis([],fullfile('world','dtu_meansealevel','DTU10MSS_1min.nc'));
    sea_surface.lat = ncread(sea_surface.fn,'lat');
    sea_surface.lon = ncread(sea_surface.fn,'lon');
    dlat = sea_surface.lat(2)-sea_surface.lat(1);
    lat_idxs = find(sea_surface.lat >= min_lat-2*dlat & sea_surface.lat <= max_lat+2*dlat);
    dlon = sea_surface.lon(2)-sea_surface.lon(1);
    rel_lon = mean_lon + angle(exp(1i*(sea_surface.lon - mean_lon)/180*pi))*180/pi;
    lon_idxs = find(rel_lon >= min_lon-2*dlon & rel_lon <= max_lon+2*dlon);
    break_idx = find(diff(lon_idxs)~=1);
    sea_surface.lat = sea_surface.lat(lat_idxs);
    sea_surface.lon = rel_lon(lon_idxs);
    % Transpose elev because "x" axis (which is longitude) must be on the
    % column dimension for interp2.
    % Convert to single because interp2 requires single or double type
    % and single is smaller yet has enough precision.
    if isempty(break_idx)
      sea_surface.elev = single(ncread(sea_surface.fn,'mss', ...
        [lon_idxs(1) lat_idxs(1)],[length(lon_idxs) length(lat_idxs)]).');
    else
      sea_surface.elev = single(ncread(sea_surface.fn,'mss', ...
        [lon_idxs(break_idx+1) lat_idxs(1)],[length(lon_idxs)-break_idx length(lat_idxs)]).');
      sea_surface.elev = [sea_surface.elev, single(ncread(sea_surface.fn,'mss', ...
        [1 lat_idxs(1)],[break_idx length(lat_idxs)]).')];
      sea_surface.lon = sea_surface.lon([break_idx+1:end,1:break_idx]);
    end
    [sea_surface.lon,unique_idxs] = unique(sea_surface.lon);
    sea_surface.elev = sea_surface.elev(:,unique_idxs);
  end
  
  % Load land DEM
  if strcmpi(param.post.ops.location,'arctic')
    dx = abs(land_surface.x_all(2)-land_surface.x_all(1));
    x_idxs = find(land_surface.x_all >= min_x-2*dx & land_surface.x_all <= max_x+2*dx);
    dy = abs(land_surface.y_all(2)-land_surface.y_all(1));
    y_idxs = find(land_surface.y_all >= min_y-2*dy & land_surface.y_all <= max_y+2*dy);
    land_surface.x = land_surface.x_all(x_idxs);
    land_surface.y = land_surface.y_all(y_idxs);
    land_surface.dem = single(land_surface.dem_all(y_idxs,x_idxs));
    
    if 0
      % PADEN: Load Arctic DEM corresponding to this segment
    end
  elseif strcmpi(param.post.ops.location,'antarctic')
    land_surface.fn = ct_filename_gis([],'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_surface.tif');
    [land_surface.dem, land_surface.R, tmp] = geotiffread(land_surface.fn);
    land_surface.proj = geotiffinfo(land_surface.fn);
  end
  
  load_surface_land_dems_day_seg = param.day_seg;
  
  if 0
    % Debug Plot
    figure(1); clf;
    imagesc(land_surface.x,land_surface.y,land_surface.dem)
    set(gca,'YDir','normal');
  end
end

%% Determine valid frames to process
load(ct_filename_support(param,'','frames'));
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

%% Load reference surface
if isfield(track,'init') && isfield(track.init,'dem_layer') ...
    && ~isempty(track.init.dem_layer)
  layer_params = track.init.dem_layer;
  
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

%% Track
orig_track = track;
for frm_idx = 1:length(param.cmd.frms)
  %% Track: Load echogram data
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
  [~,data_fn_name] = fileparts(data_fn);
  fprintf('%d of %d %s (%s)\n', frm_idx, length(param.cmd.frms), data_fn, datestr(now,'HH:MM:SS'));
  
  if ~exist(data_fn,'file')
    warning('  Missing file\n');
    continue;
  end
  track = orig_track;
  
  if strcmpi(track.method,'') && debug_level == 0
    mdata = load_L1B(data_fn,'GPS_time','Latitude','Longitude','Elevation','Time');
  else
    mdata = load_L1B(data_fn);
  end
  data = lp(mdata.Data);
  Nx = size(mdata.Data,2);
  
  %% Track: Interpolate GIMP and Geoid
  if strcmpi(track.init.method,'dem')
    
    mdata.sea_dem = interp2(sea_surface.lon,sea_surface.lat,sea_surface.elev,mdata.Longitude,mdata.Latitude);
    [mdata.x,mdata.y] = projfwd(land_surface.proj,mdata.Latitude,mdata.Longitude);
    if length(land_surface.x) > 2 && length(land_surface.y) > 2
      mdata.land_dem = interp2(land_surface.x,land_surface.y,land_surface.dem,mdata.x,mdata.y);
      mdata.land_dem(mdata.land_dem==land_surface.bad_value) = NaN;
    else
      mdata.land_dem = nan(size(mdata.GPS_time));
    end
    
    min_lat = min(mdata.Latitude);
    max_lat = max(mdata.Latitude);
    % Handle longitude in a special way because it wraps around.
    mean_lon = angle(mean(exp(1i*mdata.Longitude/180*pi)))*180/pi;
    max_lon = mean_lon + max(angle(exp(1i*(mdata.Longitude-mean_lon)/180*pi)))*180/pi;
    min_lon = mean_lon + min(angle(exp(1i*(mdata.Longitude-mean_lon)/180*pi)))*180/pi;
    min_x = min(mdata.x);
    max_x = max(mdata.x);
    min_y = min(mdata.y);
    max_y = max(mdata.y);
    
    % Restrict ocean mask features to our dataset (i.e. mask all features
    % whose bounding boxes fall outside our limits.
    if isempty(ocean_shp_day_seg)
      ocean_shp = [];
    else
      rel_min_lon = mean_lon + angle(exp(1i*(ocean_shp_bb_day_seg(1,1:2:end) - mean_lon)/180*pi))*180/pi;
      rel_max_lon = mean_lon + angle(exp(1i*(ocean_shp_bb_day_seg(2,1:2:end) - mean_lon)/180*pi))*180/pi;
      bb_good_mask = ~(ocean_shp_bb_day_seg(1,2:2:end)>max_lat | ocean_shp_bb_day_seg(2,2:2:end)<min_lat ...
        | rel_min_lon>max_lon | rel_max_lon<min_lon);
      ocean_shp = ocean_shp_day_seg(bb_good_mask);
    end
    
    % Create polygons, poly_x/poly_y, with all ocean shape features which
    % lie in the bounding box.
    % Further restrict the polygons by checking for bounding box overlap in
    % projected coordinates
    
    poly_x = cell(0);
    poly_y = cell(0);
    for shp_idx = 1:length(ocean_shp)
      % convert polygon to projected coordinates
      [x,y] = projfwd(land_surface.proj,ocean_shp(shp_idx).Y,ocean_shp(shp_idx).X);
      % if polygon is within projected bounding box
      if min(x) < max_x && max(x) > min_x ...
          && min(y) < max_y && max(y)>min_y
        % add polygon
        poly_x{end+1} = [x,nan];
        poly_y{end+1} = [y,nan];
      end
    end
    
    % Create ocean mask to determine which points lie in the ocean
    ocean_mask = true(size(mdata.Latitude));
    for poly_idx = 1:length(poly_x)
      
      % Mask showing which DEM points are in polygon (on land)
      land_mask_tmp = inpolygon(mdata.x,mdata.y,[poly_x{poly_idx}(1:100:end)],[poly_y{poly_idx}(1:100:end)]);
      ocean_mask(land_mask_tmp) = false;
    end
    
    % Merge land surface and sea surface DEMs
    track.dem = mdata.land_dem;
    track.dem(ocean_mask) = mdata.sea_dem(ocean_mask);
    track.dem = (mdata.Elevation - track.dem) / (c/2);
    track.dem = interp1(mdata.Time,1:length(mdata.Time),track.dem + track.init.dem_offset,'linear','extrap');
    track.dem = interp_finite(track.dem,1);
    track.dem(track.dem<1) = 1;
    track.dem(track.dem>length(mdata.Time)) = length(mdata.Time);
  end

  %% Track: Prefilter trim
  % Also set leading/following zeros or ~isfinite to NaN
  for rline = 1:Nx
    start_bin = find(isfinite(data(:,rline)),1);
    if ~isempty(start_bin)
      stop_bin = min(size(data,1), start_bin+track.prefilter_trim(1)-1);
      data(1:stop_bin,rline) = NaN;
    end
    stop_bin = find(isfinite(data(:,rline)),1,'last');
    if ~isempty(stop_bin)
      start_bin = max(1, stop_bin-track.prefilter_trim(2)+1);
      data(start_bin:end,rline) = NaN;
    end
  end
  
  %% Track: Feed through removal
  if isfield(track,'feedthru')
    
    % Interpolate feed through power levels on to data time axis
    feedthru_threshold = interp1(track.feedthru.time,track.feedthru.power_dB,mdata.Time);
    feedthru_threshold = interp_finite(feedthru_threshold,-inf);
    
    % Set all data to zero that does not exceed the feed through
    % threshold power
    for rline=1:Nx
      data(data(:,rline) <= feedthru_threshold,rline) = NaN;
    end
  end
  
  %% Track: Detrend
  if ischar(track.detrend)
    detrend_curve = interp_finite(interp1(detrend.time,interp_finite(detrend.min_means),mdata.Time),NaN);
    if 0
      % Debug
      rline = 200;
      figure(1); clf;
      plot(data(:,rline))
      hold on
      mean_power = nanmean(data,2);
      plot(mean_power)
      plot(detrend_curve);
      keyboard
    end
    
  elseif track.detrend > 0
    poly_x = (-size(data,1)/2+(1:size(data,1))).';
    mean_power = nanmean(data,2);
    good_mask = isfinite(mean_power);
    p = polyfit(poly_x(good_mask),mean_power(good_mask),track.detrend);
    detrend_curve = polyval(p,poly_x);
    detrend_curve(~good_mask) = NaN;
    detrend_curve = interp_finite(detrend_curve,0);
    if 0
      % Debug
      rline = 200;
      figure(1); clf;
      plot(data(:,rline))
      hold on
      plot(mean_power)
      plot(detrend_curve);
      keyboard
    end
    
  end
  data = bsxfun(@minus,data,detrend_curve);
  
  %% Track: Sidelobe
  if ~isempty(track.sidelobe_rows)
    mask = sidelobe_mask_mex(single(data),int32(track.sidelobe_rows),single(track.sidelobe_dB));
    data(mask) = NaN;
  end
  
  %% Track: Filter
  if track.filter(1) ~= 1
    % Multilooking in cross-track/fast-time
    data = lp(nan_fir_dec(10.^(data/10).',ones(1,track.filter(1))/track.filter(1),1,[],[],[],[],2.0).');
  end
  if track.filter(2) ~= 1
    % Multilooking in along-track
    data = lp(nan_fir_dec(10.^(data/10),ones(1,track.filter(2))/track.filter(2),1,[],[],[],[],2.0));
  end
  
  %% Track: Post-filter trim
  for rline = 1:Nx
    start_bin = find(isfinite(data(:,rline)),1);
    if ~isempty(start_bin)
      stop_bin = min(size(data,1), start_bin+track.filter_trim(1)-1);
      data(start_bin:stop_bin,rline) = NaN;
    end
    stop_bin = find(isfinite(data(:,rline)),1,'last');
    if ~isempty(stop_bin)
      start_bin = max(1, stop_bin-track.filter_trim(2)+1);
      data(start_bin:stop_bin,rline) = NaN;
    end
  end
  
  %% Track: min_bin/max_bin + time to bins conversions
  % Convert from two way travel time to bins
  track.min_bin = find(mdata.Time >= orig_track.min_bin, 1);
  track.max_bin = find(mdata.Time <= orig_track.max_bin, 1, 'last');
  dt = mdata.Time(2) - mdata.Time(1);
  track.init.max_diff = orig_track.init.max_diff/dt;
  track.max_rng = round(orig_track.max_rng(1)/dt) : round(orig_track.max_rng(end)/dt);
  track.threshold_noise_rng = round(orig_track.threshold_noise_rng/dt);
  data = data(track.min_bin:track.max_bin,:);
  
  %% Track: Create Initial Surface
  if strcmpi(track.init.method,'dem')
    % Correct for min_bin removal
    track.dem = track.dem - track.min_bin + 1;
  else
    if strcmp(track.init.method,'snake')
      track.dem = tracker_snake_simple(data,track.init);
    else
      [~,track.dem] = max(data,[],1);
      if strcmp(track.init.method,'medfilt')
        track.dem = medfilt1(track.dem,track.init.medfilt);
      end
    end
  end
  
  %% Track: Merge DEM and reference layer
  if isfield(track,'init') && isfield(track.init,'dem_layer') ...
      && ~isempty(track.init.dem_layer)
    % Interpolate reference layer onto radar GPS time
    ref_interp_gaps_dist = [150 75];
    ops_layer = [];
    ops_layer{1}.gps_time = layers(lay_idx).gps_time;
    ops_layer{1}.type = layers(lay_idx).type;
    ops_layer{1}.quality = layers(lay_idx).quality;
    ops_layer{1}.twtt = layers(lay_idx).twtt;
    ops_layer{1}.type(isnan(ops_layer{1}.type)) = 2;
    ops_layer{1}.quality(isnan(ops_layer{1}.quality)) = 1;
    lay = opsInterpLayersToMasterGPSTime(mdata,ops_layer,ref_interp_gaps_dist);
    dem_layer = lay.layerData{1}.value{2}.data;
    dem_layer = interp1(mdata.Time,1:length(mdata.Time),dem_layer);
    track.dem = merge_vectors(dem_layer, track.dem);
  end
  
  %% Track: Tracking
  if strcmpi(track.method,'threshold')
    [new_layer,new_quality] = tracker_threshold(data,track);
  elseif strcmpi(track.method,'viterbi')
    new_layer = tracker_viterbi(data,track);
    new_quality = ones(1,Nx);
  elseif strcmpi(track.method,'max')
    new_layer = tracker_max(data,track);
    new_quality = ones(1,Nx);
  elseif strcmpi(track.method,'snake')
    new_layer = tracker_snake_simple(data,track);
    new_quality = ones(1,Nx);
  elseif strcmpi(track.method,'fixed')
    new_layer = ones(size(mdata.GPS_time)) * track.fixed_value;
    new_quality = ones(1,Nx);
  elseif isempty(track.method)
    new_layer = track.dem;
    new_quality = ones(1,Nx);
  else
    error('Not a supported layer tracking method.');
  end

  %% Track: max_diff
  new_layer(abs(new_layer  - track.dem) > track.init.max_diff) = NaN;
  if any(strcmpi(track.init.method,{'dem','lidar'}))
    new_layer = merge_vectors(new_layer, track.dem);
  else
    new_layer = interp_finite(new_layer,NaN);
  end
  
  %% Track: Median filtering
  if isfield(track,'medfilt') && ~isempty(track.medfilt)
    % OLD METHOD: new_layer = medfilt1(new_layer,track.medfilt);
    for rline=1:Nx
      rlines = rline + (-track.medfilt:track.medfilt);
      rlines = rlines(rlines>=1 & rlines<=length(new_layer));
      if abs(new_layer(rline) - nanmedian(new_layer(rlines))) > track.medfilt_threshold
        new_layer(rline) = nanmedian(new_layer(rlines));
      end
    end
  end
  
  %% Track: Max search
  if (length(track.max_rng) > 1 || track.max_rng ~= 0)
    % Find the next peak after the threshold
    for rline = 1:Nx
      search_bins = round(new_layer(rline)) + track.max_rng;
      search_bins = search_bins(find(search_bins >= 1 & search_bins <= size(data,1)));
      [~,offset] = max(data(search_bins,rline));
      if ~isempty(offset)
        new_layer(rline) = search_bins(offset);
      end
    end
  end
  
  %% Track: Convert bins to twtt
  Surface = interp1(1:length(mdata.Time), mdata.Time, new_layer + track.min_bin);
  
  %% Track: Debug plot
  if debug_level > 0
    figure(1); clf;
    imagesc([],mdata.Time,lp(mdata.Data));
    colormap(1-gray(256));
    hold on;
    plot(find(new_quality==1),Surface(new_quality==1),'g.');
    plot(find(new_quality==3),Surface(new_quality==3),'r.');
    if any(strcmpi(track.init.method,{'dem','lidar'}))
      plot(interp1(1:length(mdata.Time),mdata.Time,track.dem),'m--')
      plot(interp1(1:length(mdata.Time),mdata.Time,track.dem-track.max_diff),'r--')
      plot(interp1(1:length(mdata.Time),mdata.Time,track.dem+track.max_diff),'b--')
    end
    hold off;
    ylim([max(mdata.Time(1),min(Surface)-debug_time_guard) min(mdata.Time(end),max(Surface)+debug_time_guard)])
    title(sprintf('%s',regexprep(data_fn_name,'_','\\_')));
    keyboard
  end
  
  %% Track: Save
  for layer_idx = 1:length(layer_params)
    layer_param = layer_params(layer_idx);
    
    if strcmpi(layer_param.source,'echogram')
      if isempty(layer_param.echogram_source)
        layer_param.echogram_source = echogram_source;
      end
      
      data_fn = fullfile(ct_filename_out(param,layer_param.echogram_source,''), ...
        sprintf('Data_%s_%03d.mat', param.day_seg, frm));
      fprintf('  Saving %s (%s)\n', data_fn, datestr(now));
      save(data_fn,'-append','Surface');
    end
    
    if strcmpi(layer_param.source,'layerdata')
      layer_fn = fullfile(ct_filename_out(param,layer_param.layerdata_source,''), ...
        sprintf('Data_%s_%03d.mat', param.day_seg, frm));
      if ~exist(layer_fn,'file')
        fprintf('  Create  %s (%s)\n', layer_fn, datestr(now));
        
        lay.GPS_time = mdata.GPS_time;
        lay.Latitude = mdata.Latitude;
        lay.Longitude = mdata.Longitude;
        lay.Elevation = mdata.Elevation;
        
        lay.layerData{1}.quality = interp1(mdata.GPS_time,quality,lay.GPS_time,'nearest');
        lay.layerData{1}.value{1}.data = nan(size(lay.GPS_time));
        lay.layerData{1}.value{2}.data = interp1(mdata.GPS_time,Surface,lay.GPS_time);
        lay.layerData{1}.value{2}.data = interp_finite(lay.layerData{1}.value{2}.data);
        
        lay.layerData{2}.quality = ones(size(lay.GPS_time));
        lay.layerData{2}.value{1}.data = nan(size(lay.GPS_time));
        lay.layerData{2}.value{2}.data = nan(size(lay.GPS_time));
        
        layer_fn_dir = fileparts(layer_fn);
        if ~exist(layer_fn_dir,'dir')
          mkdir(layer_fn_dir);
        end
        save(layer_fn,'-struct','lay');
        
      else
        % Load the layerData file
        lay = load(layer_fn);
        % Update the surface auto picks
        lay.layerData{1}.quality = interp1(mdata.GPS_time,new_quality,lay.GPS_time,'nearest');
        lay.layerData{1}.value{2}.data = interp1(mdata.GPS_time,Surface,lay.GPS_time);
        lay.layerData{1}.value{2}.data = interp_finite(lay.layerData{1}.value{2}.data,NaN);
        % Append the new results back to the layerData file
        fprintf('  Saving %s (%s)\n', layer_fn, datestr(now));
        save(layer_fn,'-append','-struct','lay','layerData');
      end
    end
    
    if strcmpi(layer_param.source,'ops')
      % Get all the frames for this segment
      if any(strcmpi({layer_params.source},'ops'))
        opsAuthenticate(param,false);
        sys = ct_output_dir(param.radar_name);
        ops_param = struct('properties',[]);
        ops_param.properties.season = param.season_name;
        ops_param.properties.segment = param.day_seg;
        [status,ops_seg_data] = opsGetSegmentInfo(sys,ops_param);
      end
      
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
      ops_param.properties.quality = interp1(mdata.GPS_time,new_quality,data.properties.gps_time,'nearest');
      ops_param.properties.lyr_name = layer_param.name;
      
      opsCreateLayerPoints(sys,ops_param);
    end
    
  end
  
end
