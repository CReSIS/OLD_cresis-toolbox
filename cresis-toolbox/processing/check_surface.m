function check_surface(param,param_override)
% check_surface(param,param_override)
%
% 1. Loads coincident LIDAR data if it exists
% 2. Loads DTU sea surface DEM and arctic/antarctica land DEM, combines
%    these two DEMS taking land DEM over sea surface DEM.
% 3. Combines LIDAR data and DEM data, taking LIDAR data over DEM data.
%    Uses elevation to interpolate where data are not available.
% 4. Estimates Tadc_adjust or t_ref error by comparing radar surface from
%    the specified layer source and the LIDAR/DEM combination.
%    This error should be subtracted from param.radar.wfs.Tadc_adjust for pulsed systems.
%    This error should be added to param.radar.wfs.t_ref for deramp systems.
% 5. Estimates GPS offset by comparing radar surface and LIDAR/DEM. This offset
%    should be added to param.records.gps.time_offset.
% 6. For deramp systems, uses the LIDAR/DEM data to determine the Nyquist
%    zone and sets the records.settings.nyquist_zone based on this.
%    The second decimal mask in frames.proc_mode is also set to one for
%    frames that will be outside max_nyquist_zone.
%
% See run_check_surface.m for how to run.
%
% cat /N/dcwan/projects/cresis/output/ct_tmp/check_surface/snow/2017_Greenland_P3/*.txt
%
%
% Author: John Paden

% =====================================================================
%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

% =====================================================================
%% Input checks
% =====================================================================

physical_constants();

if ~isfield(param,'check_surface') || isempty(param.check_surface)
  param.check_surface = [];
end

if ~isfield(param.check_surface,'use_lidar_data') || isempty(param.check_surface.use_lidar_data)
  param.check_surface.use_lidar_data = true;
end
use_lidar_data = param.check_surface.use_lidar_data;

if ~isfield(param.check_surface,'combine_elev_lidar_en') || isempty(param.check_surface.combine_elev_lidar_en)
  param.check_surface.combine_elev_lidar_en = true;
end
combine_elev_lidar_en = param.check_surface.combine_elev_lidar_en;

if ~isfield(param.check_surface,'combine_surface_land_dems') || isempty(param.check_surface.combine_surface_land_dems)
  param.check_surface.combine_surface_land_dems = true;
end
combine_surface_land_dems = param.check_surface.combine_surface_land_dems;

if ~isfield(param.check_surface,'max_twtt_diff') || isempty(param.check_surface.max_twtt_diff)
  param.check_surface.max_twtt_diff = 200e-9;
end
max_twtt_diff = param.check_surface.max_twtt_diff;

if ~isfield(param.check_surface,'debug_level') || isempty(param.check_surface.debug_level)
  param.check_surface.debug_level = 0;
end
debug_level = param.check_surface.debug_level;

if ~isfield(param.check_surface,'save_records_en') || isempty(param.check_surface.save_records_en)
  param.check_surface.save_records_en = false;
end
save_records_en = param.check_surface.save_records_en;

if ~isfield(param.check_surface,'refine_Tsys_en') || isempty(param.check_surface.refine_Tsys_en)
  param.check_surface.refine_Tsys_en = false;
end
refine_Tsys_en = param.check_surface.refine_Tsys_en;

if ~isfield(param.check_surface,'lidar_interp_gaps_dist') || isempty(param.check_surface.lidar_interp_gaps_dist)
  param.check_surface.lidar_interp_gaps_dist = [150 75];
end
lidar_interp_gaps_dist = param.check_surface.lidar_interp_gaps_dist;

if ~isfield(param.check_surface,'radar_twtt_ratio') || isempty(param.check_surface.radar_twtt_ratio)
  % Default one: this value will be multiplied with the radar twtt. Used to
  % test incorrect sampling frequency (e.g. deramp on receive with wrong
  % f0/f1/Tpd parameters).
  param.check_surface.radar_twtt_ratio = 1.0;
end
radar_twtt_ratio = param.check_surface.radar_twtt_ratio;

if ~isfield(param.check_surface,'radar_twtt_offset') || isempty(param.check_surface.radar_twtt_offset)
  % Default zero: this value will be added to the radar twtt. Used to test
  % offsets.
  param.check_surface.radar_twtt_offset = 0.0;
end
radar_twtt_offset = param.check_surface.radar_twtt_offset;

if ~isfield(param.check_surface,'radar_layer_params') || isempty(param.check_surface.radar_layer_params)
  param.check_surface.radar_layer_params.name = 'surface';
  param.check_surface.radar_layer_params.source = 'layerdata';
end

if ~isfield(param.check_surface,'ref_layer_params') || isempty(param.check_surface.ref_layer_params)
  param.check_surface.ref_layer_params.name = 'surface';
  param.check_surface.ref_layer_params.source = 'lidar';
  param.check_surface.ref_layer_params.lidar_source = 'atm';
end

layer_params = cat_structs(2,param.check_surface.ref_layer_params,param.check_surface.radar_layer_params);

% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);

% Load frames file
frames_fn = ct_filename_support(param,'','frames');
load(frames_fn);

% =========================================================================
%% Load in ocean mask, land DEM, and sea surface DEM
% =========================================================================
global load_surface_land_dems_finished;
global load_surface_land_dems_day_seg;
global ocean_shp_all;
global ocean_shp_bb;
global ocean_shp_day_seg;
global ocean_shp_bb_day_seg;
global land_surface;
global sea_surface;
load_surface_land_dems = false;
if isempty(load_surface_land_dems_finished) ...
    || ~load_surface_land_dems_finished
  load_surface_land_dems = true;
  load_surface_land_dems_day_seg = '';
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

if ~strcmpi(param.day_seg,load_surface_land_dems_day_seg)
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
    hold on;
    plot(records.x,records.y,'r','LineWidth',2)
    set(gca,'YDir','normal');
    xlabel('X (m)');
    ylabel('Y (m)');
    legend('Flightline','location','best');
  end
end

% =====================================================================
%% Load layer data
% =====================================================================

ref_idx = 1; % Make the radar the master "slow" time axis
radar_idx = 2; % Make the radar the master "slow" time axis

% Load radar surface (default layerdata) and reference surface (default ATM
% lidar)
layers = opsLoadLayers(param,layer_params);

% Ensure that layer gps times are monotonically increasing
for lay_idx = 1:length(layers)
  layers_fieldnames = fieldnames(layers(lay_idx));
  [~,unique_idxs] = unique(layers(lay_idx).gps_time);
  for field_idx = 1:length(layers_fieldnames)-1
    if ~isempty(layers(lay_idx).(layers_fieldnames{field_idx}))
      layers(lay_idx).(layers_fieldnames{field_idx}) = layers(lay_idx).(layers_fieldnames{field_idx})(unique_idxs);
    end
  end
end

% Throw out low quality radar data
layers(radar_idx).twtt(layers(radar_idx).quality==3) = NaN;

% Interpolate ref layer to radar gps time
master = [];
master.GPS_time = layers(radar_idx).gps_time;
master.Latitude = layers(radar_idx).lat;
master.Longitude = layers(radar_idx).lon;
master.Elevation = layers(radar_idx).elev;

ops_layer = [];
ops_layer{1}.gps_time = layers(ref_idx).gps_time;
ops_layer{1}.type = layers(ref_idx).type;
ops_layer{1}.quality = layers(ref_idx).quality;
ops_layer{1}.twtt = layers(ref_idx).twtt;
ops_layer{1}.type(isnan(ops_layer{1}.type)) = 2;
ops_layer{1}.quality(isnan(ops_layer{1}.quality)) = 1;
lay = opsInterpLayersToMasterGPSTime(master,ops_layer,lidar_interp_gaps_dist);
layers(ref_idx).twtt_ref = lay.layerData{1}.value{2}.data;

% Stretch and offset the twtt
layers(radar_idx).twtt = layers(radar_idx).twtt * radar_twtt_ratio;
layers(radar_idx).twtt = layers(radar_idx).twtt + radar_twtt_offset;

% =====================================================================
%% Merge sea-surface/land DEM and elevation information with layer info
% =====================================================================

mdata.Longitude = layers(radar_idx).lon;
mdata.Latitude = layers(radar_idx).lat;
mdata.Elevation = layers(radar_idx).elev;

mdata.sea_dem = interp2(sea_surface.lon,sea_surface.lat,sea_surface.elev,mdata.Longitude,mdata.Latitude);
[mdata.x,mdata.y] = projfwd(land_surface.proj,mdata.Latitude,mdata.Longitude);
if length(land_surface.x) > 2 && length(land_surface.y) > 2
  mdata.land_dem = interp2(land_surface.x,land_surface.y,land_surface.dem,mdata.x,mdata.y);
  mdata.land_dem(mdata.land_dem==land_surface.bad_value) = NaN;
else
  mdata.land_dem = nan(size(layers(radar_idx).gps_time));
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
surf.dem = mdata.land_dem;
surf.dem(ocean_mask) = mdata.sea_dem(ocean_mask);
surf.dem_twtt = (mdata.Elevation - surf.dem) / (c/2);

% Merge land/sea surface with reference (usually lidar) layer
surf.dem_twtt = merge_vectors(layers(ref_idx).twtt_ref, surf.dem_twtt);

mdata.land_dem_twtt = (mdata.Elevation - mdata.land_dem) / (c/2);
mdata.sea_dem_twtt = (mdata.Elevation - mdata.sea_dem) / (c/2);

% =====================================================================
%% Check surface: System time delay
% =====================================================================

lay_idx = 1;
twtt_error_all = surf.dem_twtt - layers(radar_idx).twtt;
twtt_error = twtt_error_all;
twtt_error(abs(twtt_error) > max_twtt_diff) = NaN;

mean_offset = nanmean(twtt_error);
origin = layers(radar_idx).gps_time(1);

h_fig = figure(1); clf(h_fig);
h_axes = axes('parent',h_fig);
plot(h_axes,layers(radar_idx).gps_time - origin, twtt_error_all*1e9)
hold(h_axes,'on');
plot(h_axes,layers(radar_idx).gps_time - origin, twtt_error*1e9,'g.')
xlabel(h_axes,sprintf('Relative GPS time (sec from %s)', datestr(epoch_to_datenum(origin))));
ylabel(h_axes,'TWTT error (ns)');
grid(h_axes,'on');
if max(twtt_error)>min(twtt_error)
  ylim([min(twtt_error) max(twtt_error)]*1e9);
end

fig_fn = ct_filename_ct_tmp(param,'','check_surface','twtt_error');
fprintf('Saving %s\n', fig_fn);
fig_fn_dir = fileparts(fig_fn);
if ~exist(fig_fn_dir,'dir')
  mkdir(fig_fn_dir);
end
saveas(h_fig,[fig_fn '.fig']);
saveas(h_fig,[fig_fn '.jpg']);

% =====================================================================
%% Check surface: GPS Offset
% =====================================================================

% Find the longest contiguous section of small twtt_error
% 1. Try LIDAR-only first
mask = ~isnan(twtt_error) & ~isnan(layers(ref_idx).twtt_ref);
mask_length = zeros(size(mask));
mask_length(1) = mask(1);
for idx=2:length(mask)
  mask_length(idx) = mask(idx)*mask_length(idx-1) + mask(idx);
end
[corr_len,corr_idx] = max(mask_length);
if corr_len==0
  recs = [];
else
  gpstime_coords = layers(radar_idx).gps_time(corr_idx+[-corr_len+1,0]);
  recs = find(layers(radar_idx).gps_time >= min(gpstime_coords) ...
    & layers(radar_idx).gps_time <= max(gpstime_coords));
end
dem_source = 'lidar';
% 2. Try LIDAR+Land next
if length(recs)<1000
  mask = ~isnan(twtt_error) & ~isnan(mdata.land_dem);
  mask_length = zeros(size(mask));
  mask_length(1) = mask(1);
  for idx=2:length(mask)
    mask_length(idx) = mask(idx)*mask_length(idx-1) + mask(idx);
  end
  [corr_len,corr_idx] = max(mask_length);
  if corr_len==0
    recs = [];
  else
    gpstime_coords = layers(radar_idx).gps_time(corr_idx+[-corr_len+1,0]);
    recs = find(layers(radar_idx).gps_time >= min(gpstime_coords) ...
      & layers(radar_idx).gps_time <= max(gpstime_coords));
  end
  dem_source = 'land';
end
% 3. Try LIDAR+Land+Sea next
if length(recs)<1000
  mask = ~isnan(twtt_error);
  mask_length = zeros(size(mask));
  mask_length(1) = mask(1);
  for idx=2:length(mask)
    mask_length(idx) = mask(idx)*mask_length(idx-1) + mask(idx);
  end
  [corr_len,corr_idx] = max(mask_length);
  if corr_len==0
    recs = [];
  else
    gpstime_coords = layers(radar_idx).gps_time(corr_idx+[-corr_len+1,0]);
    recs = find(layers(radar_idx).gps_time >= min(gpstime_coords) ...
      & layers(radar_idx).gps_time <= max(gpstime_coords));
  end
  dem_source = 'sea';
end

debug_gps_offset = 0;
debug_Tsys_offset = 0;
debug_Tsys_ratio = 1;

if any(~isnan(mdata.land_dem(recs)))
end
if any(~isnan(layers(ref_idx).twtt_ref(recs)))
end

clf(h_fig);
h_axes = axes('parent',h_fig);
origin = layers(radar_idx).gps_time(1);
h_plot = [];
frms = interp1([records.gps_time(frames.frame_idxs), records.gps_time(end)+diff(records.gps_time(end-1:end))], ...
  [1:length(frames.frame_idxs), length(frames.frame_idxs)+1], layers(radar_idx).gps_time);
h_plot(end+1) = plot(h_axes,frms, mdata.land_dem_twtt);
hold(h_axes,'on');
h_plot(end+1) = plot(h_axes,frms, mdata.sea_dem_twtt);
h_plot(end+1) = plot(h_axes,frms, layers(ref_idx).twtt_ref);
h_plot(end+1) = plot(h_axes,frms, surf.dem_twtt);
if ~isempty(recs)
  h_plot(end+1) = plot(h_axes,frms(recs), mdata.land_dem_twtt(recs), 'x');
  h_plot(end+1) = plot(h_axes,frms(recs), mdata.sea_dem_twtt(recs), 'o');
  h_plot(end+1) = plot(h_axes,frms(recs), layers(ref_idx).twtt_ref(recs), '<');
  h_plot(end+1) = plot(h_axes,frms(recs), surf.dem_twtt(recs), '.');
  for idx=1:4
    set(h_plot(idx+4),'Color',get(h_plot(idx),'Color'));
  end
end
frms = interp1([records.gps_time(frames.frame_idxs), records.gps_time(end)+diff(records.gps_time(end-1:end))] + debug_gps_offset, ...
  [1:length(frames.frame_idxs), length(frames.frame_idxs)+1], layers(radar_idx).gps_time);
h_plot(9) = plot(h_axes,frms, layers(radar_idx).twtt*debug_Tsys_ratio + debug_Tsys_offset,'k','LineWidth',2);
legend(h_axes,h_plot([1:4 9]),'Land','Sea','Ref','Combined','Radar','location','best');
grid(h_axes,'on');
xlabel(h_axes,'Frame');
ylabel(h_axes,'TWTT (sec)');
fig_fn = ct_filename_ct_tmp(param,'','check_surface','twtt_frm');
fprintf('Saving %s\n', fig_fn);
saveas(h_fig,[fig_fn '.fig']);
saveas(h_fig,[fig_fn '.jpg']);

clf(h_fig);
h_axes = axes('parent',h_fig);
origin = layers(radar_idx).gps_time(1);
h_plot = [];
h_plot(end+1) = plot(h_axes,layers(radar_idx).gps_time - origin, mdata.land_dem_twtt);
hold(h_axes,'on');
h_plot(end+1) = plot(h_axes,layers(radar_idx).gps_time - origin, mdata.sea_dem_twtt);
h_plot(end+1) = plot(h_axes,layers(radar_idx).gps_time - origin, layers(ref_idx).twtt_ref);
h_plot(end+1) = plot(h_axes,layers(radar_idx).gps_time - origin, surf.dem_twtt);
if ~isempty(recs)
  h_plot(end+1) = plot(h_axes,layers(radar_idx).gps_time(recs) - origin, mdata.land_dem_twtt(recs), 'x');
  h_plot(end+1) = plot(h_axes,layers(radar_idx).gps_time(recs) - origin, mdata.sea_dem_twtt(recs), 'o');
  h_plot(end+1) = plot(h_axes,layers(radar_idx).gps_time(recs) - origin, layers(ref_idx).twtt_ref(recs), '<');
  h_plot(end+1) = plot(h_axes,layers(radar_idx).gps_time(recs) - origin, surf.dem_twtt(recs), '.');
  for idx=1:4
    set(h_plot(idx+4),'Color',get(h_plot(idx),'Color'));
  end
end
h_plot(9) = plot(h_axes,layers(radar_idx).gps_time + debug_gps_offset - origin, layers(radar_idx).twtt*debug_Tsys_ratio + debug_Tsys_offset,'k','LineWidth',2);
legend(h_axes,h_plot([1:4 9]),'Land','Sea','Ref','Combined','Radar','location','best');
grid(h_axes,'on');
xlabel(h_axes,sprintf('Relative GPS time (sec from %s)', datestr(epoch_to_datenum(origin))));
ylabel(h_axes,'TWTT (sec)');
fig_fn = ct_filename_ct_tmp(param,'','check_surface','twtt');
fprintf('Saving %s\n', fig_fn);
saveas(h_fig,[fig_fn '.fig']);
saveas(h_fig,[fig_fn '.jpg']);

if 0
  % Correlate to determine GPS offsets
  figure(2);
  fprintf('Pick two points on figure 1 to constrain the cross correlation to a section of lidar elevation and radar elevation that are correct (but possibly offset in radar.wfs.Tsys and gps.time_offset).\n');
  fprintf('For each click hold mouse button still after click until cross hairs re-appear\n');
  [gpstime_coords,tmp] = ginput(2);
end

% Uniformly time sample the two signals
dt = median(diff(layers(radar_idx).gps_time(recs)));
if isempty(recs)
  ref_corr = NaN;
  lags = NaN;
  peak_idx = 1;
else
  t0 = layers(radar_idx).gps_time(recs(1));
  gps_time = t0:dt:layers(radar_idx).gps_time(recs(end));
  radar_layer = interp1(layers(radar_idx).gps_time, ...
    layers(radar_idx).twtt, gps_time);
  if length(layers(radar_idx).gps_time) < 2
    ref_layer = NaN(size(gps_time));
  else
    ref_layer = interp1(layers(radar_idx).gps_time, ...
      surf.dem_twtt, gps_time);
  end
  
  radar_layer = radar_layer - mean(radar_layer);
  ref_layer = ref_layer - mean(ref_layer);
  
  if 0
    [ref_corr,lags] = xcorr(radar_layer,ref_layer);
  else
    max_lag = round(10/dt);
    lags = -max_lag:max_lag;
    ref_corr = zeros(1,length(lags));
    for lag_idx = 1:length(lags)
      lag = lags(lag_idx);
      ref_corr(lag_idx) = sum(abs(radar_layer(1+max_lag:end-max_lag)-ref_layer(lag + [1+max_lag:end-max_lag])));
    end
    ref_corr = 1./ref_corr;
  end
  [peak_val,peak_idx] = max(ref_corr);
  if peak_idx == 1 || isnan(peak_val)
    lags(peak_idx) = NaN;
  end
    
end

clf(h_fig);
h_axes = axes('parent',h_fig);
plot(h_axes,-lags*dt,ref_corr)
xlabel(h_axes,'GPS lag (sec)');
ylabel(h_axes,'Cross correlation');
grid(h_axes,'on');
fig_fn = ct_filename_ct_tmp(param,'','check_surface','gps');
fprintf('Saving %s\n', fig_fn);
saveas(h_fig,[fig_fn '.fig']);
saveas(h_fig,[fig_fn '.jpg']);


% =====================================================================
%% Check surface: Nyquist Zone
% =====================================================================
[~,radar_type] = ct_output_dir(param.radar_name);
if strcmpi(radar_type,'deramp')
  
  param.load.imgs = {[1 1]};
  [wfs,~] = data_load_wfs(param,records);
  
  % Calculate Nyquist zone based on above ground level (AGL) altitude
  wf = 1;
  nz_twtt = param.radar.fs/2 / wfs(wf).chirp_rate;
  
  nz = floor((surf.dem_twtt+wfs(wf).Tsys-wfs(wf).t_ref) / nz_twtt);
  interp_nz = isnan(nz);
  nz = round(interp_finite(nz,NaN));
  
  default_nz = mode(nz(nz <=  max(param.radar.nz_valid) & ~isnan(nz)));
  
  nz(isnan(nz)) = default_nz;
  nz(nz<min(param.radar.nz_valid)) = min(param.radar.nz_valid);
  nz(nz>max(param.radar.nz_valid)) = max(param.radar.nz_valid);
  
  if isfield(records.settings,'nyquist_zone')
    original_nz = records.settings.nyquist_zone;
  else
    original_nz = nan(size(records.gps_time));
  end
  records.settings.nyquist_zone = interp1(layers(radar_idx).gps_time,nz,records.gps_time,'nearest','extrap');
  
  if param.check_surface.save_records_en
    save(records_fn,'-append','-struct','records','settings');
    create_records_aux_files(records_fn,false);
  end
  
  clf(h_fig);
  h_axes = axes('parent',h_fig);
  plot(h_axes,layers(radar_idx).gps_time - origin, nz,'o');
  hold(h_axes,'on');
  plot(layers(radar_idx).gps_time(~interp_nz) - origin, nz(~interp_nz), 'r.');
  plot(h_axes,records.gps_time - origin, original_nz,'g-');
  title(h_axes,'Nyquist Zone');
  ylim(h_axes,[-0.1+min(param.radar.nz_valid), max(param.radar.nz_valid)+0.1]);
  
  fig_fn = ct_filename_ct_tmp(param,'','check_surface','nz');
  fprintf('Saving %s\n', fig_fn);
  saveas(h_fig,[fig_fn '.fig']);
  saveas(h_fig,[fig_fn '.jpg']);
else
  default_nz = NaN;
end

% =====================================================================
%% Check surface: Text file
% =====================================================================
if strcmpi(radar_type,'deramp')
  wf = 1;
  BW = diff(param.radar.wfs(wf).BW_window);
  dt = 1/BW;
  t_ref_new = param.radar.wfs(wf).t_ref + round(nanmedian(twtt_error)/dt)*dt;
else
  t_ref_new = 0;
end

txt_fn = [ct_filename_ct_tmp(param,'','check_surface','time') '.txt'];
fprintf('Saving %s\n', txt_fn);
fid = fopen(txt_fn,'wb');
fprintf(fid,'%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%.1f\t%.0f\t%.12g\t%s\n', ...
  param.day_seg, 1e9*mean_offset, ...
  1e9*nanmedian(twtt_error), ...
  1e9*nanstd(twtt_error), ...
  1e9*nanmax(abs(twtt_error-mean_offset)), ...
  1e9*nanmean(twtt_error_all), ...
  1e9*nanmedian(twtt_error_all), numel(recs), -lags(peak_idx)*dt, default_nz, t_ref_new, dem_source);
fprintf(1,'%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%.1f\t%.0f\t%.12g\t%s\n', ...
  param.day_seg, 1e9*mean_offset, ...
  1e9*nanmedian(twtt_error), ...
  1e9*nanstd(twtt_error), ...
  1e9*nanmax(abs(twtt_error-mean_offset)), ...
  1e9*nanmean(twtt_error_all), ...
  1e9*nanmedian(twtt_error_all), numel(recs), -lags(peak_idx)*dt, default_nz, t_ref_new, dem_source);
fclose(fid);

% =====================================================================
%% Check surface: Tsys Refinement
% =====================================================================
if refine_Tsys_en
  % Uses specular leads to check Tsys
  deconv_fn = fullfile(ct_filename_out(param, 'noise', '', 1), sprintf('specular_%s.mat', param.day_seg));
  fprintf('Loading %s (%s)\n', deconv_fn, datestr(now))
  spec = load(deconv_fn);
  spec_gps_time = spec.gps_time(spec.peakiness > 40);
  fprintf('  %d specularity records\n', length(spec_gps_time));
  records = load(ct_filename_support(param,'','records'));
  along_track = geodetic_to_along_track(records.lat,records.lon);
  load(ct_filename_support(param,'','frames'));
  
  record = [];
  spec_frm = [];
  spec_gps = [];
  spec_fn = {};
  spec_radar_surf = [];
  spec_radar_peak = [];
  spec_atm = [];
  spec_atm_slope = [];
  spec_radar_slope = [];
  old_spec_fn = '';
  for idx = 1:length(spec_gps_time)
    record(idx) = find(records.gps_time > spec_gps_time(idx), 1);
    spec_gps(idx) = records.gps_time(record(idx));
    spec_frm(idx) = find(frames.frame_idxs < record(idx), 1, 'last');
    
    spec_fn{idx} = fullfile(ct_filename_out(param, 'deconv', ''), sprintf('Data_%s_%03d.mat', param.day_seg, spec_frm(idx)));
    if ~exist(spec_fn{idx},'file')
      spec_radar_peak(idx) = NaN;
      spec_radar_surf(idx) = NaN;
      spec_atm(idx) = NaN;
      spec_atm_slope(idx) = NaN;
      spec_radar_slope(idx) = NaN;
      continue;
    end
    if ~strcmpi(old_spec_fn,spec_fn{idx})
      mdata = load(spec_fn{idx});
    end
    old_spec_fn = spec_fn{idx};
    [time_offset,rline] = min(abs(mdata.GPS_time - spec_gps_time(idx)));
    if time_offset > 1
      warning('Time offset too large (%d sec)\n', time_offset);
    end
    
    [~,spec_radar_peak(idx)] = max(mdata.Data(:,rline));
    spec_radar_peak(idx) = mdata.Time(spec_radar_peak(idx));
    spec_radar_surf(idx) = mdata.Surface(rline);
    warning off 'MATLAB:interp1:NaNinY'
    spec_atm(idx) = interp1(layers(radar_idx).gps_time, layers(ref_idx).twtt_ref, spec_gps_time(idx));
    warning on 'MATLAB:interp1:NaNinY'
    
    % Find the time 50 m before and after
    slope_along_track = 50;
    start_record = record(idx);
    while start_record > 1 && along_track(record(idx)) - along_track(start_record) < slope_along_track
      start_record = start_record - 1;
    end
    stop_record = record(idx);
    while stop_record < length(along_track) && along_track(stop_record) - along_track(record(idx)) < slope_along_track
      stop_record = stop_record + 1;
    end
    dx = along_track(stop_record) - along_track(start_record);
    elev_time = records.elev(start_record:stop_record);
    delev = (elev_time - mean(elev_time)) / (c/2);
    atm_times = layers(ref_idx).twtt_ref(start_record:stop_record) - delev;
    radar_times = records.surface(start_record:stop_record) - delev;
    spec_atm_slope(idx) = max(atm_times) - min(atm_times);
    spec_radar_slope(idx) = max(radar_times) - min(radar_times);
    
    if debug_level > 0 && ~isnan(spec_atm(idx))
      fprintf('%d of %d (%s)\n', idx, length(spec_gps_time), datestr(now));
      spec_error = spec_atm(idx) - spec_radar_surf(idx)
      figure(1); clf;
      imagesc(mdata.GPS_time, mdata.Time, lp(mdata.Data))
      hold on
      plot(spec_gps_time(idx), spec_radar_peak(idx), 'bx');
      plot(spec_gps_time(idx), spec_radar_surf(idx), 'kx');
      plot(layers(radar_idx).gps_time, surf.dem_twtt, 'k-');
      plot([spec_gps_time(idx) spec_gps_time(idx)], mdata.Time([1 end]), 'k-');
      hold off;
      xlim(mdata.GPS_time(rline) + [-3 3]);
      ylim(spec_atm(idx) + [-10e-9 10e-9]);
      keyboard
    end
  end
  
  save(deconv_fn,'-append','spec_gps_time','spec_radar_surf','spec_radar_peak','spec_atm','spec_frm','spec_atm_slope','spec_radar_slope');
  
  spec_error = spec_radar_surf - spec_atm;
  spec_error = spec_error(~isnan(spec_error));
  fprintf('  Median error (ns)\t%.2f\n', median(spec_error)*1e9);
  figure(1); clf;
  title(param.day_seg,'interpreter','none')
  plot(sort(spec_error))
  if ~isempty(spec_error)
    ylim(median(spec_error) + [-3e-9 3e-9])
  end
  grid on;
  drawnow;
end

% =====================================================================
%% Clean up
% =====================================================================

fprintf('Done %s (%s)\n', param.day_seg, datestr(now));

if debug_level > 0
  keyboard
end
