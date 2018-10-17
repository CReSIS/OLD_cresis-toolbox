% script fmcw_set_nyquist_zone_Tsys_GPSoffset
%
% Updates the records and frames files b
% 1. Loads coincident LIDAR data if it exists
% 2. Loads DTU sea surface DEM and arctic/antarctica land DEM, combines
%    these two DEMS taking land DEM over sea surface DEM.
% 3. Combines LIDAR data and DEM data, taking LIDAR data over DEM data
% 4. Estimates Tsys error by comparing radar surface from records file and
%    LIDAR. This error should be subtracted from param.radar.wfs.Tsys.
% 5. Estimates GPS offset by comparing radar surface and LIDAR. This offset
%    should be added to param.vectors.gps.time_offset.
% 6. Uses these data to determine the Nyquist zone. Suggests a default
%    Nyquist zone for the segment and then sets frames file and records
%    file based on this. The second decimal mask in frames.proc_mode is
%    also set to one for frames that will be outside max_nyquist_zone.
%
% See run_fmcw_set_nyquist_zone_Tsys_GPSoffset.m for how to run.
%
% Author: John Paden

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

% =========================================================================
%% Automated Section
% =========================================================================

physical_constants();

%% Load in ocean mask, land DEM, and sea surface DEM
global load_surface_land_dems_finished;
global load_surface_land_dems_day_seg;
global ocean_shp_all;
global ocean_shp_bb;
global land_surface;
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
    hold on;
    plot(records.x,records.y,'r','LineWidth',2)
    set(gca,'YDir','normal');
    xlabel('X (m)');
    ylabel('Y (m)');
    legend('Flightline','location','best');
  end
end

%% Load Data

layer_params = [];

ref_idx = 2; % Make the radar the reference "slow" time axis
idx = 0;

idx = idx + 1;
layer_params(idx).name = 'surface';
layer_params(idx).source = 'lidar';
layer_params(idx).lidar_source = lidar_source;
idx = idx + 1;
layer_params(idx).name = 'surface';
layer_params(idx).source = 'layerdata';

global gRadar;
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  param = merge_structs(param,gRadar);
  param.cmd.frms = [];
  
  %fprintf('\nSet nyquist zone %s (%s)\n', param.day_seg, datestr(now));
  
  % Load radar surface from records and lidar surface
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

  % Interpolate LIDAR onto RADAR time
  
  % Interpolate onto reference
  lay_idxs = [1:ref_idx-1 ref_idx+1:length(layers)];
  
  layers(ref_idx).twtt_ref = layers(ref_idx).twtt;
  
  master = [];
  master.GPS_time = layers(ref_idx).gps_time;
  master.Latitude = layers(ref_idx).lat;
  master.Longitude = layers(ref_idx).lon;
  master.Elevation = layers(ref_idx).elev;
  for lay_idx = lay_idxs
    ops_layer = [];
    ops_layer{1}.gps_time = layers(lay_idx).gps_time;
    ops_layer{1}.type = layers(lay_idx).type;
    ops_layer{1}.quality = layers(lay_idx).quality;
    ops_layer{1}.twtt = layers(lay_idx).twtt;
    ops_layer{1}.type(isnan(ops_layer{1}.type)) = 2;
    ops_layer{1}.quality(isnan(ops_layer{1}.quality)) = 1;
    lay = opsInterpLayersToMasterGPSTime(master,ops_layer,lidar_interp_gaps_dist);
    layers(lay_idx).twtt_ref = lay.layerData{1}.value{2}.data;
  end
  
  if ~use_lidar_data
    % Set all LIDAR data to NaN
    lay_idx = 1;
    layers(lay_idx).twtt_ref(:) = NaN;
    layers(lay_idx).twtt(:) = NaN;
  end
  
  layers(ref_idx).twtt = layers(ref_idx).twtt * radar_twtt_ratio;
  layers(ref_idx).twtt = layers(ref_idx).twtt + radar_twtt_offset;
  
  lay_idx = 1;
  twtt_error = layers(lay_idx).twtt_ref - layers(ref_idx).twtt;
  twtt_error(abs(twtt_error) > 100e-9) = NaN;
  if debug_level > 0
    figure(1); clf;
    plot(layers(ref_idx).gps_time, twtt_error*1e9)
    xlabel('GPS time');
    ylabel('TWTT error (ns)');
  end
  
  mean_offset = nanmean(twtt_error);
  fprintf('%s\tTsys error (ns):\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
    param.day_seg, 1e9*mean_offset, ...
    1e9*nanmedian(twtt_error), ...
    1e9*nanstd(twtt_error), ...
    1e9*nanmax(abs(twtt_error- mean_offset)));
  
  % Find the longest contiguous section of small twtt_error
  mask = ~isnan(twtt_error);
  mask_length = zeros(size(mask));
  mask_length(1) = mask(1);
  for idx=2:length(mask)
    mask_length(idx) = mask(idx)*mask_length(idx-1) + mask(idx);
  end
  [corr_len,corr_idx] = max(mask_length);
  gpstime_coords = layers(ref_idx).gps_time(corr_idx+[-corr_len+1,0]);
  recs = find(layers(ref_idx).gps_time >= min(gpstime_coords) ...
    & layers(ref_idx).gps_time <= max(gpstime_coords));
  
  
  if combine_elev_lidar_en
    % Combine surface/lidar twtt with elevation to fill in gaps in
    % surface/lidar twtt
    % 1. Always use LIDAR when it is available
    % 2. Use elevation when LIDAR is not available, but correct it based on available LIDAR data
    elev_twtt = layers(2).elev/(3e8/2);
    twtt_error = layers(1).twtt_ref - elev_twtt;
    if ~isempty(regexpi(param.cmd.notes,'sea'))
      twtt_error = interp_finite(twtt_error,NaN);
    else
      good_mask = isfinite(twtt_error);
      if sum(good_mask) >= 2
        twtt_error(~good_mask) = interp1(find(good_mask),twtt_error(good_mask),find(~good_mask),'linear');
      end
    end
    combined_twtt = layers(1).twtt_ref;
    combined_twtt(isnan(combined_twtt)) = elev_twtt(isnan(combined_twtt)) ...
      + twtt_error(isnan(combined_twtt));
  else
    combined_twtt = layers(1).twtt_ref;
  end
  
  if combine_surface_land_dems
    % Interpolate to find sea surface elevation along flight path
    records.sea_dem = interp2(sea_surface.lon,sea_surface.lat,sea_surface.elev,mod(layers(2).lon,360),layers(2).lat);
    
    % Interpolate to find land surface elevation along flight path
    [records.x,records.y] = projfwd(land_surface.proj,layers(2).lat,layers(2).lon);
    if 0
      % Debug Plot
      figure(1); clf;
      land_surface.x = land_surface.R(3,1) + land_surface.R(2,1)*(0:size(land_surface.dem,2)-1);
      land_surface.y = land_surface.R(3,2) + land_surface.R(1,2)*(0:size(land_surface.dem,1)-1);
      imagesc(land_surface.x,land_surface.y,land_surface.dem)
      set(gca,'YDir','normal');
      hold on;
      plot(records.x,records.y);
    end
    
    if length(land_surface.x) > 2 && length(land_surface.y) > 2
      records.land_dem = interp2(land_surface.x,land_surface.y,land_surface.dem,records.x,records.y);
      records.land_dem(records.land_dem==land_surface.bad_value) = NaN;
    else
      records.land_dem = nan(size(records.x));
    end
    
    if 0
      % Debug Plot
      figure(20); clf;
      plot(records.land_dem+records.sea_dem,'g')
      hold on
      plot(layers(2).elev,'r')
      plot(records.sea_dem,'b')
      plot(layers(2).elev - layers(1).twtt_ref*3e8/2,'m')
      plot(layers(2).elev - layers(2).twtt_ref*3e8/2,'k')
    end
    
    % Merge land and surface DEMs
    surf_elev = records.land_dem;
    surf_elev(isnan(surf_elev)) = records.sea_dem(isnan(surf_elev));
    
    % Combine surface/lidar twtt with elevation to fill in gaps in
    % surface/lidar twtt
    % 1. Always use LIDAR when it is available
    % 2. Use elevation when LIDAR is not available
    elev_twtt = (layers(2).elev - surf_elev)/(3e8/2);

    combined_twtt = merge_vectors(layers(1).twtt_ref,elev_twtt);
    
  end
  
  if debug_level > 0
    %% DEBUG: Plot reference layer and then interpolated layers
    figure(2); clf;
    debug_gps_offset = 0;
    debug_Tsys_offset = 0;
    debug_Tsys_ratio = 1;
    origin = layers(ref_idx).gps_time(1);
    %debug_gps_offset = -107.5; figure(2); hold on; plot(layers(ref_idx).gps_time + debug_gps_offset,layers(ref_idx).twtt + debug_Tsys_offset, '.')
    plot(layers(ref_idx).gps_time + debug_gps_offset - origin, ...
      layers(ref_idx).twtt*debug_Tsys_ratio + debug_Tsys_offset, 'b')
    hold on;
    plot(layers(ref_idx).gps_time - origin,combined_twtt, 'g')
    for lay_idx = lay_idxs
      hold on;
      plot(layers(ref_idx).gps_time - origin, layers(lay_idx).twtt_ref, 'r')
      plot(layers(ref_idx).gps_time(recs) - origin, layers(lay_idx).twtt_ref(recs), 'r','LineWidth',4)
      grid on
    end
    title('Are the Radar and LIDAR synchronized?');
    xlabel('Relative GPS time (sec)');
    ylabel('TWTT (sec)');
    legend('RADAR','LIDAR+ELEV','LIDAR','location','best');
  end
  
  if 0
    % Correlate to determine GPS offsets
    figure(2);
    fprintf('Pick two points on figure 1 to constrain the cross correlation to a section of lidar elevation and radar elevation that are correct (but possibly offset in radar.wfs.Tsys and gps.time_offset).\n');
    fprintf('For each click hold mouse button still after click until cross hairs re-appear\n');
    [gpstime_coords,tmp] = ginput(2);
  end
  
  recs = find(layers(ref_idx).gps_time >= min(gpstime_coords) ...
    & layers(ref_idx).gps_time <= max(gpstime_coords));
  
  % Uniformly time sample the two signals
  dt = median(diff(layers(ref_idx).gps_time(recs)));
  t0 = layers(ref_idx).gps_time(recs(1));
  gps_time = t0:dt:layers(ref_idx).gps_time(recs(end));
  ref_layer = interp1(layers(ref_idx).gps_time, ...
    layers(ref_idx).twtt, gps_time);
  if length(layers(1).gps_time) < 2
    lidar_layer = NaN(size(gps_time));
  else
    lidar_layer = interp1(layers(1).gps_time, ...
      layers(1).twtt, gps_time);
  end
  
  ref_layer = ref_layer - mean(ref_layer);
  lidar_layer = lidar_layer - mean(lidar_layer);
  
  if 0
    [lidar_corr,lags] = xcorr(ref_layer,lidar_layer);
  else
    max_lag = round(10/dt);
    lags = -max_lag:max_lag;
    lidar_corr = zeros(1,length(lags));
    for lag_idx = 1:length(lags)
      lag = lags(lag_idx);
      lidar_corr(lag_idx) = sum(abs(ref_layer(1+max_lag:end-max_lag)-lidar_layer(lag + [1+max_lag:end-max_lag])));
    end
    lidar_corr = 1./lidar_corr;
  end
  [peak_val,peak_idx] = max(lidar_corr);
  
  if debug_level > 0
    figure(3); clf;
    plot(-lags*dt,lidar_corr)
    xlabel('GPS lag (sec)');
    ylabel('Cross correlation');
    grid on;
  end
  
  fprintf('%s\tGPS offset (sec) using %d records:\t%.1f\n', param.day_seg, numel(recs), -lags(peak_idx)*dt);
  
  if save_records_en
    nz = fmcw_set_nyquist_zone_from_elev(param, layers(2).gps_time, combined_twtt, max_nz);
    
    if debug_level > 0
      figure(4); clf;
      plot(layers(ref_idx).gps_time - origin, nz);
      title('Nyquist Zone');
      ylim([-0.1, max_nz+0.1]);
      figure(2); h_axes=gca;
      figure(4); h_axes(2)=gca;
      linkaxes(h_axes,'x');
    end
  end
  
  if refine_Tsys_en
    
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
      spec_atm(idx) = interp1(layers(2).gps_time, layers(1).twtt_ref, spec_gps_time(idx));
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
      atm_times = layers(1).twtt_ref(start_record:stop_record) - delev;
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
        plot(layers(1).gps_time, layers(1).twtt, 'k-');
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
  
  if debug_level > 0
    keyboard
  end
end

return











