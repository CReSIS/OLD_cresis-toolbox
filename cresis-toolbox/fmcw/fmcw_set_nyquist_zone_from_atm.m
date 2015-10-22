
% script fmcw_set_nyquist_zone_from_atm
%
% Author: John Paden

%% User Settings
% param_fn = ct_filename_param('kuband_param_2009_Antarctica_DC8.xls'); max_nz = 2;
% param_fn = ct_filename_param('snow_param_2009_Greenland_P3.xls'); max_nz = 2;
% param_fn = ct_filename_param('snow_param_2010_Greenland_DC8.xls'); max_nz = 2;
% param_fn = ct_filename_param('snow_param_2011_Greenland_P3.xls'); max_nz = 2;
param_fn = ct_filename_param('snow_param_2012_Greenland_P3.xls'); max_nz = 4;

combine_elev_lidar_en = true;
combine_elev_lidar_max_diff = 400e-9;
combine_gimp_and_geoid = true;

debug_level = 0; % Set to zero to remove stops

save_records_en = true;

refine_Tsys_en = false;

%lidar_interp_gaps_dist = [1000 250];
lidar_interp_gaps_dist = [150 75];

params = read_param_xls(param_fn,'','post');

% =========================================================================
%% Automated Section
% =========================================================================

physical_constants();

%% Load Data

layer_params = [];

ref_idx = 2; % Make the radar the reference "slow" time axis
idx = 0;

idx = idx + 1;
layer_params(idx).name = 'surface';
layer_params(idx).source = 'lidar';
layer_params(idx).echogram_source = 'qlook';
layer_params(idx).layerdata_source = 'layerData';
idx = idx + 1;
layer_params(idx).name = 'surface';
layer_params(idx).source = 'records';
layer_params(idx).echogram_source = 'qlook';
layer_params(idx).layerdata_source = 'layerData';

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  param = merge_structs(param,gRadar);
  param.cmd.frms = [];
  
  %fprintf('\nSet nyquist zone %s (%s)\n', param.day_seg, datestr(now));
  
  layers = opsLoadLayers(param,layer_params);
  
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
  
  if combine_gimp_and_geoid
    % Interpolate to find geoid along flight path
    fn = ct_filename_gis([],'world\egm96_geoid\WW15MGH.DAC');
    [lat,lon,egm96] = egm96_loader(fn);
    records.geoid_elev = interp2(lon,lat,egm96,mod(layers(2).lon,360),layers(2).lat);
    
    % Interpolate to find GIMP elevation along flight path
    [gimp_dem, R, tmp] = geotiffread('/cresis/snfs1/dataproducts/GIS_data/greenland/DEM/GIMP/gimpdem_90m.tif');
    gimp_dem = double(gimp_dem);
    proj = geotiffinfo('/cresis/snfs1/dataproducts/GIS_data/greenland/DEM/GIMP/gimpdem_90m.tif');
    
    [records.x,records.y] = projfwd(proj,layers(2).lat,layers(2).lon);
    if 0
      figure(1); clf;
      imagesc(gimp_x,gimp_y,gimp_dem)
      set(gca,'YDir','normal');
      hold on;
      plot(records.x,records.y);
    end
    records.gimp_dem = interp2(gimp_dem,(records.x-R(3,1))/R(2,1)+1,(records.y-R(3,2))/R(1,2)+1);

    if 0
      figure(2); clf;
      plot(records.gimp_dem,'g')
      hold on
      plot(layers(2).elev,'r')
      plot(records.geoid_elev,'b')
      plot(layers(2).elev - layers(1).twtt_ref*3e8/2,'m')
    end
    
    % Merge GIMP and Geoid
    surf_elev = records.gimp_dem;
    surf_elev(isnan(surf_elev)) = records.geoid_elev(isnan(surf_elev));
    
    % Combine surface/lidar twtt with elevation to fill in gaps in
    % surface/lidar twtt
    % 1. Always use LIDAR when it is available
    % 2. Use elevation when LIDAR is not available
    elev_twtt = (layers(2).elev - surf_elev)/(3e8/2);

    combined_twtt = layers(1).twtt_ref;
    combined_twtt(isnan(combined_twtt)) = elev_twtt(isnan(combined_twtt));
    
  end
  
  if debug_level > 0
    % Plot reference layer and then interpolated layers
    figure(2); clf;
    debug_gps_offset = 0;
    debug_Tsys_offset = 0;
    %debug_gps_offset = -107.5; figure(2); hold on; plot(layers(ref_idx).gps_time + debug_gps_offset,layers(ref_idx).twtt + debug_Tsys_offset, '.')
    plot(layers(ref_idx).gps_time + debug_gps_offset, ...
      layers(ref_idx).twtt + debug_Tsys_offset, '.')
    hold on;
    plot(layers(ref_idx).gps_time,combined_twtt, 'g')
    for lay_idx = lay_idxs
      hold on;
      plot(layers(ref_idx).gps_time, layers(lay_idx).twtt_ref, 'r')
      plot(layers(ref_idx).gps_time(recs), layers(lay_idx).twtt_ref(recs), 'r.')
      grid on
    end
    title('TWTT');
    legend('RADAR','LIDAR+ELEV','LIDAR');
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
  
  [lidar_corr,lags] = xcorr(ref_layer,lidar_layer);
  lag_times = -120:120;
  
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
    nz = fmcw_set_nyquist_zone_from_elev(param, combined_twtt, max_nz);
    
    if debug_level > 0
      figure(4); clf;
      plot(nz);
      title('Nyquist Zone');
      ylim([-0.1, max_nz+0.1])
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











