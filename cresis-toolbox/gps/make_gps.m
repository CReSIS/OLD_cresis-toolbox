% script make_gps
%
% Makes GPS files (called from each mission specific routine)
% Variables set up in make_gps_SEASON_SOURCE
%
% in_fn = 1 by N cell array of strings containing filenames
%  Full path of input GPS file
%  Can also do a cell array of strings in place of one of the strings and
%  this will cause all the file to be concatenated into one output GPS file
% file_type = 1 by N cell array of strings
%  File format for each instance (e.g. 'Applanix', 'ATM', 'Litton', etc)
%  Determines which read_gps_? function will be called
% params = 1 by N cell array of structures
%  Parameters for each instance
%  If the option to pass a cell array of filename strings for a single output GPS
%  file is used, then the params structure must also be cell array of
%  structures which correspond to each of the files in the cell array of
%  filename strings.
% gps_path = string
%  Output directory
% out_fn = 1 by N cell array of strings
%  Output file name for each instance
% gps_source = string describing GPS processing location and the
%  version of the processed file. This is done by placing a hyphen
%  after the location identifier. The location identifier is used
%  with the season ID and the lever arm function to determine which
%  lever arm to use. Examples:
%  gravimeter-field: processed to gravimeter INS, field version of GPS/INS
%    product
%  atm-field: processed to ATM DGPS antenna, field version of GPS/INS
%    product
%  atm-final_20110826: processed to ATM DGPS antenna, final version of GPS/INS
%    product as of 20110826
%
%
% Author: John Paden

% ======================================================================
% Read and translate files according to user settings
% ======================================================================

for file_idx = 1:length(in_fns)
  in_fn = in_fns{file_idx};
  out_fn = fullfile(gps_path,out_fns{file_idx});
  if isempty(in_fn)
    error('No input GPS files found in in_fns{%d} variable. Usually this is because the file path is not setup.\n', file_idx);
  end
  if iscell(in_fn)
    fprintf('Input files (%.1f sec)\n', toc);
    for in_fn_idx = 1:length(in_fn)
      fprintf('  %s\n', in_fn{in_fn_idx});
    end
  else
    fprintf('Input file %s (%.1f sec)\n', in_fn, toc);
  end
  fprintf('  Output file %s\n', out_fn);
  
  if exist('sync_flag','var') && sync_flag{file_idx}
    %% Load Radar-Synchronization NMEA format 3 sync files (mcrds, accum2)
    if ~iscell(sync_fns{file_idx})
      sync_fns{file_idx} = {sync_fns{file_idx}};
    end
    clear sync_gps;
    for sync_fn_idx = 1:length(sync_fns{file_idx})
      fprintf('  Sync file %s\n', sync_fns{file_idx}{sync_fn_idx});
      if sync_fn_idx == 1
        sync_gps = read_gps_nmea(sync_fns{file_idx}{sync_fn_idx},sync_params{file_idx});
      else
        sync_gps_tmp = read_gps_nmea(sync_fns{file_idx}{sync_fn_idx},sync_params{file_idx});
        sync_gps.gps_time = [sync_gps.gps_time, sync_gps_tmp.gps_time];
        sync_gps.lat = [sync_gps.lat, sync_gps_tmp.lat];
        sync_gps.lon = [sync_gps.lon, sync_gps_tmp.lon];
        sync_gps.elev = [sync_gps.elev, sync_gps_tmp.elev];
        sync_gps.comp_time = [sync_gps.comp_time, sync_gps_tmp.comp_time];
        if isfield(sync_gps,'radar_time')
          sync_gps.radar_time = [sync_gps.radar_time, sync_gps_tmp.radar_time];
        end
      end
    end
  end
  
  %% Load GPS files
  if ischar(in_fn)
    in_fn = {in_fn};
  end
  if isstruct(params{file_idx})
    params{file_idx} = repmat({params{file_idx}},size(in_fn));
  end
  clear gps;
  for in_fn_idx = 1:length(in_fn)
    if strcmpi(file_type{file_idx},'Applanix')
      gps_tmp = read_gps_applanix(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'awi_netcdf')
      gps_tmp = read_gps_netcdf(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'awi_netcdf+awi_netcdf')
      gps_tmp = read_gps_netcdf(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'cresis')
      gps_tmp = read_gps_cresis(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'DMSraw')
      gps_tmp = read_gps_dmsraw(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'General_ASCII')
      gps_tmp = read_gps_general_ascii(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'Litton')
      gps_tmp = read_ins_litton(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'Litton_DGPS')
      gps_tmp = read_gps_litton(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'NMEA')
      gps_tmp = read_gps_nmea(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'Novatel')
      gps_tmp = read_gps_novatel(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'Reveal')
      gps_tmp = read_gps_reveal(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'Novatel_RPYGGA')
      gps_tmp = read_gps_novatel_rpygga(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif any(strcmpi(file_type{file_idx},{'Traj','Traj+Litton','Traj+Litton_DGPS','Traj+General_ASCII'}))
      gps_tmp = read_gps_traj(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'TXT')
      gps_tmp = read_gps_txt(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    elseif strcmpi(file_type{file_idx},'csv')
      gps_tmp = read_gps_csv(in_fn{in_fn_idx},params{file_idx}{in_fn_idx});
    else
      error('Unrecognized GPS file type %s', file_type{file_idx});
    end
    
    if in_fn_idx == 1
      gps = gps_tmp;
    else
      gps.gps_time = [gps.gps_time, gps_tmp.gps_time];
      gps.lat = [gps.lat, gps_tmp.lat];
      gps.lon = [gps.lon, gps_tmp.lon];
      gps.elev = [gps.elev, gps_tmp.elev];
      gps.roll = [gps.roll, gps_tmp.roll];
      gps.pitch = [gps.pitch, gps_tmp.pitch];
      gps.heading = [gps.heading, gps_tmp.heading];
    end
  end
  
  %% Remove records with NaN
  good_mask = ~(isnan(gps.gps_time(1:length(gps.gps_time))) | isnan(gps.lat) ...
    | isnan(gps.lon) | isnan(gps.elev) ...
    | isnan(gps.roll) | isnan(gps.pitch) | isnan(gps.heading));
  gps.gps_time = gps.gps_time(good_mask);
  gps.lat = gps.lat(good_mask);
  gps.lon = gps.lon(good_mask);
  gps.elev = gps.elev(good_mask);
  gps.roll = gps.roll(good_mask);
  gps.pitch = gps.pitch(good_mask);
  gps.heading = gps.heading(good_mask);
  gps.gps_source = gps_source{file_idx};
  
  %% Fabricating a heading now
  along_track = geodetic_to_along_track(gps.lat,gps.lon,gps.elev);
  rlines = get_equal_alongtrack_spacing_idxs(along_track,10);
  physical_constants;
  est_heading = size(gps.heading);
  clear origin heading east north;
  for rline_idx = 1:length(rlines)
    rline = rlines(rline_idx);
    if rline_idx < length(rlines)
      rline_end = rlines(rline_idx+1);
    else
      rline_end = length(along_track);
    end
    [origin(1),origin(2),origin(3)] = geodetic2ecef(gps.lat(rline)/180*pi,gps.lon(rline)/180*pi,gps.elev(rline),WGS84.ellipsoid);
    [heading(1),heading(2),heading(3)] = geodetic2ecef(gps.lat(rline_end)/180*pi,gps.lon(rline_end)/180*pi,gps.elev(rline_end),WGS84.ellipsoid);
    heading = heading - origin;
    % Determine east vector
    [east(1) east(2) east(3)] = lv2ecef(1,0,0,gps.lat(rline)/180*pi,gps.lon(rline)/180*pi,gps.elev(rline),WGS84.ellipsoid);
    east = east - origin;
    % Determine north vector
    [north(1) north(2) north(3)] = lv2ecef(0,1,0,gps.lat(rline)/180*pi,gps.lon(rline)/180*pi,gps.elev(rline),WGS84.ellipsoid);
    north = north - origin;
    % Determine heading
    est_heading(rline:rline_end) = atan2(dot(east,heading),dot(north,heading));
  end
  
  %% Load INS data for special case where it is separate from GPS
  if any(strcmpi(file_type{file_idx},{'Traj+Litton','Traj+Litton_DGPS','Traj+General_ASCII','awi_netcdf+awi_netcdf'}))
    if ischar(in_fns_ins{file_idx})
      in_fns_ins{file_idx} = {in_fns_ins{file_idx}};
    end
    if isstruct(params_ins{file_idx})
      params_ins{file_idx} = repmat({params_ins{file_idx}},size(in_fns_ins{file_idx}));
    end
    clear ins;
    if isempty(in_fns_ins{file_idx})
      error('No input INS files found in in_fns{%d} variable. Usually this is because the file path is not setup.\n', file_idx);
    end
    for in_fn_idx = 1:length(in_fns_ins{file_idx})
      fprintf('  INS file %s\n', in_fns_ins{file_idx}{in_fn_idx});
      if strcmpi(file_type{file_idx},'Traj+Litton')
        ins_tmp = read_ins_litton(in_fns_ins{file_idx}{in_fn_idx},params_ins{file_idx}{in_fn_idx});
      elseif strcmpi(file_type{file_idx},'Traj+Litton_DGPS')
        ins_tmp = read_gps_litton(in_fns_ins{file_idx}{in_fn_idx},params_ins{file_idx}{in_fn_idx});
      elseif strcmpi(file_type{file_idx},'Traj+General_ASCII')
        ins_tmp = read_gps_general_ascii(in_fns_ins{file_idx}{in_fn_idx},params_ins{file_idx}{in_fn_idx});
      elseif  strcmpi(file_type{file_idx},'awi_netcdf+awi_netcdf')
        ins_tmp = read_gps_netcdf(in_fns_ins{file_idx}{in_fn_idx},params_ins{file_idx}{in_fn_idx});
      end
      if in_fn_idx == 1
        ins = ins_tmp;
      else
        ins.gps_time = [ins.gps_time, ins_tmp.gps_time];
        ins.lat = [ins.lat, ins_tmp.lat];
        ins.lon = [ins.lon, ins_tmp.lon];
        ins.elev = [ins.elev, ins_tmp.elev];
        ins.roll = [ins.roll, ins_tmp.roll];
        ins.pitch = [ins.pitch, ins_tmp.pitch];
        ins.heading = [ins.heading, ins_tmp.heading];
      end
    end
    
    bad_mask = ins.heading == 0;
    if sum(bad_mask) > 0
      fprintf('Found %d bad records based on heading == 0\n', sum(bad_mask));
      ins.gps_time = ins.gps_time(~bad_mask);
      ins.lat = ins.lat(~bad_mask);
      ins.lon = ins.lon(~bad_mask);
      ins.elev = ins.elev(~bad_mask);
      ins.roll = ins.roll(~bad_mask);
      ins.pitch = ins.pitch(~bad_mask);
      ins.heading = ins.heading(~bad_mask);
    end
    
    % Check for large time gaps in INS data and fill with level flight
    if 0
      figure(1); clf;
      bad_mask = diff(ins.gps_time) > 10;
      bad_mask = grow(bad_mask,1);
      time_skips = diff(ins.gps_time);
      along_track = geodetic_to_along_track(gps.lat,gps.lon);
      vel = diff(along_track) ./ diff(gps.gps_time);
      plot(gps.gps_time(1:end-1),vel)
      hold on
      plot(ins.gps_time(bad_mask),time_skips(bad_mask),'r');
      hold off;
      
      figure(2); clf;
      plot(ins.gps_time,ins.roll*180/pi,'.')
      hold on
      plot(ins.gps_time(bad_mask),time_skips(bad_mask),'r');
      hold off;
      ylim([-20 20])
      grid on
      
      figure(3); clf;
      plot(gps.gps_time,est_heading*180/pi,'g.');
      hold on
      plot(ins.gps_time,ins.heading*180/pi,'.')
      plot(ins.gps_time(bad_mask),time_skips(bad_mask),'r');
      hold off;
      grid on
    end
    
    bad_idxs = find(diff(ins.gps_time) > 20);
    if ins.gps_time(1) - gps.gps_time(1) > 10
      bad_idxs = [0 bad_idxs];
    end
    if gps.gps_time(end) - ins.gps_time(end) > 10
      bad_idxs = [bad_idxs length(ins.gps_time)];
    end
    for bad_idxs_idx = 1:length(bad_idxs)
      idx = bad_idxs(bad_idxs_idx);
      
      % Find gps indices in this gap
      if idx == 0
        new_gps_times = gps.gps_time(1):ins.gps_time(idx+1)-10;
        edge_gps_times = [gps.gps_time(1) ins.gps_time(idx+1)];
        est_heading_error = [0 ins.heading(idx+1) - interp1(gps.gps_time,est_heading,ins.gps_time(idx+1))];
      elseif idx == length(ins.gps_time)
        new_gps_times = ins.gps_time(idx)+10:gps.gps_time(end);
        edge_gps_times = [ins.gps_time(idx) gps.gps_time(end)];
        est_heading_error = [ins.heading(idx) - interp1(gps.gps_time,est_heading,ins.gps_time(idx)) 0];
      else
        new_gps_times = ins.gps_time(idx)+10:ins.gps_time(idx+1)-10;
        edge_gps_times = ins.gps_time(idx:idx+1);
        est_heading_error = ins.heading(idx:idx+1) - interp1(gps.gps_time,est_heading,ins.gps_time(idx:idx+1));
      end
      bad_idxs(bad_idxs_idx+1:end) = bad_idxs(bad_idxs_idx+1:end) + length(new_gps_times);
      fprintf('  Inserting %d level flight records due to missing INS data\n', length(new_gps_times));
      
      ins.roll = [ins.roll(1:idx) zeros(size(new_gps_times)) ins.roll(idx+1:end)];
      ins.pitch = [ins.pitch(1:idx) zeros(size(new_gps_times)) ins.pitch(idx+1:end)];
      ins.lat = [ins.lat(1:idx) NaN*zeros(size(new_gps_times)) ins.lat(idx+1:end)];
      ins.lon = [ins.lon(1:idx) NaN*zeros(size(new_gps_times)) ins.lon(idx+1:end)];
      ins.elev = [ins.elev(1:idx) NaN*zeros(size(new_gps_times)) ins.elev(idx+1:end)];
      
      inserted_heading = interp1(gps.gps_time,est_heading,new_gps_times) ...
        + interp1(edge_gps_times,est_heading_error,new_gps_times);
      ins.heading = [ins.heading(1:idx) inserted_heading ins.heading(idx+1:end)];
      
      ins.gps_time = [ins.gps_time(1:idx) new_gps_times ins.gps_time(idx+1:end)];
      
    end
    
    if 0
      figure(4); clf;
      plot(gps.gps_time,est_heading*180/pi,'g.');
      hold on
      plot(ins.gps_time,ins.heading*180/pi,'.')
      hold off;
      grid on
      keyboard
    end
    
    % Insert filler points from gps (ensures that sparse ins data does not
    % create a sparse final dataset when gps data is available)
    new_gps_time = union(gps.gps_time,ins.gps_time);
    ins.roll = interp1(ins.gps_time,ins.roll,new_gps_time);
    ins.pitch = interp1(ins.gps_time,ins.pitch,new_gps_time);
    ins.heading = interp1(ins.gps_time,ins.heading,new_gps_time);
    ins.gps_time = new_gps_time;
    
    ins.lat = interp1(gps.gps_time,gps.lat,ins.gps_time);
    ins.lon = interp1(gps.gps_time,gps.lon,ins.gps_time);
    ins.elev = interp1(gps.gps_time,gps.elev,ins.gps_time);
    gps = ins;
    
    %% Remove records with NaN
    good_mask = ~(isnan(gps.gps_time) | isnan(gps.lat) ...
      | isnan(gps.lon) | isnan(gps.elev) ...
      | isnan(gps.roll) | isnan(gps.pitch) | isnan(gps.heading));
    gps.gps_time = gps.gps_time(good_mask);
    gps.lat = gps.lat(good_mask);
    gps.lon = gps.lon(good_mask);
    gps.elev = gps.elev(good_mask);
    gps.roll = gps.roll(good_mask);
    gps.pitch = gps.pitch(good_mask);
    gps.heading = gps.heading(good_mask);
    gps.gps_source = gps_source{file_idx};
  end

  %% Check that GPS time is monotonically increasing
  make_gps_monotonic(gps);
  
  %% Fabricate a heading from the trajectory if it is all zeros
  if all(gps.heading == 0)
    warning('This file has heading(:) == 0. Using the estimated heading.');
    gps.heading = est_heading;
  end
  
  %% Check for day wraps
  % Find jumps in the GPS time that are probably due to day interval
  % 86400 seconds.
  day_jumps = find(diff(gps.gps_time) < -60000);
  if ~isempty(day_jumps)
    warning('Found a day wrap... correcting');
  end
  for jump_idx = day_jumps
    gps.gps_time(jump_idx+1:end) = gps.gps_time(jump_idx+1:end) + 86400;
  end
  if exist('sync_flag','var') && sync_flag{file_idx}
    day_jumps = find(diff(sync_gps.gps_time) < -60000);
    if ~isempty(day_jumps)
      warning('Found a day wrap in sync gps... correcting');
    end
    for jump_idx = day_jumps
      sync_gps.gps_time(jump_idx+1:end) = sync_gps.gps_time(jump_idx+1:end) + 86400;
    end
  end
  
  %% Add software revision information
  gps.sw_version = current_software_version;

  %% Add the Radar Synchronization variables for mcrds, accum2, acords
  if exist('sync_flag','var') && sync_flag{file_idx}
    % Synchronize computer time and radar time from NMEA file to gps time
    gps.sync_gps_time = sync_gps.gps_time;
    gps.sync_lat = sync_gps.lat;
    gps.sync_lon = sync_gps.lon;
    gps.sync_elev = sync_gps.elev;
    gps.comp_time = sync_gps.comp_time;
    if isfield(sync_gps,'radar_time')
      gps.radar_time = sync_gps.radar_time;
    end

    %% Fabricating a heading for sync GPS data
    along_track = geodetic_to_along_track(gps.sync_lat,gps.sync_lon,gps.sync_elev);
    rlines = get_equal_alongtrack_spacing_idxs(along_track,10);
    physical_constants;
    est_heading = size(gps.sync_lat);
    clear origin heading east north;
    for rline_idx = 1:length(rlines)
      rline = rlines(rline_idx);
      if rline_idx < length(rlines)
        rline_end = rlines(rline_idx+1);
      else
        rline_end = length(along_track);
      end
      [origin(1),origin(2),origin(3)] = geodetic2ecef(gps.sync_lat(rline)/180*pi,gps.sync_lon(rline)/180*pi,gps.sync_elev(rline),WGS84.ellipsoid);
      [heading(1),heading(2),heading(3)] = geodetic2ecef(gps.sync_lat(rline_end)/180*pi,gps.sync_lon(rline_end)/180*pi,gps.sync_elev(rline_end),WGS84.ellipsoid);
      heading = heading - origin;
      % Determine east vector
      [east(1) east(2) east(3)] = lv2ecef(1,0,0,gps.sync_lat(rline)/180*pi,gps.sync_lon(rline)/180*pi,gps.sync_elev(rline),WGS84.ellipsoid);
      east = east - origin;
      % Determine north vector
      [north(1) north(2) north(3)] = lv2ecef(0,1,0,gps.sync_lat(rline)/180*pi,gps.sync_lon(rline)/180*pi,gps.sync_elev(rline),WGS84.ellipsoid);
      north = north - origin;
      % Determine heading
      est_heading(rline:rline_end) = atan2(dot(east,heading),dot(north,heading));
    end
    gps.sync_heading = est_heading;
    
    
     % If sync gps data vectors are longer than gps vectors or sync_gps_time
    % extends beyond gps_time, merge the two data sets. Check for longitude
    % 360 deg offset and elevation differences.
    if gps.gps_time(end) < gps.sync_gps_time(end) || gps.gps_time(1) > gps.sync_gps_time(1)
      gps = merge_sync_gps(gps,2);
    end

    if isfield(gps,'radar_time')
      save(out_fn,'-v7.3','-STRUCT','gps','gps_time','lat','lon','elev','roll','pitch','heading','gps_source','sync_gps_time','sync_lat','sync_lon','sync_elev','sync_heading','comp_time','radar_time','sw_version');
    else
      save(out_fn,'-v7.3','-STRUCT','gps','gps_time','lat','lon','elev','roll','pitch','heading','gps_source','sync_gps_time','sync_lat','sync_lon','sync_elev','sync_heading','comp_time','sw_version');
    end
  else
    save(out_fn,'-v7.3','-STRUCT','gps','gps_time','lat','lon','elev','roll','pitch','heading','gps_source','sw_version');
  end
  
  if debug_level >= 2
    plot_gps(out_fn);
  end
end

return;

