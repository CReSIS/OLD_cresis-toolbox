function my_struct = sync_radar_to_gps(param,my_struct,radar_time,comp_time)
% my_struct = sync_radar_to_gps(param,my_struct,radar_time,comp_time)
%
% param = struct from read_param_xls
%  param.radar_name
%  param.season_name
%  param.day_seg
%  param.vectors.gps.fn
%  param.vectors.gps.time_offset
% my_struct = vector or records struct to add GPS fields to
% radar_time = time recorded to raw data files by radar (contents depend on radar)
%   mcrds,accum2: radar_time is a free running 64 bit counter
%   others: radar_time is utc_time_sod or seconds of day
% comp_time = computer time recorded to raw data files by radar
%   only required for mcrds and accum2
%
% my_struct = Adds lat, lon, elev, roll, pitch, heading, gps_time, gps_source
%  to the my_struct struct
%
% Works for the following radars:
%  mcrds, accum2, mcords, mcords2, mcords3, fmcw, fmcw2, fmcw3, accum
% Important usage note for all radars besides mcrds,accum2:
%  The absolute radar time is determined by taking the date from the
%  segment name (e.g. 20111023), adding in the gps time offset from
%  the vectors worksheet and the utc_time_sod read from the radar data
%  files. Data in places like McMurdo with +12 hours UTC time often
%  require the vectors gps time offset to be -86400 because the UTC time
%  will be one day behind local time (which is what should be used for
%  the segment name).
%
% Called from create_vectors* or create_records*

if ~isfield(param.vectors.gps,'fn')
  param.vectors.gps.fn = '';
end

%% Load the GPS data
gps = load(ct_filename_support(param,param.vectors.gps.fn,'gps',true));

%% Check for non-monotonically increasing gps time
if any(diff(gps.gps_time) <= 0)
  error('GPS times are not monotonically increasing');
end

%% Check for repeat values
if length(unique(gps.gps_time)) ~= length(gps.gps_time)
  error('GPS file has repeat values in it');
end

%% Check for NaN in any of the fields
if any(isnan(gps.lat)) ...
    || any(isnan(gps.lon)) ...
    || any(isnan(gps.elev)) ...
    || any(isnan(gps.roll)) ...
    || any(isnan(gps.pitch)) ...
    || any(isnan(gps.heading)) ...
    || any(isnan(gps.gps_time))
  error('GPS file has NaN');
end

if any(strcmpi(param.radar_name,{'mcrds','accum2'}))
  %% Isolate the section of radar time from gps.radar_time that will be
  % used to interpolate with.
  if strcmpi(param.season_name,'2013_Antarctica_P3') & strcmpi(param.radar_name,'accum2') & strcmpi(param.day_seg(1:8),'20131119')
    % no gps sync files, set radar_gps_time = comp_time(1) + radar_time -
    % radar_time(1) + comp_time_offset (-6*3600), the computer time offset
    % from gps time was 6 hours late
    radar_gps_time =  comp_time(1) + radar_time-radar_time(1) -6*3600 + param.vectors.gps.time_offset;
  else
    guard_time = 5;
    good_idxs = find(gps.comp_time >= comp_time(1)-guard_time ...
      & gps.comp_time <= comp_time(end)+guard_time);
    good_radar_time = gps.radar_time(good_idxs);
    good_sync_gps_time = gps.sync_gps_time(good_idxs);
    
    %% From these good indexes, remove any repeat radar times (usually caused
    % by there being more than one NMEA string every 1 PPS
    good_idxs = 1+find(diff(good_radar_time) ~= 0);
    good_radar_time = good_radar_time(good_idxs);
    good_sync_gps_time = good_sync_gps_time(good_idxs);
    
    %% Interpolate gps.sync_gps_time to radar gps_time using gps.radar_time
    % and radar_time
    radar_gps_time = interp1(good_radar_time, good_sync_gps_time, ...
      radar_time,'linear','extrap');
  end
  
  %% DO NOT Apply GPS sync correction to radar time (this is done already
  % in create_records for these radars)
  
elseif any(strcmpi(param.radar_name,{'icards'}))% there's a minor inacurracy (1e-7)of first
  % when read the csv file. This may cause radar_gps_time start earlier than gps.gps_time(first file)
  % or later than gps.gps_time(last file). This phenomenon will further
  % cause NaN when using interp1 to sync radar and gps!!This problem
  % cannot be solved in other scripts----qishi
  utc_time_sod = radar_time;
  
  %% Check for seconds of day roll over and unwrap (assume jump backward
  % of more than 23 hours is a roll over)
  wrap_idxs = find(abs(diff(utc_time_sod) + 86400) < 3600);
  for wrap_idx = wrap_idxs
    utc_time_sod(wrap_idx+1:end) = utc_time_sod(wrap_idx+1:end) + 86400;
  end
  
  %% Apply GPS sync correction to radar time
  utc_time_sod = utc_time_sod + param.vectors.gps.time_offset;
  
  %% Determine absolute radar time and convert from UTC to GPS
  year = str2double(param.day_seg(1:4));
  month = str2double(param.day_seg(5:6));
  day = str2double(param.day_seg(7:8));
  radar_gps_time = datenum_to_epoch(datenum(year,month,day,0,0,utc_time_sod)) + utc_leap_seconds(gps.gps_time(1));
  if radar_gps_time(1)<gps.gps_time(1)
    radar_gps_time=radar_gps_time(find(radar_gps_time>=gps.gps_time(1)));
  end
  if radar_gps_time(end)>gps.gps_time(end)
    radar_gps_time=radar_gps_time(find(radar_gps_time<=gps.gps_time(end)));
  end
  
else
  utc_time_sod = radar_time;
  
  %% Check for seconds of day roll over and unwrap (assume jump backward
  % of more than 23 hours is a roll over)
  wrap_idxs = find(abs(diff(utc_time_sod) + 86400) < 3600);
  for wrap_idx = wrap_idxs
    utc_time_sod(wrap_idx+1:end) = utc_time_sod(wrap_idx+1:end) + 86400;
  end
  
  %% Apply GPS sync correction to radar time
  utc_time_sod = utc_time_sod + param.vectors.gps.time_offset;
  
  %% Determine absolute radar time and convert from UTC to GPS
  year = str2double(param.day_seg(1:4));
  month = str2double(param.day_seg(5:6));
  day = str2double(param.day_seg(7:8));
  radar_gps_time = datenum_to_epoch(datenum(year,month,day,0,0,utc_time_sod)) + utc_leap_seconds(gps.gps_time(1));
end

%% Synchronize times to get positions and absolute time
my_struct.lat = double(interp1(gps.gps_time,gps.lat,radar_gps_time));
my_struct.lon = double(mod(interp1(gps.gps_time,unwrap(gps.lon/180*pi),radar_gps_time)*180/pi+180, 360)-180);
my_struct.elev = double(interp1(gps.gps_time,gps.elev,radar_gps_time));
my_struct.roll = double(interp1(gps.gps_time,gps.roll,radar_gps_time));
my_struct.pitch = double(interp1(gps.gps_time,gps.pitch,radar_gps_time));
my_struct.heading = double(mod(interp1(gps.gps_time,unwrap(gps.heading),radar_gps_time)+pi,2*pi)-pi);
my_struct.gps_time = radar_gps_time;
my_struct.gps_source = gps.gps_source;

%% Check for bad synchronization (NaN in synced GPS data)
nan_detected = false;
if any(isnan(my_struct.lat))
  fprintf('There are NaN in synced GPS lat.\n');
  nan_detected = true;
end
if any(isnan(my_struct.lon))
  fprintf('There are NaN in synced GPS lon.\n');
  nan_detected = true;
end
if any(isnan(my_struct.elev))
  fprintf('There are NaN in synced GPS elev.\n');
  nan_detected = true;
end
if any(isnan(my_struct.roll))
  fprintf('There are NaN in synced GPS roll.\n');
  nan_detected = true;
end
if any(isnan(my_struct.pitch))
  fprintf('There are NaN in synced GPS pitch.\n');
  nan_detected = true;
end
if any(isnan(my_struct.heading))
  fprintf('There are NaN in synced GPS heading.\n');
  nan_detected = true;
end
if any(isnan(my_struct.gps_time))
  fprintf('There are NaN in synced GPS time.\n');
  nan_detected = true;
end
if nan_detected
  warning('NaN found in GPS data');
  if any(strcmpi(param.radar_name,{'mcrds','accum2'}))
    fprintf('GPS COMP TIME: %s to %s\n', datestr(epoch_to_datenum(gps.comp_time(1))), datestr(epoch_to_datenum(gps.comp_time(end))));
    fprintf('RADAR COMP TIME: %s to %s\n', datestr(epoch_to_datenum(comp_time(1))), datestr(epoch_to_datenum(comp_time(end))));
  else
    fprintf('GPS TIME: %s to %s\n', datestr(epoch_to_datenum(gps.gps_time(1))), datestr(epoch_to_datenum(gps.gps_time(end))));
    fprintf('RADAR GPS TIME: %s to %s\n', datestr(epoch_to_datenum(my_struct.gps_time(1))), datestr(epoch_to_datenum(my_struct.gps_time(end))));
  end
  keyboard
end

return;

