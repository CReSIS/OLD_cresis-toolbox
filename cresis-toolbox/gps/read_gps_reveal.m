function gps = read_gps_reveal(in_fn,param)
% gps = read_gps_reveal(in_fn,param)
%
% Created for 2009 Antarctica DC8 field season.  Reads in text files
% produced by the Reveal system.
%   Reveal time is UTC time and is converted to GPS time.
%
% Input Args:
%   in_fn (string) input Reveal filename
%
% Output Args:
% gps = output structure with fields
%  .time = GPS time in seconds since Jan 1, 1970 epoch (sec)
%  .lat = latitude (deg)
%  .lon = longitude (deg)
%  .elev = elevation (m)
%  .roll = roll (rad)
%  .pitch = pitch (rad)
%  .heading = true heading (rad)
%
% Example:
%   fn = '/cresis/data2/MCoRDS/2009_Chile/GPS/IWG1_110209.log';
%   gps = read_gps_reveal(fn);
%   gps_plot(gps);
%   datestr(epoch_to_datenum(gps.gps_time(1)));
%   gps.utc_time = gps.gps_time - utc_leap_seconds(gps.gps_time(1))
%
% Author: William Blake, John Paden
%
% See also read_gps_*.m, gps_plot.m, gps_create.m

debug_level = 1;

[ident, ...
      raw.time_stamp, ...
      gps.lat, ...
      gps.lon, ...
      raw.gps_alt_msl, ...
      gps.elev, ...
      raw.pressure_alt, ...
      raw.radar_alt, ...
      raw.ground_speed, ...
      raw.true_air_speed, ...
      raw.indicated_air_speed, ...
      raw.mach, ...
      raw.vertical_speed, ...
      gps.heading, ...
      raw.track_angle, ...
      raw.drift_angle, ...
      gps.pitch, ...
      gps.roll, ...
      raw.slip_angle, ...
      raw.attack_angle, ...
      raw.static_air_temp, ...
      raw.dew_point, ...
      raw.total_air_temp, ...
      raw.static_pressure, ...
      raw.dynamic_pressure, ...
      raw.cabin_pressure, ...
      raw.wind_speed, ...
      raw.wind_direction, ...
      raw.vert_wind_speed, ...
      raw.solar_zenith_angle, ...
      raw.aircraft_sun_elevation, ...
      raw.sun_azimuth, ...
      raw.aircraft_sun_azimuth] ...
      = textread(in_fn, ...
 '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',','emptyvalue',NaN);


%       raw_vert_speed, ...
%       ADC_baro_alt, ...
%       potential_temp, ...
%       cabin_temp, ...
%       cabin_humidity, ...
%       acceleration_x, ...
%       acceleration_y, ...
%       acceleration_z, ...
%       ozone, ...
%       co, ...
%       ch4, ...
%       dlh_wv, ...
%       scattering, ...
%       cn, ...
%       no, ...
%       noy]

epoch = datenum(1970,1,1,0,0,0);

% Reader only outputs a string for the time stamp. This
% converts the time stamps to date numbers. The format is:
if any(raw.time_stamp{1}=='.')
  %    2009-11-02T11:35:54.011
  [year month day hour minute sec] = datevec(raw.time_stamp, 'yyyy-mm-ddTHH:MM:SS.FFF');
else
  %    2009-11-02T11:35:54
  [year month day hour minute sec] = datevec(raw.time_stamp, 'yyyy-mm-ddTHH:MM:SS');
end
n = datenum(year,month,day,hour,minute,sec).';
gps.gps_time = (n - epoch) * 86400;

if strcmpi(param.time_reference,'utc')
  % UTC time stored in file, so need to add leap seconds back in
  gps.gps_time = gps.gps_time + utc_leap_seconds(gps.gps_time(1));
else
  warning('NMEA files are usually always UTC time, but GPS time has been specified.\n');
end

gps.lat = gps.lat.';
gps.lon = gps.lon.';
gps.elev = gps.elev.';

gps.roll = gps.roll.' / 180*pi;
gps.pitch = gps.pitch.' / 180*pi;
gps.heading = gps.heading.' / 180*pi;
% Take mod 2*pi in uniform way:
gps.heading = angle(exp(j*gps.heading));

% NaN removal
good_idxs = find(isfinite(gps.lat) & isfinite(gps.lon) & isfinite(gps.elev) & gps.elev > -1e4);
gps.gps_time = gps.gps_time(good_idxs);
gps.lat = gps.lat(good_idxs);
gps.lon = gps.lon(good_idxs);
gps.elev = gps.elev(good_idxs);
gps.roll = gps.roll(good_idxs);
gps.pitch = gps.pitch(good_idxs);
gps.heading = gps.heading(good_idxs);

if isfield(param,'filter') && ~isempty(param.filter)
  param.filter = 0.02;
  [B,A] = butter(2,param.filter);
  gps.lat = filtfilt(B,A,gps.lat);
  gps.lon = filtfilt(B,A,gps.lon);
  gps.elev = filtfilt(B,A,gps.elev);
end

return;
