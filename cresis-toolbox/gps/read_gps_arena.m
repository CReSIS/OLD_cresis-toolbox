function gps = read_gps_arena(fn, param)
% gps = read_gps_arena(fn, param)
%
% Reads in Arena GPS files that include NMEA and radar time.
%
% nmea:$GPGGA,144752.011,3857.134780,N,09515.862427,W,1,12,0.76,317.502,M,-29.504,M,,*5D
% relTimeCntr:15345172719345668
% profileCntr:2968
% ppsCntr:1534517272
%
% Input Args:
%   fn = string containing input Arena GPS filename
%     e.g. 20180817_094746_ARENA-CTU-ctu-gps.txt
%   param = tells the file GPS type and the year, month, day to determine
%     absolute time (GGA NMEA files just give the time of day)
%     .year
%     .month
%     .day
%     .time_reference: 'gps' or 'utc' (should always be 'utc')
%     .nmea_tag: NEMA string to identify good lines (e.g. '$GPGGA')
%     .clk: clock for radar_time counter, defaults to 10e6
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
%  .relTimeCntr: radar time, 64 bit counter free-running at 10 MHz
%  .profileCntr: pulse counter
%  .ppsCntr: PPS counter
%
% Example:
%   fn = 'E:\tmp\2018_Antarctica_TObas\20180817\logs\20180817_094746_ARENA-CTU-ctu-gps.txt';
%   gps = read_gps_arena(fn, struct('year',2018,'month',8,'day',17,'time_reference','utc'));
%   plot(gps.lon,gps.lat);
%   datestr(epoch_to_datenum(gps.gps_time(1)));
%   gps.utc_time = gps.gps_time - utc_leap_seconds(gps.gps_time(1))
%   plot_gps(gps)
%
% Author: John Paden
%
% See also read_gps_applanix, read_gps_atm, read_gps_csv, read_gps_litton,
%   read_gps_nmea, read_gps_novatel, read_gps_reveal, read_gps_traj, 
%   read_gps_txt, plot_gps

if ~exist('param','var') || isempty(param)
  error('Year, month, day must be specified in param struct');
end
if ~isfield(param,'clk') || isempty(param.clk)
  param.clk = 10e6;
end
  
[fid,msg] = fopen(fn,'r');
if fid < 0
  error('Error opening %s: %s', fn, msg);
end

format_str = '%s%f%f%c%f%c%u%u%f%f%c%f%c%s%s';
nmea_idx = 1;
while ~feof(fid)
  str = fgets(fid);
  [token,remain] = strtok(str,':');
  if strcmpi(token,'nmea')
    C = textscan(remain(2:end),format_str,'delimiter',', ','emptyvalue',NaN);
    [tag,UTC_time_file(nmea_idx),latitude(nmea_idx),N_S(nmea_idx),longitude(nmea_idx),E_W(nmea_idx),fix,NoSatelite,dilution,...
      altitude(nmea_idx),alt_unit,geode_ref,geode_unit,dgps,checksum] = deal(C{:});
    nmea_idx = nmea_idx + 1;
  elseif strcmpi(token,'relTimeCntr')
    relTimeCntr(nmea_idx) = str2double(remain(2:end));
  elseif strcmpi(token,'profileCntr')
    profileCntr(nmea_idx) = str2double(remain(2:end));
  elseif strcmpi(token,'ppsCntr')
    ppsCntr(nmea_idx) = str2double(remain(2:end));
  end
end
fclose(fid);

%   CONVERT LATITUDE AND LONGITUDE TO [DD.DDD] FROM [DDDMM.MMM]
lat_MM = mod(latitude,100);
lat_DD = (latitude - lat_MM)./100;
lat = lat_DD + lat_MM./60;
lon_MM = mod(longitude,100);
lon_DD = (longitude - lon_MM)./100;
lon = lon_DD + lon_MM./60;

%   IMPLY NEGATIVE LATITUDE AND LONGITUDE TO SOUTH AND WEST COORDINATES
lat(N_S == 'S') = -1 * lat(N_S == 'S');
lon(E_W == 'W') = -1 * lon(E_W == 'W');

%   CREATE NEW ELEVATION VARIABLE
elev = altitude;

% Convert HHMMSS format to ANSI-C standard, seconds since Jan 1 1970
sec = mod(UTC_time_file,100);
min = mod([UTC_time_file-sec]./100,100);
hour = [[UTC_time_file-sec]./100 - min]./100;
UTC_time = datenum_to_epoch(datenum(param.year,param.month,param.day,hour,min,sec));

% ENSURE ALL VECTORS IN 1xN FORMAT
UTC_time = reshape(UTC_time,[1 length(UTC_time)]);
lat = reshape(lat,[1 length(lat)]);
lon = reshape(lon,[1 length(lon)]);
elev = reshape(elev,[1 length(elev)]);

goodIdxs = find(~isnan(lat));
UTC_time = UTC_time(goodIdxs);
lat = lat(goodIdxs);
lon = lon(goodIdxs);
elev = elev(goodIdxs);

% ===================================================================
% Find jumps in the GPS time that are probably due to day interval
% 86400 seconds.
day_jumps = find(diff(UTC_time) < -60000);
for jump_idx = day_jumps
  UTC_time(jump_idx+1:end) = UTC_time(jump_idx+1:end) + 86400;
end

% ===================================================================
% Store outputs in structure
% ===================================================================
if strcmpi(param.time_reference,'utc')
  % UTC time stored in file, so need to add leap seconds back in
  if ~isempty(UTC_time)
    gps.gps_time = UTC_time + utc_leap_seconds(UTC_time(1));
  else
    gps.gps_time = [];
  end
else
  warning('NMEA files are usually always UTC time, but GPS time has been specified.');
  gps.gps_time = UTC_time;
end

gps.lat = lat;
gps.lon = lon;
gps.elev = elev;

gps.roll = zeros(size(gps.lat));
gps.pitch = zeros(size(gps.lat));
gps.heading = zeros(size(gps.lat));

gps.radar_time = relTimeCntr/param.clk;
gps.profileCntr = profileCntr;
gps.comp_time = ppsCntr;

% If GPS is bad, NMEA may be constant, return empty arrays in this case
if all(gps.lat(1)==gps.lat)
  warning('GPS data is constant. Not using this file.');
  gps.gps_time = [];
  gps.lat = [];
  gps.lon = [];
  gps.elev = [];
  gps.roll = [];
  gps.pitch = [];
  gps.heading = [];
  gps.radar_time = [];
  gps.profileCntr = [];
  gps.comp_time = [];
  return;
end

[comp_time_year,comp_time_month,comp_time_day] = datevec(epoch_to_datenum(gps.comp_time));
if any(comp_time_year ~= param.year)
  warning('comp_time year is %d and does not match param.year %d.', comp_time_year, param.year);
end
if any(comp_time_month ~= param.month)
  warning('comp_time month is %d and does not match param.month %d.', comp_time_month, param.month);
end
if any(comp_time_day ~= param.day)
  comp_time_day = comp_time_day(find(comp_time_day ~= param.day,1));
  warning('comp_time day is %d and does not match param.day %d.', comp_time_day, param.day);
end

return;
