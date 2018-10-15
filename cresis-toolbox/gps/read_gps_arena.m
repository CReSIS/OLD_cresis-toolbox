function gps = read_gps_arena(fn, param)
% gps = read_gps_arena(fn, param)
%
% Reads in Arena GPS files that include NMEA and radar time.
%
% Ignores NMEA strings besides GPRMC and GPGGA which show up in the file as
% the following. GPGGA is required. GPRMC is optional and is used to
% estimate the date and heading.
%
% relTimeCntr:15393734057677147
% profileCntr:7678
% ppsCntr:1539373405
% nmea:$GPRMC,194325.010,A,3857.135037,N,09515.861757,W,0.000,0.00,121018,,,A*46
% relTimeCntr:15393734057677147
% profileCntr:7678
% ppsCntr:1539373405
% nmea:$GPGGA,194325.010,3857.135037,N,09515.861757,W,1,9,0.88,311.412,M,-29.504,M,,*65
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

UTC_time_file = [];
latitude = [];
N_S = [];
longitude = [];
E_W = [];
altitude = [];
relTimeCntr = [];
profileCntr = [];
ppsCntr = [];
gps_date = [];
heading = [];

format_str_GPRMC = '%s%f%s%f%c%f%c%f%f%f%f%s%s';
format_str_GPGGA = '%s%f%f%c%f%c%u%u%f%f%c%f%c%s%s';
nmea_idx = 1;
relTimeCntrTmp = NaN;
profileCntrTmp = NaN;
ppsCntrTmp = NaN;
gps_date_tmp = NaN;
heading_tmp = NaN;
while ~feof(fid)
  str = fgets(fid);
  [token,remain] = strtok(str,':');
  if strcmpi(token,'nmea')
    if numel(str)>=11 && strcmp(str(7:11),'GPGGA')
      C = textscan(remain(2:end),format_str_GPGGA,'delimiter',', ','emptyvalue',NaN);
      [tag,UTC_time_file_tmp,latitude_tmp,N_S_tmp,longitude_tmp,E_W_tmp,fix,NoSatelite,dilution,...
        altitude_tmp,alt_unit,geode_ref,geode_unit,dgps,checksum] = deal(C{:});
      
      if ~isnan(relTimeCntrTmp) && ~isnan(profileCntrTmp) && ~isnan(ppsCntrTmp)
        UTC_time_file(nmea_idx) = UTC_time_file_tmp;
        latitude(nmea_idx) = latitude_tmp;
        N_S(nmea_idx) = N_S_tmp;
        longitude(nmea_idx) = longitude_tmp;
        E_W(nmea_idx) = E_W_tmp;
        altitude(nmea_idx) = altitude_tmp;
        relTimeCntr(nmea_idx) = relTimeCntrTmp;
        profileCntr(nmea_idx) = profileCntrTmp;
        ppsCntr(nmea_idx) = ppsCntrTmp;
        gps_date(nmea_idx) = gps_date_tmp;
        heading(nmea_idx) = heading_tmp;
        nmea_idx = nmea_idx + 1;
        relTimeCntrTmp = NaN;
        profileCntrTmp = NaN;
        ppsCntrTmp = NaN;
        gps_date_tmp = NaN;
        heading_tmp = NaN;
      end
    elseif numel(str)>=11 && strcmp(str(7:11),'GPRMC')
      C = textscan(remain(2:end),format_str_GPRMC,'delimiter',', ','emptyvalue',NaN);
      [tag,UTC_time_file_tmp,nav_rx_warning,latitude_tmp,N_S_tmp,longitude_tmp,E_W_tmp,speed,heading_tmp,...
        gps_date_tmp,mag_Var,mag_var_E_W,checksum] = deal(C{:});
    end
  elseif strcmpi(token,'relTimeCntr')
    relTimeCntrTmp = str2double(remain(2:end));
  elseif strcmpi(token,'profileCntr')
    profileCntrTmp = str2double(remain(2:end));
  elseif strcmpi(token,'ppsCntr')
    ppsCntrTmp = str2double(remain(2:end));
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
  % gps_date = HHMMSS format
sec = mod(UTC_time_file,100);
min = mod((UTC_time_file-sec)./100,100);
hour = ((UTC_time_file-sec)./100 - min)./100;
if ~all(isnan(gps_date))
  % gps_date = DDMMYY format
  gps_date = interp_finite(gps_date);
  year = mod(gps_date,100);
  month = mod((gps_date-year)./100,100);
  day = ((gps_date-year)./100 - month)./100;
  year = year+2000; % Add in century information
  if isfield(param,'year')
    year(:) = param.year(:);
  else
    param.year = year;
  end
  if isfield(param,'month')
    month(:) = param.month;
  else
    param.month = month;
  end
  if isfield(param,'day')
    day(:) = param.day;
  else
    param.day = day;
  end
  UTC_time = datenum_to_epoch(datenum(year,month,day,hour,min,sec));
else
  UTC_time = datenum_to_epoch(datenum(param.year,param.month,param.day,hour,min,sec));
end

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

if all(isnan(gps_date))
  % ===================================================================
  % Find jumps in the GPS time that are probably due to day interval
  % 86400 seconds.
  day_jumps = find(diff(UTC_time) < -60000);
  for jump_idx = day_jumps
    UTC_time(jump_idx+1:end) = UTC_time(jump_idx+1:end) + 86400;
  end
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
gps.heading = interp_finite(heading,0);

gps.radar_time = relTimeCntr/param.clk;
gps.profileCntr = profileCntr;
gps.comp_time = ppsCntr;

% If GPS is bad, NMEA may be constant, return empty arrays in this case
if all(gps.lat(1)==gps.lat)
  warning('GPS data is constant. This generally means invalid data so not loading this file.');
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

if all(gps.radar_time(1)==gps.radar_time)
  warning('gps.radar_time is constant. GPS 1 PPS may not have been received. Radar synchronization to GPS will not be possible with this file.');
end

if all(gps.comp_time(1)==gps.comp_time)
  warning('gps.comp_time is constant. GPS 1 PPS may not have been received.');
  
else
  [comp_time_year,comp_time_month,comp_time_day] = datevec(epoch_to_datenum(gps.comp_time));
  if any(comp_time_year ~= param.year)
    warning('comp_time year values are %s and does not always match values param.year %s.', ...
      mat2str_generic(unique(comp_time_year)), mat2str_generic(unique(param.year)));
  end
  if any(comp_time_month ~= param.month)
    warning('comp_time month values are %s and does not always match values param.month %s.', ...
      mat2str_generic(unique(comp_time_month)), mat2str_generic(unique(param.month)));
  end
  if any(comp_time_day ~= param.day)
    warning('comp_time day values are %s and does not always match values param.day %s.', ...
      mat2str_generic(unique(comp_time_day)), mat2str_generic(unique(param.day)));
  end
end

