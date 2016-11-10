function gps = read_gps_nmea(in_fn, param)
% gps = read_gps_nmea(in_fn, param)
%
% Reads in NMEA files and modified NMEA files that include the MCRDS
% radar time stamps.  GPS file must contain only GPGGA strings!!! If 
% non-compatible strings are present you must remove them. In Linux:
%   grep GPGGA OLD_FILENAME >NEW_FILENAME
%
% $GPGGA,120448.00,6952.4649165,N,03256.5105286,W,1,10,0.80,3135.1254,M,55.3383,M,,*7B
% $GPGGA,221343.00,7928.1873712,S,11203.3163302,W,0,10,   ,1769.932,M,     , , ,*5e
%
% Input Args:
%   in_fn = string containing input NMEA filename
%     OR cell array of strings containging input NMEA filenames
%     If the length of the cell array is one, then it operates on that one
%     file as if just a string was passed in. If the length is more than
%     one, then param.combine is set to true.
%   param = tells the file GPS type and the year, month, day to determine
%     absolute time (GGA NMEA files just give the time of day)
%     .format = scalar integer from 1 to 4, default is 1
%       1. Standard NMEA
%       2. Standard NMEA comp_time_sec comp_time_usec
%          comp_time_sec + comp_time_usec/1e6
%           = seconds since the epoch Jan 1, 1970 00:00:00
%       3. Standard NMEA comp_time_sec comp_time_usec radar_32MSB radar_32LSB
%          radar_32MSB*2^32 + radar_32LSB (64 bit MCRDS radar time, 10 MHz clock)
%       4. From reveal system
%     .year
%     .month
%     .day
%     .time_reference = 'gps' or 'utc' (should always be 'utc')
%     .combine = logical (default is false)
%        When enabled, all input files will be concatenated together
%        in the temp directory and this file will be used.
%        If in_fn is a string, then a search is done for all files in the
%        directory specified by in_fn which have the string YYYYMMDD in
%        their filename.
%        If in_fn is a cell array, then these specific files will be used.
%     .nmea_tag = NEMA string to identify good lines (e.g. '$GPGGA')
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
%  .comp_time = computer time in seconds since Jan 1, 1970 epoch (sec)
%  .radar_time = radar time, 64 bit counter free-running at 10 MHz
%
% Example:
%   fn = '/cresis/data2/MCoRDS/2010_Antarctica/GPS_new/GGA/RevealGPS_20100103A';
%   gps = read_gps_nmea(fn, struct('year',2010,'month',1,'day',3));
%   plot(gps.lon,gps.lat);
%   datestr(epoch_to_datenum(gps.gps_time(1)));
%   gps.utc_time = gps.gps_time - utc_leap_seconds(gps.gps_time(1))
%   plot_gps(gps)
%
% Author: William Blake, John Paden, Anthony Hoch, Logan Smith
%
% See also read_gps_applanix, read_gps_atm, read_gps_csv, read_gps_litton,
%   read_gps_nmea, read_gps_novatel, read_gps_reveal, read_gps_traj, 
%   read_gps_txt, plot_gps

global gRadar

debug_level = 1;

if ~exist('param','var') || isempty(param)
  error('Year, month, day must be specified in param struct');
end
if ~isfield(param,'format')
  param.format = 1;
end
if ~isfield(param,'combine')
  param.combine = 0;
end
if iscell(in_fn)
  % If in_fn is a cell array of more than one file name, then we need
  % to combined all of these files
  if length(in_fn) == 1
    in_fn = in_fn{1};
  else
    param.combine = 1;
  end
end

if param.combine
  date_str = sprintf('%04d%02d%02d',param.year,param.month,param.day);
  if iscell(in_fn)
    % in_fn is already of cell array of filenames
    in_fns = in_fn;
  else
    in_fns = get_filenames(in_fn,'',date_str,'');
  end
  tmp_nmea_path = fullfile(gRadar.tmp_path,sprintf('tmp_nmea_%s.txt',date_str));
  syscmd = sprintf('cat %s > %s',in_fns{1},tmp_nmea_path);
  system(syscmd);
  for idx=2:length(in_fns)
    syscmd = sprintf('cat %s >> %s',in_fns{idx},tmp_nmea_path);
    system(syscmd);
  end
  in_fn = tmp_nmea_path;
end

if ~exist(in_fn,'file')
  error('File does not exist %s\n', in_fn);
end

%   LOAD NMEA FILE

switch param.format
  case 1
    format_str = '%s%f%f%c%f%c%u%u%f%f%c%f%c%s%s';
  case 2
    format_str = '%s%f%f%c%f%c%u%u%f%f%c%f%c%s%s%f%f';
  case 3
    format_str = '%s%f%f%c%f%c%u%u%f%f%c%f%c%s%s%f%f%f%f';
  case 4
    format_str = '%s%f%f%c%f%c%u%u%f%f%f%f%s';
  case 5
    format_str = '%s%f%f%c%f%c%u%u%f%f%c%f%c%s%s%f%f%f%f';
end

[fid,msg] = fopen(in_fn,'r');
if fid < 0
  error('Error opening %s: %s', in_fn, msg);
end
finfo = dir(in_fn);
C_final = {};
num_lines = 0;
while ftell(fid) < finfo.bytes
  % This reader uses textscan to read in the file. If textscan finds a bad
  % line, it stops reading before the end of the file
  C = textscan(fid,format_str,'delimiter',', ','emptyvalue',NaN);
  good_lines = length(C{end});
  num_lines = num_lines + good_lines+1;
  if ftell(fid) < finfo.bytes
    readchar = fread(fid,1,'uint8');
    while readchar ~= 10 && ftell(fid) > 1
      fseek(fid,-2,0);
      readchar = fread(fid,1,'uint8');
    end
    readline = fgets(fid);
    if readline(end) ~= 10
      readline(end+1) = 10;
    end
    fprintf('Bad line %d: %s', num_lines, readline);
  end
  Clength = cellfun(@length,C);
  if good_lines > 0
    for field_idx = 1:length(C)
      if length(C_final) < field_idx
        C_final{field_idx} = C{field_idx}(1:good_lines,1);
      else
        C_final{field_idx} = cat(1,C_final{field_idx}, C{field_idx}(1:good_lines));
      end
    end
  end
end
fclose(fid);

if isempty(C_final)
  C_final = cell(1,sum(format_str=='%'));
end

switch param.format
  case 1
    [tag,UTC_time_file,latitude,N_S,longitude,E_W,fix,NoSatelite,dilution,...
      altitude,alt_unit,geode_ref,geode_unit,dgps,checksum] = deal(C_final{:});
  case 2
    [tag,UTC_time_file,latitude,N_S,longitude,E_W,fix,NoSatelite,dilution,...
      altitude,alt_unit,geode_ref,geode_unit,dgps,checksum,time1,time2] = deal(C_final{:});
  case 3
    [tag,UTC_time_file,latitude,N_S,longitude,E_W,fix,NoSatelite,dilution,...
      altitude,alt_unit,geode_ref,geode_unit,dgps,checksum,time1,time2,time3,time4] = deal(C_final{:});
  case 4
    [tag,UTC_time_file,latitude,N_S,longitude,E_W,fix,NoSatelite,dilution,...
      altitude,alt_unit,geode_ref,checksum] = deal(C_final{:});
  case 5
    [tag,UTC_time_file,latitude,N_S,longitude,E_W,fix,NoSatelite,dilution,...
      altitude,alt_unit,geode_ref,geode_unit,dgps,checksum,time1,time2,time3,time4] = deal(C_final{:});
end

if param.format == 2 || param.format == 3
  comp_time = time1 + time2/1e6;
end
if param.format == 3
  radar_time = (time3*2^32 + time4)/10e6;
end

if isfield(param,'nmea_tag')
  good_idxs = strmatch(param.nmea_tag,tag);
  tag = tag(good_idxs);
  UTC_time_file = UTC_time_file(good_idxs);
  latitude = latitude(good_idxs);
  N_S = N_S(good_idxs);
  longitude = longitude(good_idxs);
  E_W = E_W(good_idxs);
  altitude = altitude(good_idxs);
  if param.format == 2 || param.format == 3
    comp_time = comp_time(good_idxs);
  end
  if param.format == 3
    radar_time = radar_time(good_idxs);
  end
end

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
if param.format == 2 || param.format == 3
  comp_time = reshape(comp_time,[1 length(comp_time)]);
end
if param.format == 3
  radar_time = reshape(radar_time,[1 length(radar_time)]);
end


goodIdxs = find(~isnan(lat));
UTC_time = UTC_time(goodIdxs);
lat = lat(goodIdxs);
lon = lon(goodIdxs);
elev = elev(goodIdxs);
if param.format == 2 || param.format == 3
  comp_time = comp_time(goodIdxs);
end
if param.format == 3
  radar_time = radar_time(goodIdxs);
end

if param.format == 4
  % Reveal file tends to have a lot of problems... this catches some
  bad_idxs = find(abs(diff(lat)) > 0.005);
  if ~isempty(bad_idxs)
    plot(diff(lat(:)));
    fprintf('Some likely bad indices found.  Set bad_idxs in the code and run.\n')
    keyboard
    % For example, two bad ranges: bad_idxs = [12690:12768, 16608:16643];
    bad_idxs = [];  % SET THIS TO APPROPRIATE VALUE EACH TIME
    good_idxs = setdiff(1:length(lat), bad_idxs);
    plot(diff(lat(good_idxs)),'.');
    keyboard
    UTC_time = UTC_time(good_idxs);
    lat = lat(good_idxs);
    lon = lon(good_idxs);
    elev = elev(good_idxs);    
  end
end

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

if param.format == 2 || param.format == 3
  gps.comp_time = comp_time;
end
if param.format == 3
  gps.radar_time = radar_time;
end

return;
