function gps = read_gps_general_ascii(fn,param)
% function gps = read_gps_general_ascii(fn,param)
%
% Generalized ASCII file reader (CSV, tab delimited, etc)
%
% Example 1: CSV File with the following content
%
% Year, DayOfYear, SecondsOfDay(UTC), Pitch(deg), Roll(deg), Heading(deg)
% 2009, 117, 39196.000,   2.587,   0.055,  68.313
%
% fn = '/cresis/snfs1/dataproducts/metadata/2009_Greenland_P3/20090427_ATM_vldInsExtract.txt';
% param = [];
% param.format_str = '%f%f%f%f%f%f';
% param.types = {'year','day','sec','pitch_deg','roll_deg','heading_deg'};
% param.textscan = {'delimiter',',','emptyvalue',NaN};
% param.headerlines = 1;
% param.time_reference = 'utc';
% gps = read_gps_general_ascii(fn,param);
%
% Example 2: CSV File with the following content
%
% NAME          , -- LATITUDE --, -- LONGITUDE --,   HEIGHT,Q,StDev,-- VE --,-- VN --,-- VZ -- (velocity m/s)
% 000327074.000 , 67 05 35.91995, -50 17 04.50389,  245.120,5,0.506,  -0.038,  -0.022,   0.059,

% fn = 'P:\HF_Sounder\Greenland2016\ProcessedGPS\04132016_UAV.csv';
% param = [];
% param.format_str = '%f%s%s%f%f%f%f%f%f';
% param.types = {'sow','lat_dsm','lon_dsm','elev_m'};
% param.textscan = {'delimiter',',','emptyvalue',NaN};
% param.headerlines = 1;
% param.time_reference = 'gps';
% param.year = 2016;
% param.month = 4;
% param.day = 13;
% gps = read_gps_general_ascii(fn,param);
%
% Author: John Paden

[fid,msg] = fopen(fn,'r');
if fid < 0
  error('Error opening %s: %s', fn, msg);
end

% Skip header lines
for tmp=1:param.headerlines
  fgets(fid);
end

% Read in lines of file
finfo = dir(fn);
C_final = {};
while ftell(fid) < finfo.bytes
  % This reader uses textscan to read in the file. If textscan finds a bad
  % line, it stops reading before the end of the file
  C = textscan(fid,param.format_str,param.textscan{:});
  if ftell(fid) < finfo.bytes
    good_lines = length(C{end});
    fprintf('Bad line %d: %s\n', good_lines+1, fgets(fid));
  else
    good_lines = length(C{end});
  end
  for field_idx = 1:length(C)
    if length(C_final) < field_idx
      C_final{field_idx} = C{field_idx}(1:good_lines);
    else
      C_final{field_idx} = cat(1,C_final{field_idx}, C{field_idx}(1:good_lines));
    end
  end
end
fclose(fid);

% Convert from cell to struct/fieldnames for better readability below
tmp_gps = [];
for idx = 1:min(length(C_final),length(param.types))
  tmp_gps.(param.types{idx}) = reshape(C_final{idx},[1 length(C_final{idx})]);
end

% Interpret fields that were read in and create the output "gps" struct
num_rows = -1;
gps = [];

% Create gps time field
if isfield(tmp_gps,'year')
  year = tmp_gps.year;
elseif isfield(param,'year')
  year = param.year;
else
  year = 0;
end
if isfield(tmp_gps,'month')
  month = tmp_gps.month;
elseif isfield(param,'month')
  month = param.month;
else
  month = 0;
end
if isfield(tmp_gps,'day')
  day = tmp_gps.day;
elseif isfield(param,'day')
  day = param.day;
else
  day = 0;
end
if isfield(tmp_gps,'hour')
  hour = tmp_gps.hour;
else
  hour = 0;
end
if isfield(tmp_gps,'minute')
  minute = tmp_gps.minute;
else
  minute = 0;
end
if isfield(tmp_gps,'sec')
  sec = tmp_gps.sec;
else
  sec = 0;
end
gps.gps_time = datenum_to_epoch(datenum(year,month,day,hour,minute,sec));

if isfield(tmp_gps,'sow')
  if num_rows == -1
    num_rows = length(tmp_gps.sow);
  elseif length(tmp_gps.sow) < num_rows
    tmp_gps.sow(end+1:num_rows) = NaN;
  elseif length(tmp_gps.sow) < num_rows
    tmp_gps.sow = tmp_gps.sow(1:num_rows);
  end
  gps.gps_time = gps_sow_to_epoch(tmp_gps.sow,param);
end

% Look for trajectory fields
if isfield(tmp_gps,'lat_deg')
  if num_rows == -1
    num_rows = length(tmp_gps.lat_deg);
  elseif length(tmp_gps.lat_deg) < num_rows
    tmp_gps.lat_deg(end+1:num_rows) = NaN;
  elseif length(tmp_gps.lat_deg) < num_rows
    tmp_gps.lat_deg = tmp_gps.lat_deg(1:num_rows);
  end
  gps.lat = tmp_gps.lat_deg;
end
if isfield(tmp_gps,'lat_dsm')
  if num_rows == -1
    num_rows = length(tmp_gps.lat_dsm);
  elseif length(tmp_gps.lat_dsm) < num_rows
    tmp_gps.lat_dsm(end+1:num_rows) = deal(NaN);
  elseif length(tmp_gps.lat_dsm) < num_rows
    tmp_gps.lat_dsm = tmp_gps.lat_dsm(1:num_rows);
  end
  gps.lat = zeros(1,num_rows);
  for idx = 1:num_rows
    tmp = sscanf(tmp_gps.lat_dsm{idx},'%f');
    gps.lat(idx) = tmp(1) + sign(tmp(1))*(tmp(2)/60 + tmp(3)/3600);
  end
end
if isfield(tmp_gps,'lon_deg')
  if num_rows == -1
    num_rows = length(tmp_gps.lon_deg);
  elseif length(tmp_gps.lon_deg) < num_rows
    tmp_gps.lon_deg(end+1:num_rows) = NaN;
  elseif length(tmp_gps.lon_deg) < num_rows
    tmp_gps.lon_deg = tmp_gps.lon_deg(1:num_rows);
  end
  gps.lon = tmp_gps.lon_deg;
end
if isfield(tmp_gps,'lon_dsm')
  if num_rows == -1
    num_rows = length(tmp_gps.lon_dsm);
  elseif length(tmp_gps.lon_dsm) < num_rows
    tmp_gps.lon_dsm(end+1:num_rows) = deal(NaN);
  elseif length(tmp_gps.lon_dsm) < num_rows
    tmp_gps.lon_dsm = tmp_gps.lon_dsm(1:num_rows);
  end
  gps.lon = zeros(1,num_rows);
  for idx = 1:num_rows
    tmp = sscanf(tmp_gps.lon_dsm{idx},'%f');
    gps.lon(idx) = tmp(1) + sign(tmp(1))*(tmp(2)/60 + tmp(3)/3600);
  end
end
if isfield(tmp_gps,'elev_m')
  if num_rows == -1
    num_rows = length(tmp_gps.elev_m);
  elseif length(tmp_gps.elev_m) < num_rows
    tmp_gps.elev_m(end+1:num_rows) = NaN;
  elseif length(tmp_gps.elev_m) < num_rows
    tmp_gps.elev_m = tmp_gps.elev_m(1:num_rows);
  end
  gps.elev = tmp_gps.elev_m;
end

% Look for attitude fields
if isfield(tmp_gps,'roll_deg')
  if num_rows == -1
    num_rows = length(tmp_gps.roll_deg);
  elseif length(tmp_gps.roll_deg) < num_rows
    tmp_gps.roll_deg(end+1:num_rows) = NaN;
  elseif length(tmp_gps.roll_deg) < num_rows
    tmp_gps.roll_deg = tmp_gps.roll_deg(1:num_rows);
  end
  gps.roll = tmp_gps.roll_deg/180*pi;
end
if isfield(tmp_gps,'pitch_deg')
  if num_rows == -1
    num_rows = length(tmp_gps.pitch_deg);
  elseif length(tmp_gps.pitch_deg) < num_rows
    tmp_gps.pitch_deg(end+1:num_rows) = NaN;
  elseif length(tmp_gps.pitch_deg) < num_rows
    tmp_gps.pitch_deg = tmp_gps.pitch_deg(1:num_rows);
  end
  gps.pitch = tmp_gps.pitch_deg/180*pi;
end
if isfield(tmp_gps,'heading_deg')
  if num_rows == -1
    num_rows = length(tmp_gps.heading_deg);
  elseif length(tmp_gps.heading_deg) < num_rows
    tmp_gps.heading_deg(end+1:num_rows) = NaN;
  elseif length(tmp_gps.heading_deg) < num_rows
    tmp_gps.heading_deg = tmp_gps.heading_deg(1:num_rows);
  end
  gps.heading = tmp_gps.heading_deg/180*pi;
end

% Final clean up to make sure all fields are present
if isfield(gps,'gps_time')
  if strcmpi(param.time_reference,'utc')
    % UTC time stored in file, so need to add leap seconds back in
    gps.gps_time = gps.gps_time + utc_leap_seconds(gps.gps_time(1));
  end
else
  gps.gps_time = NaN*zeros(1,num_rows);
end

if ~isfield(gps,'lat')
  gps.lat = NaN*zeros(1,num_rows);
end

if ~isfield(gps,'lon')
  gps.lon = NaN*zeros(1,num_rows);
end

if ~isfield(gps,'elev')
  gps.elev = NaN*zeros(1,num_rows);
end

if ~isfield(gps,'roll')
  gps.roll = zeros(1,num_rows);
end

if ~isfield(gps,'pitch')
  gps.pitch = zeros(1,num_rows);
end

if ~isfield(gps,'heading')
  gps.heading = zeros(1,num_rows);
end

return;
