function data = read_seaice_kurtz(fn)
% data = read_seaice_kurtz(fn)
%
% Example:
% fn = 'c:\tmp\snow_thickness\2014_GR_NASA\OIB_20140326_IDCSI2.txt';
% data = read_seaice_kurtz(fn);
% plot(data.snow_depth); % Snow Thickness
% plot(data.n_atm); % Number of ATM elevations used
% plot(data.elev); % ATM elevation
% plot(data.ssh); % Sea surface height
% plot(data.surface_roughness); % Standard deviation of ATM elevations
% plot(data.sa_int_elev); % Height of snow air interface
% plot(data.si_int_elev); % Height of snow ice interface
%
% Author: John Paden

fid = fopen(fn);
header_strings = fgets(fid);
header_strings = textscan(header_strings,'%s','Delimiter',',');
header_strings = header_strings{1};
header_format = '';
for idx = 1:length(header_strings)
  if strcmpi(header_strings{idx},'ATM_File_Name')
    header_format = cat(2,header_format,'%s');
  else
    header_format = cat(2,header_format,'%f');
  end
end
data_vals = textscan(fid,header_format,'Delimiter',',');
data = [];
for idx = 1:length(header_strings)
  if ~strcmpi(header_strings{idx},'ATM_File_Name')
    data_vals{idx}(data_vals{idx} == -99999) = NaN;
  end
  data.(header_strings{idx}) = data_vals{idx};
end
fclose(fid);

data.gps_time = datenum(floor(data.date/10000), ...
  mod(floor(data.date/100),100), ...
  mod(data.date,100), 0, 0, data.elapsed);
data.gps_time = datenum_to_epoch(data.gps_time);
data.gps_time = data.gps_time + utc_leap_seconds(data.gps_time(1));

return

