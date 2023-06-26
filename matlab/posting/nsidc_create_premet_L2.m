function nsidc_create_premet_L2(in_fn,out_fn,param)
% nsidc_create_premet_L2(in_fn,out_fn,param)
%
% Generate a .premet file from an existing .csv data file. Same filename, 
% different extension.
% Called from nsidc_delivery_script.m
%
% in_fn = string containing the full path of given .csv file
% out_fn = string containing the full path of the output .premet file
% param = structure containing premet field information
%
% This function aims to extract basic platform information from a L2 data in 
% .csv file to a .premet file, including filename, version ID, begin and
% end date, begin and end time, theme ID, aircraft ID, platform
% shortname, and so on.
%
% Author: Yi Zhu
%
% See also: type "nsidc_help.m"

%% Load the data from .csv file to a 1x9 cell array
fid_csv = fopen(in_fn);
csv_data = textscan(fid_csv,'%f%f%f%f%f%s%f%f%f','delimiter',',','headerlines',1);
fclose(fid_csv);

%% Write data into .premet file
out_fn_dir = fileparts(out_fn);
if ~exist(out_fn_dir,'dir')
  mkdir(out_fn_dir);
end
[fid,msg] = fopen(out_fn,'w+');
if fid < 1
  fprintf('Could not open file %s\n', out_fn);
  error(msg);
end

% Data_FileName
fprintf(fid, 'Data_FileName=%s\n',param.data_fn_name);

% VersionID_local
fprintf(fid, 'VersionID_local=%s\n',param.version_id);

year = str2double(csv_data{6}{1}(1:4));
month = str2double(csv_data{6}{1}(5:6));
day = str2double(csv_data{6}{1}(7:8));
gps_time(1) = datenum(year,month,day,0,0,csv_data{3}(1));
year = str2double(csv_data{6}{end}(1:4));
month = str2double(csv_data{6}{end}(5:6));
day = str2double(csv_data{6}{end}(7:8));
gps_time(2) = datenum(year,month,day,0,0,csv_data{3}(end));

% Begin_date
fprintf(fid, 'Begin_date=%s\n', datestr(gps_time(1),'yyyy-mm-dd'));
% End_date
fprintf(fid, 'End_date=%s\n', datestr(gps_time(2),'yyyy-mm-dd'));

% Begin_time
fprintf(fid, 'Begin_time=%s\n', datestr(gps_time(1),'HH:MM:SS.FFF'));
% End_time
fprintf(fid, 'End_time=%s\n', datestr(gps_time(2),'HH:MM:SS.FFF'));

% Manually input platform information
fprintf(fid, 'Container=AdditionalAttributes\n');
fprintf(fid, 'AdditionalAttributeName=ThemeID\n');
fprintf(fid, 'ParameterValue=%s\n', param.nsidc_season_name);
fprintf(fid, 'Container=AdditionalAttributes\n');
fprintf(fid, 'AdditionalAttributeName=AircraftID\n');
fprintf(fid, 'ParameterValue=%s\n', param.nsidc_aircraft_id);
fprintf(fid, 'Container=AssociatedPlatformInstrumentSensor\n');
fprintf(fid, 'AssociatedPlatformShortName=%s\n', param.nsidc_platform_short_name);
fprintf(fid, 'AssociatedInstrumentShortName=%s\n', param.nsidc_instrument_short_name);
fprintf(fid, 'AssociatedSensorShortName=%s', param.nsidc_sensor_short_name);
  
fclose(fid);

end