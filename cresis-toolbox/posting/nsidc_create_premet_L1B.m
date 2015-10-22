function nsidc_create_premet_L1B(in_fn,out_fn,param)
% nsidc_create_premet_L1B(in_fn,out_fn,param)
%
% Generate a .premet file from an existing .mat data file. Same filename, 
% different extension.
% Called from nsidc_delivery_script.m
%
% file_path = the full path of given .mat file
%
% This function aims to extract basic platform information from a L1B data in 
% .mat file to a .premet file, including filename, version ID, begin and
% end date, begin and end time, theme ID, aircraft ID, platform
% shortname, and so on.
%
% Author: Yi Zhu
%
% See also: type "nsidc_help.m"

%% Generate the .spatial file. Obtain filename from original .mat data file
metaData = load(in_fn);

%% Write data into .spatial file
out_fn_dir = fileparts(out_fn);
if ~exist(out_fn_dir,'dir')
  mkdir(out_fn_dir);
end

[fid,msg] = fopen(out_fn,'w+');

if fid < 1
  fprintf('Could not open file %s\n', out_fn);
  error(msg);
end

%% Write data into .premet file
% Data_FileName
[~,Filename,~] = fileparts(out_fn);
fprintf(fid, 'Data_FileName=%s\n',Filename);

% VersionID_local
fprintf(fid, 'VersionID_local=%s\n',param.version_id);

% Handle with GPS_time
begin_date = epoch_to_datenum(metaData.GPS_time(1));
beginDateDetail = datestr(begin_date,'yyyy-mm-dd HH:MM:SS.FFF');
yymmdd_begin = beginDateDetail(1:10);
time_begin = beginDateDetail(12:end);

end_date = epoch_to_datenum(metaData.GPS_time(end));
endDateDetail = datestr(end_date,'yyyy-mm-dd HH:MM:SS.FFF');
yymmdd_end = endDateDetail(1:10);
time_end = endDateDetail(12:end);

% Begin_date
fprintf(fid, 'Begin_date=%s\n', yymmdd_begin);
% End_date
fprintf(fid, 'End_date=%s\n', yymmdd_end);

% Begin_time
fprintf(fid, 'Begin_time=%s\n', time_begin);
% End_time
fprintf(fid, 'End_time=%s\n', time_end);

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
