function nsidc_create_premet_0(in_fn,out_fn,param)
% nsidc_create_premet_L0(in_fn,out_fn,param)
%
% Generate a .premet file from an existing raw data .dat or .bin file. 
% Called from nsidc_delivery_rawdata_script.m
%
% file_path = the full path of given .dat or .bn file
%
% This function aims to extract basic platform information from a raw data in 
% .dat or .bin file to a .premet file, including filename, version ID, begin and
% end date, begin and end time, theme ID, aircraft ID, platform
% shortname, and so on.
%
% Author: Yi Zhu
%
% See also: type "nsidc_help.m"

%% Generate the .premet file. Obtain date from output file name


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

%% Write data into .premet file
% Data_FileName
[~,Filename,~] = fileparts(out_fn);
fprintf(fid, 'Data_FileName=%s\n',Filename);

% VersionID_local
fprintf(fid, 'VersionID_local=%s\n',param.version_id);

% Handle with date
yymmdd_begin = sprintf('%s-%s-%s',Filename(14:17),Filename(18:19),Filename(20:21));
time_begin = '00:00:00.000';

yymmdd_end = yymmdd_begin;
time_end = '23:59:59.999';

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
