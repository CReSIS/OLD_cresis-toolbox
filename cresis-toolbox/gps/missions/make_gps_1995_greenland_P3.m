% script make_gps_2002_greenland_P3
% Makes the DGPSwINS????? files for 2002 Greenland P3 field season
%see icards_gps_missinNASA_csv.m to get csv files for days without
%trajectory data. (check time reference: should be gps)
tic;

global gRadar;

support_path = '';
data_support_path = '';

if isempty(support_path)
  support_path = gRadar.support_path;
end

gps_path = fullfile(support_path,'gps','1995_Greenland_P3');
if ~exist(gps_path,'dir')
  fprintf('Making directory %s\n', gps_path);
  fprintf('  Press a key to proceed\n');
  pause;
  mkdir(gps_path);
end

if isempty(data_support_path)
  data_support_path = gRadar.data_support_path;
end

% ======================================================================
% User Settings
% ======================================================================
debug_level = 1;

in_base_path = fullfile(data_support_path,'1995_Greenland_P3');
gps_path = fullfile(support_path,'gps','1995_Greenland_P3');

file_idx = 0; in_fns = {}; out_fns = {}; file_type = {}; params = {};


file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19950518_nmea.csv');
out_fns{file_idx} = 'gps_19950518.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19950518'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19950519_nmea.csv');
out_fns{file_idx} = 'gps_19950519.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19950519'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19950520_nmea.csv');
out_fns{file_idx} = 'gps_19950520.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19950520'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19950522_nmea.csv');
out_fns{file_idx} = 'gps_19950522.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19950522'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19950523_nmea.csv');
out_fns{file_idx} = 'gps_19950523.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19950523'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19950524_nmea.csv');
out_fns{file_idx} = 'gps_19950524.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19950524'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19950526_nmea.csv');
out_fns{file_idx} = 'gps_19950526.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19950526'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19950527_nmea.csv');
out_fns{file_idx} = 'gps_19950527.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19950527'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19950530_nmea.csv');
out_fns{file_idx} = 'gps_19950530.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19950530'; 

make_gps;

