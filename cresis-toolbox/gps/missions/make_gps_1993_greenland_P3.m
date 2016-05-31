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

in_base_path = fullfile(data_support_path,'1993_Greenland_P3');
gps_path = fullfile(support_path,'gps','1993_Greenland_P3');

file_idx = 0; in_fns = {}; out_fns = {}; file_type = {}; params = {};


file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19930623_nmea.csv');
out_fns{file_idx} = 'gps_19930623.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19930623'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19930624_nmea.csv');
out_fns{file_idx} = 'gps_19930624.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19930624'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19930627_nmea.csv');
out_fns{file_idx} = 'gps_19930627.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19930627'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19930628_nmea.csv');
out_fns{file_idx} = 'gps_19930628.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19930628'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19930701_nmea.csv');
out_fns{file_idx} = 'gps_19930701.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19950523'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19930702_nmea.csv');
out_fns{file_idx} = 'gps_19930702.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19930702'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19930703_nmea.csv');
out_fns{file_idx} = 'gps_19930703.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19930703'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19930707_nmea.csv');
out_fns{file_idx} = 'gps_19930707.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19930707'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19930708_nmea.csv');
out_fns{file_idx} = 'gps_19930708.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19930708'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'19930709_nmea.csv');
out_fns{file_idx} = 'gps_19930709.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_19930709'; 

make_gps;

