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

gps_path = fullfile(support_path,'gps','2002_Greenland_P3');
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

in_base_path = fullfile(data_support_path,'2002_Greenland_P3');
gps_path = fullfile(support_path,'gps','2002_Greenland_P3');

file_idx = 0; in_fns = {}; out_fns = {}; file_type = {}; params = {};

% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'20020518_nmea.csv');
% out_fns{file_idx} = 'gps_20020518.mat';
% file_type{file_idx} = 'csv';
% params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
% gps_source{file_idx} = 'atm-final_20020518'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'20020520_nmea.csv');
out_fns{file_idx} = 'gps_20020520.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_20020520'; 

file_idx = file_idx + 1;
in_fns{file_idx} = fullfile(in_base_path,'20020522_nmea.csv');
out_fns{file_idx} = 'gps_20020522.mat';
file_type{file_idx} = 'csv';
params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
gps_source{file_idx} = 'atm-final_20020522'; 
% 
% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'20020524_nmea.csv');
% out_fns{file_idx} = 'gps_20020524.mat';
% file_type{file_idx} = 'csv';
% params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
% gps_source{file_idx} = 'atm-final_20020524'; 
% 
% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'20020528_nmea.csv');
% out_fns{file_idx} = 'gps_20020528.mat';
% file_type{file_idx} = 'csv';
% params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
% gps_source{file_idx} = 'atm-final_20020528'; 
% 
% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'20020529_nmea.csv');
% out_fns{file_idx} = 'gps_20020529.mat';
% file_type{file_idx} = 'csv';
% params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
% gps_source{file_idx} = 'atm-final_20020529';
% 
% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'20020530_nmea.csv');
% out_fns{file_idx} = 'gps_20020530.mat';
% file_type{file_idx} = 'csv';
% params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
% gps_source{file_idx} = 'atm-final_20020530';
% 
% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'20020531_nmea.csv');
% out_fns{file_idx} = 'gps_20020531.mat';
% file_type{file_idx} = 'csv';
% params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
% gps_source{file_idx} = 'atm-final_20020531';
% 
% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'20020601_nmea.csv');
% out_fns{file_idx} = 'gps_20020601.mat';
% file_type{file_idx} = 'csv';
% params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
% gps_source{file_idx} = 'atm-final_20020601';
% 
% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'20020604_nmea.csv');
% out_fns{file_idx} = 'gps_20020604.mat';
% file_type{file_idx} = 'csv';
% params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
% gps_source{file_idx} = 'atm-final_20020604';

make_gps;

