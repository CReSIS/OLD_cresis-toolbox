% script make_gps_2019_southdakota_CESSNA
%
% Makes the GPS files for 2018 South Dakota CESSNA field season

tic;

global gRadar;

support_path = '';
data_support_path = '';

if isempty(support_path)
  support_path = gRadar.support_path;
end

gps_path = fullfile(support_path,'gps','2019_SouthDakota_CESSNA');
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

in_base_path = fullfile(data_support_path,'2019_SouthDakota_CESSNA');

file_idx = 0; in_fns = {}; out_fns = {}; file_type = {}; params = {}; gps_source = {};
sync_fns = {}; sync_params = {};

gps_source_to_use = 'NMEA';
if strcmpi(gps_source_to_use,'NMEA')
%     year = 2019; month = 12; day = 11 ;
%     file_idx = file_idx + 1;
%     in_fns{file_idx} = get_filenames(fullfile(in_base_path,sprintf('%04d%02d%02d',year,month,day)),'nmea','','.gps');
%     out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat', year, month, day);
%     file_type{file_idx} = 'NMEA';
%     params{file_idx} = struct('year',year,'month',month,'day',day,'format',3,'time_reference','utc');
%     gps_source{file_idx} = 'nmea-field';
%     sync_flag{file_idx} = 0;
    
    year = 2020; month = 01; day = 16;
    file_idx = file_idx + 1;
    in_fns{file_idx} = get_filenames(fullfile(in_base_path,sprintf('%04d%02d%02d',year,month,day)),'nmea','','.gps');
    out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat', year, month, day);
    file_type{file_idx} = 'NMEA';
    params{file_idx} = struct('year',year,'month',month,'day',day,'format',3,'time_reference','utc');
    gps_source{file_idx} = 'nmea-field';
    sync_flag{file_idx} = 0;
end

% ======================================================================
% Read and translate files according to user settings
% ======================================================================
make_gps;