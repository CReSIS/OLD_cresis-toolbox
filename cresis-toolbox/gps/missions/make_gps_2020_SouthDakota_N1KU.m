% script make_gps_2020_SouthDakota_N1KU
%
% Makes the GPS files for 2020_South_Dakota_N1KU field season

tic;

%% Set Season Information
% ======================================================================
season_name = '2020_SouthDakota_N1KU';

%% Construct default paths/make output directory
% ======================================================================

global gRadar;

support_path = '';
data_support_path = '';

if isempty(support_path)
  support_path = gRadar.support_path;
end

gps_path = fullfile(support_path,'gps',season_name);
if ~exist(gps_path,'dir')
  fprintf('Making directory %s\n', gps_path);
  fprintf('  Press a key to proceed\n');
  pause;
  mkdir(gps_path);
end

if isempty(data_support_path)
  data_support_path = gRadar.data_support_path;
end

%% Set GPS File Information
% ======================================================================
debug_level = 1;

in_base_path = fullfile(data_support_path,season_name);

file_idx = 0; in_fns = {}; out_fns = {}; file_type = {}; params = {}; gps_source = {};
sync_fns = {}; sync_params = {};

gps_source_to_use = 'NMEA';
% gps_source_to_use = 'novatel';

if strcmpi(gps_source_to_use,'NMEA')
  %% NMEA
  % ======================================================================
  
%   year = 2021; month = 01; day = 30;
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = get_filenames(fullfile(in_base_path,sprintf('%04d%02d%02d',year,month,day)),'GPS','','.txt');
%   out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat', year, month, day);
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',year,'month',month,'day',day,'format',1,'time_reference','utc');
%   gps_source{file_idx} = 'nmea-field';
%   sync_flag{file_idx} = 0;
%   date_str{file_idx} = sprintf('%04d%02d%02d',year,month,day);
  
  year = 2021; month = 02; day = 01;
  file_idx = file_idx + 1;
  in_fns{file_idx} = get_filenames(fullfile(in_base_path,sprintf('%04d%02d%02d',year,month,day)),'GPS','','.txt');
  out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat', year, month, day);
  file_type{file_idx} = 'NMEA';
  params{file_idx} = struct('year',year,'month',month,'day',day,'format',1,'time_reference','utc');
  gps_source{file_idx} = 'nmea-field';
  sync_flag{file_idx} = 0;
  date_str{file_idx} = sprintf('%04d%02d%02d',year,month,day);
  
elseif strcmpi(gps_source_to_use,'novatel')
  %% NOVATEL
  % ======================================================================
  
  year = 2020; month = 1; day = 28;
  file_idx = file_idx + 1;
  in_fns{file_idx} = get_filenames(in_base_path,'',sprintf('%04d%02d%02d',year,month,day),'ver2.txt');
  out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat', year, month, day);
  file_type{file_idx} = 'General_ASCII';
  params{file_idx} = struct('time_reference','gps','headerlines',21,'format_str','%s%s%f%f%f%f%f%f%f%f%f');
  params{file_idx}.types = {'date_MDY','time_HMS','lat_deg','lon_deg','elev_m','roll_deg','pitch_deg','heading_deg','tmp_1','tmp_2','tmp_3'};
  params{file_idx}.textscan = {};
  gps_source{file_idx} = 'cresis-final20200601';
  sync_flag{file_idx} = 0;
  date_str{file_idx} = sprintf('%04d%02d%02d',year,month,day);
  
end

%% Make GPS files
% ======================================================================
gps_make;

%% Manual fixes to GPS files
% ======================================================================
for idx = 1:length(file_type)
  out_fn = fullfile(gps_path,out_fns{idx});
  
  gps = load(out_fn);
  if regexpi(gps.gps_source,'nmea')
    warning('Making monotonic gps time: %s', out_fn);
    [gps,error_flag] = gps_make_monotonic(gps);
    
    warning('Smoothing elevation and heading data: %s', out_fn);
    gps.elev = sgolayfilt(gps.elev,2,101); % Adjust filter length as needed to remove high frequency noise
    heading_x = cos(gps.heading);
    heading_y = sin(gps.heading);
    heading_x  = sgolayfilt(heading_x,2,101); % Adjust filter length as needed to remove high frequency noise
    heading_y  = sgolayfilt(heading_y,2,101); % Adjust filter length as needed to remove high frequency noise
    gps.heading = atan2(heading_y,heading_x);
    
    save(out_fn,'-append','-struct','gps','elev','heading');
  end
end
