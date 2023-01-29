% script gps_create_2022_Antarctica_GroundGHOST
%
% Makes the GPS files for 2022_Antarctica_GroundGHOST field season

%% Setup
% =========================================================================
tic;

global gRadar;

support_path = '';
data_support_path = '';

if isempty(support_path)
  support_path = gRadar.support_path;
end

season_name = '2022_Antarctica_GroundGHOST';

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

in_base_path = fullfile(data_support_path,season_name);

file_idx = 0; in_fns = {}; out_fns = {}; file_type = {}; params = {}; gps_source = {};
sync_flag = {}; sync_fns = {}; sync_file_type = {}; sync_params = {};

%% <== CHOOSE WHICH GPS SOURCE TO PROCESS
gps_source_to_use = 'arena';
% gps_source_to_use = 'cresis';

if strcmpi(gps_source_to_use,'arena')
  %% ARENA GPS SOURCE
  % =======================================================================

%   year = 2022; month = 12; day = 25;
%   datestr_year = 2022; datestr_month = 12; datestr_day = 24; % <--- UPDATE TO MATCH WHAT PREPROCESS PRINTS OUT
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = get_filenames(fullfile(in_base_path,sprintf('%04d%02d%02d',year,month,day)),'','','gps.txt');
%   out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat', datestr_year, datestr_month, datestr_day);
%   %date_str{file_idx} = sprintf('%04d%02d%02d', year, month, day);
%   date_str{file_idx} = sprintf('%04d%02d%02d', datestr_year, datestr_month, datestr_day);
%   file_type{file_idx} = 'arena';
%   params{file_idx} = struct('time_reference','utc');
%   gps_source{file_idx} = 'arena-field';
%   sync_flag{file_idx} = 1;
%   sync_fns{file_idx} = get_filenames(fullfile(in_base_path,sprintf('%04d%02d%02d',year,month,day)),'','','gps.txt');
%   sync_file_type{file_idx} = 'arena';
%   sync_params{file_idx} = struct('time_reference','utc');

  year = 2023; month = 1; day = 29;
  datestr_year = 2023; datestr_month = 1; datestr_day = 29; % <--- UPDATE TO MATCH WHAT PREPROCESS PRINTS OUT
  file_idx = file_idx + 1;
  in_fns{file_idx} = get_filenames(fullfile(in_base_path,sprintf('%04d%02d%02d',year,month,day)),'','','gps.txt');
  out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat', datestr_year, datestr_month, datestr_day);
  %date_str{file_idx} = sprintf('%04d%02d%02d', year, month, day);
  date_str{file_idx} = sprintf('%04d%02d%02d', datestr_year, datestr_month, datestr_day);
  file_type{file_idx} = 'arena';
  params{file_idx} = struct('time_reference','utc');
  gps_source{file_idx} = 'arena-field';
  sync_flag{file_idx} = 1;
  sync_fns{file_idx} = get_filenames(fullfile(in_base_path,sprintf('%04d%02d%02d',year,month,day)),'','','gps.txt');
  sync_file_type{file_idx} = 'arena';
  sync_params{file_idx} = struct('time_reference','utc');
  
elseif strcmpi(gps_source_to_use,'cresis')
  %% CReSIS GPS SOURCE
  % =======================================================================

  % DO NOT USE IN FIELD
end

%% gps_create
% Read and translate files according to user settings
% =========================================================================
gps_create;

%% custom fixes
% =========================================================================
for idx = 1:length(file_type)
  out_fn = fullfile(gps_path,out_fns{idx});
  
  load(out_fn,'gps_source');
  if ~isempty(regexpi(gps_source,'arena'))
    % Extrapolation is necessary because GPS data starts after/stops before
    % the beginning/end of the radar data.
    warning('Extrapolating and filtering elevation for arena GPS data: %s', out_fn);
    gps = load(out_fn);
    
    if length(gps.lat) >= 2
      new_gps_time = [gps.gps_time(1)-10, gps.gps_time,gps.gps_time(end)+10];
      gps.lat = interp1(gps.gps_time,gps.lat,new_gps_time,'linear','extrap');
      gps.lon = interp1(gps.gps_time,gps.lon,new_gps_time,'linear','extrap');
      gps.elev = interp1(gps.gps_time,gps.elev,new_gps_time,'linear','extrap');
      gps.roll = interp1(gps.gps_time,gps.roll,new_gps_time,'linear','extrap');
      gps.pitch = interp1(gps.gps_time,gps.pitch,new_gps_time,'linear','extrap');
      gps.heading = interp1(gps.gps_time,gps.heading,new_gps_time,'linear','extrap');
      gps.gps_time = new_gps_time;
      
      gps.elev = fir_dec(gps.elev,ones(1,101)/101,1);
      
      good_idxs = gps.radar_time~=0;
      gps.radar_time = gps.radar_time(good_idxs);
      gps.comp_time = gps.comp_time(good_idxs);
      gps.sync_elev = gps.sync_elev(good_idxs);
      gps.sync_gps_time = gps.sync_gps_time(good_idxs);
      gps.sync_lat = gps.sync_lat(good_idxs);
      gps.sync_lon = gps.sync_lon(good_idxs);
      
      save(out_fn,'-struct','gps');
    end
  end
  
  if ~isempty(regexpi(gps_source,'arena')) && ~isempty(regexpi(out_fn,'20221209'))
    [out_fn_dir,out_fn_name,out_fn_ext] = fileparts(out_fn);
    new_out_fn = fullfile(out_fn_dir,[out_fn_name(1:end-2),'10',out_fn_ext]);
    warning('Copying 20221209 to 20221210:\n  %s\n  %s\n', out_fn, new_out_fn);
    copyfile(out_fn,new_out_fn);
    new_out_fn = fullfile(out_fn_dir,[out_fn_name(1:end-2),'11',out_fn_ext]);
    warning('Copying 20221209 to 20221211:\n  %s\n  %s\n', out_fn, new_out_fn);
    copyfile(out_fn,new_out_fn);
    new_out_fn = fullfile(out_fn_dir,[out_fn_name(1:end-2),'12',out_fn_ext]);
    warning('Copying 20221209 to 20221212:\n  %s\n  %s\n', out_fn, new_out_fn);
    copyfile(out_fn,new_out_fn);
  end
 
  if ~isempty(regexpi(out_fn,'20220419'))
    % Fix bad radar time 
    warning('Fixing bad radar_time: %s', out_fn);
    gps = load(out_fn);
    
    if gps.radar_time(1) == 0
      gps.radar_time = gps.radar_time(2:end);
      gps.sync_elev = gps.sync_elev(2:end);
      gps.sync_gps_time = gps.sync_gps_time(2:end);
      gps.sync_lat = gps.sync_lat(2:end);
      gps.sync_lon = gps.sync_lon(2:end);
    end
    save(out_fn,'-append','-struct','gps','radar_time','sync_elev','sync_gps_time','sync_lat','sync_lon');
  end

  if ~isempty(regexpi(out_fn,'201910XX'))
    % Fake GPS for testing
    warning('Faking GPS data: %s', out_fn);
    gps = load(out_fn);
    
    velocity = 4;
    gps.lat = -75.5 - (gps.gps_time-gps.gps_time(1))*velocity/111111;
    gps.lon(:) = -106.75;
    gps.elev(:) = 500;
    gps.heading(:) = -pi;
    
    save(out_fn,'-append','-struct','gps','gps_time','lat','lon','elev','roll','pitch','heading');
  end
  
end
