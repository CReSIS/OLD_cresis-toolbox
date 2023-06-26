function newdatacell = sync_atm_data(csvdatacell,atmpath,maxgap)
% newdatacell = sync_atm_data(csvdatacell)
%
% Syncs a CReSIS CSV file with ATM Lidar Data.
%
% Input
%   csvdatacell: Cell stucture containing fields:
%         "LAT,LON,TIME,THICK,ELEV,FRAME,SURF,BOTT,QUALITY,SEASON"
%   atmpath: Path to root directory of ATM smooth_nadir data
%   maxgap: The maximum gap in meters to interpolate ATM data (500 is good)
%
% Output
%     New cell structure with A_SURF,A_BED, and DATATYPE added.   
%
% Author: Kyle Purdon
% Center For Remote Sensing of Ice Sheets (CReSIS)
% University of Kansas
%
% Example:
% cresis_csv_fn = 'C:\SomeFolder\Data_20100412_01_008_190146.csv;
% fid = fopen(cresis_csv_fn,'r');
% csvdatacell = textscan(fid,'%f%f%f%f%f%s%f%f%d%s','delimiter',',','headerlines',1);
% fclose(fid);
% atmpath = 'P:\metadata\ATM_smooth_nadir\';
% maxgap = 500;
% newcsvdatacell = sync_atm_data(csvdatacell,atmpath,maxgap);
%


%% Deal CSV Data Structure
[LATITUDE LONGITUDE UTCTIMESOD THICKNESS ELEVATION FRAME SURFACE BOTTOM QUALITY SEASON] = deal(csvdatacell{:});

%% Check for ATM Data by Season

fprintf('Preparing for ATM Sync ... ');
% Get unique seasons list and the indexes
[season_list, season_idx_end] = unique(SEASON,'last');
[tmp, season_idx_start] = unique(SEASON,'first');

% Get a list of seasons with ATM data to sync (NOT ALL THE SAME)
atm_dirs = strcat(atmpath,season_list,filesep);
atm_sync_list = zeros(length(atm_dirs),1);
for dir_idx = 1:length(atm_dirs)
  atm_sync_list(dir_idx) = logical(exist(atm_dirs{dir_idx},'dir'));
end
fprintf('%s\n',datestr(now,'HH:MM:SS'));

fprintf('Sorting CSV Data by Season ... ');
% Split CSV variables by season (store in cell arrays)
num_seasons = length(season_list);
lat = cell(num_seasons,1);        lon = cell(num_seasons,1);
utcsod = cell(num_seasons,1);     thick = cell(num_seasons,1);
elev = cell(num_seasons,1);       frame = cell(num_seasons,1);
surf = cell(num_seasons,1);       bott = cell(num_seasons,1);
quality = cell(num_seasons,1);    season = cell(num_seasons,1);

for s_idx = 1:length(season_list)
  
  % Pre-Allocate Each run
  csl = season_idx_end(s_idx)-season_idx_start(s_idx);
  lat{s_idx} = zeros(csl,1);
  lon{s_idx} = zeros(csl,1);
  utcsod{s_idx} = zeros(csl,1);
  thick{s_idx} = zeros(csl,1);
  elev{s_idx} = zeros(csl,1);
  frame{s_idx} = zeros(csl,1);
  surf{s_idx} = zeros(csl,1);
  bott{s_idx} = zeros(csl,1);
  quality{s_idx} = zeros(csl,1);
  season{s_idx} = zeros(csl,1);
  
  % Fill each Season with Data
  lat{s_idx}= LATITUDE(season_idx_start(s_idx):season_idx_end(s_idx),1);
  lon{s_idx} = LONGITUDE(season_idx_start(s_idx):season_idx_end(s_idx),1);
  utcsod{s_idx} = UTCTIMESOD(season_idx_start(s_idx):season_idx_end(s_idx),1);
  thick{s_idx} = THICKNESS(season_idx_start(s_idx):season_idx_end(s_idx),1);
  elev{s_idx} = ELEVATION(season_idx_start(s_idx):season_idx_end(s_idx),1);
  frame{s_idx} = FRAME(season_idx_start(s_idx):season_idx_end(s_idx),1);
  surf{s_idx} = SURFACE(season_idx_start(s_idx):season_idx_end(s_idx),1);
  bott{s_idx} = BOTTOM(season_idx_start(s_idx):season_idx_end(s_idx),1);
  quality{s_idx} = QUALITY(season_idx_start(s_idx):season_idx_end(s_idx),1);
  season{s_idx} = SEASON(season_idx_start(s_idx):season_idx_end(s_idx),1);
end

% Create new, empty variable cell arrays (A_SURF,A_BED,DATATYPE)
a_surf = cell(num_seasons,1);
a_bed = cell(num_seasons,1);
data_type = cell(num_seasons,1);

%% ITERATE THROUGH EACH SEASON (SYNC ATM OR RADAR)

fprintf('%s\n',datestr(now,'HH:MM:SS'));
fprintf('May take awhile, please be patient. \n');

for season_idx_end = 1:num_seasons
  if atm_sync_list(season_idx_end)
    % ============= READ LIDAR DATA =============
    fprintf('%d/%d  ',season_idx_end,num_seasons);
    fprintf('Syncing %s with ATM data ... ',sprintf('%s',season{season_idx_end}{1}));
    
    % Get all the ATM filenames
    atm_fns = get_filenames(atm_dirs{season_idx_end},'','smooth_nadir','','recursive');
    
    % Check for Bad Filenames (Gets rid of .swf files)
    bad_file_idxs = [];
    for file_idx = 1:length(atm_fns)
      if ~isempty(strfind(atm_fns{file_idx},'.'))
        bad_file_idxs(end+1) = file_idx;
      end
    end
    atm_fns(bad_file_idxs) = [];
    
    % Read ATM lidar data
    lidar = read_lidar_atm(atm_fns);
    
    % Remove NAN's from LIDAR Data
    good_lidar_idxs = (~isnan(lidar.gps_time)) & (lidar.rms<=50);
    lidar.gps_time = lidar.gps_time(good_lidar_idxs);
    lidar.surface = lidar.surface(good_lidar_idxs);
    lidar.lat = lidar.lat(good_lidar_idxs);
    lidar.lon = lidar.lon(good_lidar_idxs);
    if ~isempty(isnan(lidar.surface))
      good_lidar_idxs = ~isnan(lidar.surface);
      lidar.gps_time = lidar.gps_time(good_lidar_idxs);
      lidar.surface = lidar.surface(good_lidar_idxs);
      lidar.lat = lidar.lat(good_lidar_idxs);
      lidar.lon = lidar.lon(good_lidar_idxs);
    end
    
    % ============= FIND LARGE GAPS IN LIDAR =============
        
    % Find Gaps Based on Along_Track Distance
    lidar.along_track = geodetic_to_along_track(lidar.lat,lidar.lon,lidar.surface);
    lidar.lidar_gaps = (diff(lidar.along_track) > maxgap);
    
    % Prepare Data for Interp1 if gaps exist
    if ~isempty(lidar.lidar_gaps)
      % Find Start Indexes of Gaps
      lidar.lidar_gap_start_idxs = strfind([0 lidar.lidar_gaps 0],[0 1]);
      lidar.lidar_gap_start_idxs = lidar.lidar_gap_start_idxs-1;
      if lidar.lidar_gap_start_idxs(1) == 0
        lidar.lidar_gap_start_idxs(1) = 1;
      end
      
      % Find End Indexes of Gaps
      lidar.lidar_gap_end_idxs = strfind([0 lidar.lidar_gaps 0],[1 0]);
      if lidar.lidar_gap_end_idxs(end) > length(lidar.lidar_gaps)
        lidar.lidar_gap_end_idxs(end) = length(lidar.lidar_gaps);
      end
      
      % Add NaN's to lidar at Gap_IDXS (Start & End)
      lidar.surface([lidar.lidar_gap_start_idxs lidar.lidar_gap_end_idxs]) = NaN;
    end
    
    % ============= SET UP RADAR DATA FOR SYNC =============
    
    % Get GPS time from UTC time for the current season
    radar.utctime=[]; radar.gpstime=[];
    radar.year = zeros(length(frame{season_idx_end}),1);
    radar.month = zeros(length(frame{season_idx_end}),1);
    radar.day = zeros(length(frame{season_idx_end}),1);
    % Populate Year/Month/Day
    for data_idx = 1:length(frame{season_idx_end})
      radar.year(data_idx) = str2double(frame{season_idx_end}{data_idx}(1:4));
      radar.month(data_idx) = str2double(frame{season_idx_end}{data_idx}(5:6));
      radar.day(data_idx) = str2double(frame{season_idx_end}{data_idx}(7:8));
    end
    % Get UTC and GPS time for radar_data
    radar.utctime = datenum_to_epoch(datenum(radar.year,radar.month,radar.day,0,0,utcsod{season_idx_end}));
    radar.gpstime = radar.utctime + utc_leap_seconds(radar.utctime(1));
    
    % Calculate real RADAR surface elevation.
    radar.real_surface = elev{season_idx_end}-surf{season_idx_end};
    
    %Test for and remove duplicates in RADAR
    if length(unique(radar.gpstime)) < length(radar.gpstime)
      [radar.gpstime radar.unique_idxs] = unique(radar.gpstime);
      radar.real_surface = radar.real_surface(radar.unique_idxs);
      % Set Unique RADAR Variables
      lat{season_idx_end} = lat{season_idx_end}(radar.unique_idxs);
      lon{season_idx_end} = lon{season_idx_end}(radar.unique_idxs);
      utcsod{season_idx_end} = utcsod{season_idx_end}(radar.unique_idxs);
      thick{season_idx_end} = thick{season_idx_end}(radar.unique_idxs);
      elev{season_idx_end} = elev{season_idx_end}(radar.unique_idxs);
      frame{season_idx_end} = frame{season_idx_end}(radar.unique_idxs);
      surf{season_idx_end} = surf{season_idx_end}(radar.unique_idxs);
      bott{season_idx_end} = bott{season_idx_end}(radar.unique_idxs);
      quality{season_idx_end} = quality{season_idx_end}(radar.unique_idxs);
      season{season_idx_end} = season{season_idx_end}(radar.unique_idxs);
    end
    
    %Test for and remove duplicates in LIDAR
    if length(unique(lidar.gps_time)) < length(lidar.gps_time)
      [lidar.gps_time lidar.unique_idxs] = unique(lidar.gps_time);
      lidar.surface = lidar.surface(lidar.unique_idxs);
    end
    
    % ============= INTERPOLATE LIDAR TO RADAR =============
    
    % Get lidar_surface for the radar points (radar.gpstime)
    warning off;
    lidar.newsurface = interp1(lidar.gps_time,lidar.surface,radar.gpstime,'linear');
    
    % ============= INSERT RADAR INTO LIDAR GAPS =============
    
    % Get the indexes if newsurface to replace.
    newsurface_replace_idxs = isnan(lidar.newsurface);
    
    % Fill in the newsurface at the replacement points with radar surface.
    lidar.newsurface(newsurface_replace_idxs)  = interp1(radar.gpstime,radar.real_surface,...
      radar.gpstime(newsurface_replace_idxs));
    
    % ============= POPULATE NEW VARIABLES =============
    
    % Populate a_surf and a_bed
    a_surf{season_idx_end} = lidar.newsurface;
    a_bed{season_idx_end} = a_surf{season_idx_end} - thick{season_idx_end};
    
    % Populate data_type
    for d_idx = 1:length(a_surf{season_idx_end})
      data_label{d_idx} = 'lidar';
    end
    data_type{season_idx_end} = data_label;
    fprintf('%s\n',datestr(now,'HH:MM:SS'));
      
  else
    % ============= SYNC WITH RADAR DATA =============
    fprintf('%d/%d ',season_idx_end,num_seasons);
    fprintf('Syncing %s with RADAR data ... ',sprintf('%s',season{season_idx_end}{1}));
    
    % Populate a_surf and a_bed
    a_surf{season_idx_end} = elev{season_idx_end} - surf{season_idx_end};
    a_bed{season_idx_end} = elev{season_idx_end} - bott{season_idx_end};
    
    % Populate data_type
    data_label = cell(length(a_surf{season_idx_end}),1);
    for d_idx = 1:length(a_surf{season_idx_end})
      data_label{d_idx} = 'radar';
    end
    data_type{season_idx_end} = data_label;
    
    % Complete
    fprintf('%s\n',datestr(now,'HH:MM:SS'));
    
    
  end
end

%% Populate New Cell Structure
newdatacell = cell(1,13);
% Pre-Allocate newdatacell
for pre_idx = [1,2,3,4,5,7,8,9,11,12]
  newdatacell{pre_idx} = zeros(length(lat),1);
end

fprintf('Filling new cell structure ... ');
newdatacell{1} = cell2mat(lat);
newdatacell{2} = cell2mat(lon);
newdatacell{3} = cell2mat(utcsod);
newdatacell{4} = cell2mat(thick);
newdatacell{5} = cell2mat(elev);
newdatacell{7} = cell2mat(surf);
newdatacell{8} = cell2mat(bott);
newdatacell{9} = cell2mat(quality);
newdatacell{11} = cell2mat(a_surf);
newdatacell{12} = cell2mat(a_bed);
% Handle cells of strings.
for s_idx = 1:length(frame)
  for c_idx = 1:length(frame{s_idx})
    newdatacell{6}{end+1,1} = frame{s_idx}{c_idx};
    newdatacell{10}{end+1,1} = season{s_idx}{c_idx};
    newdatacell{13}{end+1,1} = data_type{s_idx}{c_idx};
  end
end
fprintf('%s\n',datestr(now,'HH:MM:SS'));
end





