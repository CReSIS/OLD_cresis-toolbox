%% ICESSN ATM Lidar to XYZ CSV
%
% Reads in ICESSN data and outputs to an XYZ csv file.
%
% Output is CSV/TXT file with fields:
% 'DATE','GPSTIME','LAT','LON','SURFACE'
%
% Author: Kyle Purdon
%
% See also read_lidar_atm.m

%% User Input

% Specify an output file (Full path with filename and extention)
% Ex.('C:\Users\SomeUser\SomeFolder\SomeFilename.txt')
output_file = 'C:\Users\kpurdon\Documents\Projects\Other\jaktest\2008_G_TO_0711_02_010_atm2.txt';

% Specify the seasons of data you want to extract from.
% Ex.(season_list{end+1,1} = '1998_Greenland_P3';)
season_list = {};
season_list{end+1,1} = '2008_Greenland_TO';

% Specify the segments of data you want to extract.
% Leave this empty to extract full seasons.
% You only need to include year/month/day. (Ex. '20080725')
seg_list = {};
seg_list{end+1,1} = '20080711';
% seg_list{end+1,1} = '20010521';

%% Automated Section

% Set the ATM Data Path
if ispc
  atmpath = '\\titan\projects\metadata\ATM_smooth_nadir\';
else
  atmpath = '/cresis/projects/metadata/ATM_smooth_nadir/';
end

fprintf('=============================\n');
fprintf('         ATM TO CSV          \n');
fprintf('=============================\n');

%% Data Search

fprintf('Searching for ATM data ... ');
% Check if seasons given have ATM data
error_chk = 0;
atm_dirs = strcat(atmpath,season_list,filesep);
for s_chk = 1:length(season_list)
  if ~exist(atm_dirs{s_chk},'dir')
    fprintf('\nSeason %s does not have ATM data.\n',season_list{s_chk});
    error_chk = 1;
  end
end
if error_chk == 1
  error('Season/s with no atm data. Please remove the season.');
end
fprintf('%s\n',datestr(now,13));


%% Filename Search and Year list creation

fprintf('Parsing ATM data paths ... ');
% Get all the ATM filenames and years (Remove bad fns (.swf files))
year_list = cell(length(season_list),1);
for year_idx = 1:length(season_list)
  % Fill a "year only" list
  year_list{year_idx} = season_list{year_idx}(1:4);
end
atm_years = {};
atm_fns = {};
for dir_idx = 1:length(atm_dirs)
  % Fill a filenames array
  curr_atm_fns = get_filenames(atm_dirs{dir_idx},'','smooth_nadir','','recursive');
  atm_fns = cat(1,atm_fns,curr_atm_fns);
  % Fill a years array
  curr_year_list = cell(length(curr_atm_fns),1);
  for year_idx = 1:length(curr_atm_fns)
    curr_year_list{year_idx} = year_list{dir_idx};
  end
  atm_years = cat(1,atm_years,curr_year_list);
end
bad_file_idxs = zeros(length(atm_fns));
for file_idx = 1:length(atm_fns)
  if ~isempty(strfind(atm_fns{file_idx},'.'))
    bad_file_idxs(file_idx) = 1;
  end
end
atm_fns(logical(bad_file_idxs)) = [];
atm_years(logical(bad_file_idxs)) = [];
fprintf('%s\n',datestr(now,13));

%% Filter ATM data

fprintf('Filtering ATM data ... ');
% Filter ATM filenames by seg_list
sep_idxs = strfind(atm_fns,filesep);
atm_segs = cell(length(atm_fns),1);
for sep_idx = 1:length(sep_idxs)
  atm_segs{sep_idx} = atm_fns{sep_idx}(sep_idxs{sep_idx}(end)+1:sep_idxs{sep_idx}(end)+6);
end
good_idxs = zeros(length(atm_segs),1);
for seg_filter = 1:length(seg_list)
  for atm_filter = 1:length(atm_segs)
    if strcmp(atm_segs{atm_filter},seg_list{seg_filter}(3:end))
      good_idxs(atm_filter,1) = 1;
    end
  end
end
good_idxs = logical(good_idxs);
atm_segs = atm_segs(good_idxs);
atm_fns = atm_fns(good_idxs);
atm_years = atm_years(good_idxs);
fprintf('%s\n',datestr(now,13));

%% Prepare to read ATM data

fprintf('Prepping ATM data ... ');
fprintf('%s\n',datestr(now,13));

%% Read ATM data

fprintf('Reading ATM data ... ');
% Read ATM lidar data
lidar = read_lidar_atm(atm_fns);
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

% Get Year/Month/Day from GPS Time
date_num = epoch_to_datenum(lidar.gps_time);
date = num2str(datestr(date_num,2));

fprintf('%s\n',datestr(now,13));

%% Write Output File

fprintf('Writing output file ... ');
% Write Output File
fid = fopen(output_file,'w+');
fprintf(fid,'%s,%s,%s,%s,%s\n','DATE','GPSTIME','LAT','LON','SURFACE');
for out_idx = 1:length(lidar.gps_time)
  fprintf(fid,'%s,%f,%2.6f,%2.6f,%6.4f\n',...
    date(out_idx,:),lidar.gps_time(out_idx),lidar.lat(out_idx),lidar.lon(out_idx),lidar.surface(out_idx));
end
cid = fclose(fid);
fprintf('%s\n',datestr(now,13));

fprintf('ATM TO CSV Complete.\n')

