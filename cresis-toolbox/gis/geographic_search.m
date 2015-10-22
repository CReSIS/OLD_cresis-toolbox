% script geographic_search
%
% Function for finding data within a geographic area with optional
% automatic download.
%
% While an arbitrary set of CSV files can be used to search for
% the data, the typical usage is to point to the browse files
% on the https site. The script is already setup to do Greenland
% and Antarctica.
%
% The automatic download grabs all of the matching data off of
% the https site and stores it into a CSV or MAT file.  Only Level-2
% (layer) data is downloaded.
%
% To download large files, you may need to set the Java Heap size
% in the Matlab preferences to 256 MB or more.
%
% If you do not have access to Matlab, please send email to
% cresis_data@cresis.ku.edu with the corner coordinates of interest
% and we can run the program for you.
%
% Outputs:
%   Plots the lat/lon for all data within the bounding box
%   Prints out the frame IDs and segments IDs with any data in the
%   bounding box.
%   All other varibles are just temporary variables
%   If the download data is set, the function will download all of
%   the data into the specified output directory and populate these
%   variables:
%     Latitude, Longitude, UTC_time_sod, Thickness, Elevation,
%     Frame, Surface, Bottom, Quality
%
% Example/Instructions:
%   Download all the season browse CSV files of interest OR specify
%      their URL
%   Fill the fns cell array with each of these filenames or URLs
%   Fill the seasons cell array with the corresponding season for each
%      CSV (this is only needed if download_en is set to true)
%   Specify the lat/lon limits
%      lon are -180 to +180 (areas which bridge -180 to +180
%      must be done as two searches)
%      lat are -90 to +90
%   Specify the radar name (mcords2, mcords, acords, icards, accum
%      snow, or kuband)
%   Specify if and where you would like the data downloaded to and
%      which data source to use (csv or csv_good)
%   Run the script
%
% Author: John Paden

% ===========================================================
% User Settings
% ===========================================================

% Latitude and longitude limits in degrees (order does not matter)
lon_limits = [-50.4 -49.5];
lat_limits = [67.28 67.05];

% radar_name = string with radar name (mcords2, mcords, acords, icards, accum
%      snow, or kuband)
radar_name = 'mcords';

% location = string with location (greenland, antarctica)
location = 'greenland';

% download_en = boolean (if true, data is downloaded)
download_en = true;

% download_out_fn = filename of where to store data, if extension
%   is 'mat' it will store in -v6 Matlab file, otherwise CSV
download_out_fn = 'geographic_search_tst.csv';        % FOR CSV OUTPUT
% download_out_fn = 'geographic_search_tst.mat';  % FOR MAT OUTPUT

% download_type = string of which csv directory to use ('csv' = all data,
%   'csv_good' = data with both surface and bottom picks)
download_type = 'csv_good';


% =========================================
% DO NOT MODIFY ANYTHING BELOW THIS LINE
% =========================================

% Set fns and seasons
fns = {};
seasons = {};
if strcmpi(location,'greenland')
  % Greenland and north east Canada
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/1993_Greenland_P3/csv/Browse_1993_Greenland_P3.csv';
  seasons{end+1} = '1993_Greenland_P3';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/1995_Greenland_P3/csv/Browse_1995_Greenland_P3.csv';
  seasons{end+1} = '1995_Greenland_P3';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/1996_Greenland_P3/csv/Browse_1996_Greenland_P3.csv';
  seasons{end+1} = '1996_Greenland_P3';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/1997_Greenland_P3/csv/Browse_1997_Greenland_P3.csv';
  seasons{end+1} = '1997_Greenland_P3';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/1998_Greenland_P3/csv/Browse_1998_Greenland_P3.csv';
  seasons{end+1} = '1998_Greenland_P3';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/1999_Greenland_P3/csv/Browse_1999_Greenland_P3.csv';
  seasons{end+1} = '1999_Greenland_P3';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2001_Greenland_P3/csv/Browse_2001_Greenland_P3.csv';
  seasons{end+1} = '2001_Greenland_P3';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2002_Greenland_P3/csv/Browse_2002_Greenland_P3.csv';
  seasons{end+1} = '2002_Greenland_P3';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2003_Greenland_P3/csv/Browse_2003_Greenland_P3.csv';
  seasons{end+1} = '2003_Greenland_P3';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2005_Greenland_TO/csv/Browse_2005_Greenland_TO.csv';
  seasons{end+1} = '2005_Greenland_TO';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2006_Greenland_TO/csv/Browse_2006_Greenland_TO.csv';
  seasons{end+1} = '2006_Greenland_TO';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2007_Greenland_P3/csv/Browse_2007_Greenland_P3.csv';
  seasons{end+1} = '2007_Greenland_P3';
  % fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2008_Greenland_TO/csv/Browse_2008_Greenland_TO.csv';
  % seasons{end+1} = '2008_Greenland_TO';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2008_Greenland_Ground/csv/Browse_2008_Greenland_Ground.csv';
  seasons{end+1} = '2008_Greenland_Ground';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2009_Greenland_TO/csv/Browse_2009_Greenland_TO.csv';
  seasons{end+1} = '2009_Greenland_TO';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2010_Greenland_DC8/csv/Browse_2010_Greenland_DC8.csv';
  seasons{end+1} = '2010_Greenland_DC8';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2010_Greenland_P3/csv/Browse_2010_Greenland_P3.csv';
  seasons{end+1} = '2010_Greenland_P3';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2011_Greenland_P3/csv/Browse_2011_Greenland_P3.csv';
  seasons{end+1} = '2011_Greenland_P3';
else
  % Antarctica
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2002_Antarctica_P3chile/csv/Browse_2002_Antarctica_P3chile.csv';
  seasons{end+1} = '2002_Antarctica_P3chile';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2004_Antarctica_P3chile/csv/Browse_2004_Antarctica_P3chile.csv';
  seasons{end+1} = '2004_Antarctica_P3chile';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2009_Antarctica_DC8/csv/Browse_2009_Antarctica_DC8.csv';
  seasons{end+1} = '2009_Antarctica_DC8';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2009_Antarctica_TO/csv/Browse_2009_Antarctica_TO.csv';
  seasons{end+1} = '2009_Antarctica_TO';
  fns{end+1} = 'https://data.cresis.ku.edu/data/rds/2010_Antarctica_DC8/csv/Browse_2010_Antarctica_DC8.csv';
  seasons{end+1} = '2010_Antarctica_DC8';
end

% ===========================================================
% Automated Section
% ===========================================================
fprintf('=======================================================\n');
fprintf('Geographic Search\n');
fprintf('=======================================================\n');

lat_limits = sort(lat_limits);
lon_limits = sort(lon_limits);

frm_list = {};
season_list = [];
figure(1); clf;
fprintf('Opening sources\n');
for fn_idx = 1:length(fns)
  fn = fns{fn_idx};
  
  fprintf(' %s: %s\n', datestr(now,'HH:MM:SS'), fn);
  if ~exist(fn,'file')
    try;
      data_ptr = urlread(fn);
    catch ME;
      warning('Does not exist as file or URL\n', fn);
      continue;
    end
  else
    data_ptr = fopen(fn,'r');
  end

  % Parse CSV file (data_ptr is either a file identifier or a string,
  % but textscan works with both)
  if ~isempty(findstr('rds',radar_name))
    % Radar depth sounder format:
    C = textscan(data_ptr, '%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
    [LAT,LON,TIME,THICK,ELEVATION,FRAME,SURFACE,BOTTOM,QUALITY] = deal(C{:});
  else
    % Other radars format:
    C = textscan(data_ptr, '%f%f%f%f%s','headerlines',1,'delimiter',',');
    [LAT,LON,ELEVATION,TIME,FRAME] = deal(C{:});
  end
  if isa(data_ptr,'double')
    % This is a file id so we need to close it
    fclose(data_ptr);
  end

  % Make sure LON values go from -180 to +180 deg
  LON = mod(LON + 180,360) - 180;

  % Find limits within the bounding box
  good_idxs = find(LAT>lat_limits(1) ...
    & LAT<lat_limits(2) & LON>lon_limits(1) & LON<lon_limits(2));
  
  % Plot the results
  figure(1); hold on;
  plot(LON(good_idxs),LAT(good_idxs),'.');

  % Store the frames corresponding to these good data points
  frm_list = cat(1,frm_list,FRAME(good_idxs));
  % Tag each of these frames with the index, fn_idx, into the seasons array
  season_list = cat(1,season_list,fn_idx*ones(length(good_idxs),1));
  
end
figure(1); hold off;
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
grid on;

% Only keep the unique frames and seasons indices for just those frames
[frm_list good_idxs] = unique(frm_list);
season_list = season_list(good_idxs);

% Print the frames list out and create the segment list from that
fprintf('Frame IDs\n');
seg_list = {};
for frm_idx = 1:length(frm_list)
  fprintf('%s_%s_%s\n', frm_list{frm_idx}(1:8), ...
    frm_list{frm_idx}(9:10), frm_list{frm_idx}(11:13));
  seg_list{end+1} = sprintf('%s_%s', frm_list{frm_idx}(1:8), ...
    frm_list{frm_idx}(9:10));
end

% Only keep the unique segments and seasons indices for just those segments
[seg_list good_idxs] = unique(seg_list);
season_list = season_list(good_idxs);

% Print the segment list out
fprintf('\nSegment IDs\n');
for seg_idx = 1:length(seg_list)
  fprintf('%s\n', seg_list{seg_idx});
end

if download_en
  fprintf('\nBeginning download\n');
  
  % Load each segment CSV from the https site, read in the CSV file, select
  % just the data in the bounding box and concatenate these data together
  Latitude = [];
  Longitude = [];
  UTC_time_sod = [];
  Thickness = [];
  Elevation = [];
  Frame = {};
  Surface = [];
  Bottom = [];
  Quality = [];
  Season = {};
  for seg_idx = 1:length(seg_list)
    fprintf(' %s: %s\n', datestr(now,'HH:MM:SS'), seg_list{seg_idx});
    if ~isempty(findstr('rds',radar_name))
      url = sprintf('https://data.cresis.ku.edu/data/rds/%s/%s/Data_%s.csv', ...
        seasons{season_list(seg_idx)}, download_type, seg_list{seg_idx});
      data_ptr = urlread(url);
    else
      error('Data download not supported for non-rds radars');
      url = sprintf('https://data.cresis.ku.edu/%s/%s/%s/Data_%s.csv', ...
        radar_name, seasons{season_list(seg_idx)}, download_type, ...
        seg_list{seg_idx})
      data_ptr = urlread(url);
    end
    if isempty(data_ptr)
      continue;
    end
    C = textscan(data_ptr, '%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
    [LAT,LON,TIME,THICK,ELEVATION,FRAME,SURFACE,BOTTOM,QUALITY] = deal(C{:});
    good_idxs = find(LAT>lat_limits(1) ...
      & LAT<lat_limits(2) & LON>lon_limits(1) & LON<lon_limits(2));
    Latitude = cat(1,Latitude,LAT(good_idxs));
    Longitude = cat(1,Longitude,LON(good_idxs));
    UTC_time_sod = cat(1,UTC_time_sod,TIME(good_idxs));
    Thickness = cat(1,Thickness,THICK(good_idxs));
    Elevation = cat(1,Elevation,ELEVATION(good_idxs));
    Frame = cat(1,Frame,FRAME(good_idxs));
    Surface = cat(1,Surface,SURFACE(good_idxs));
    Bottom = cat(1,Bottom,BOTTOM(good_idxs));
    Quality = cat(1,Quality,QUALITY(good_idxs));
    Season = cat(1,Season,repmat(seasons(season_list(seg_idx)),length(LAT(good_idxs)),1));
  end
  
  % Save the data into a file depending on the file extension supplied
  % in the user settings
  [tmp tmp download_out_fn_ext] = fileparts(download_out_fn);
  if strcmpi(download_out_fn_ext,'.mat')
    % Save to MAT -v6 file
    save(download_out_fn,'-v6','Latitude','Longitude','UTC_time_sod','Thickness','Elevation','Frame','Surface','Bottom','Quality','Season');
  else
    % Save to CSV file
    fid_csv = fopen(download_out_fn,'w');
    fprintf(fid_csv,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
      'LAT','LON','UTCTIMESOD','THICK','ELEVATION','FRAME','SURFACE','BOTTOM','QUALITY','SEASON');
    for txt_idx = 1:length(Latitude)
      fprintf(fid_csv,'%2.6f,%2.6f,%5.4f,%6.2f,%4.4f,%s,%6.2f,%6.2f,%01d,%s\n',...
        Latitude(txt_idx),Longitude(txt_idx),...
        UTC_time_sod(txt_idx),Thickness(txt_idx),...
        Elevation(txt_idx),Frame{txt_idx},Surface(txt_idx),Bottom(txt_idx),...
        Quality(txt_idx),Season{txt_idx});
    end
  end
end

return;
