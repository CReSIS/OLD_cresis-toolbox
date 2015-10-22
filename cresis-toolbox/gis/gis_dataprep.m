function gis_dataprep(param)
%
% gis_dataprep(param)
%
% Prepares a CSV file created using geographic_search_gui.m for the ArcGIS gridding process.
% The script preforms three main functions.
% (1) Removes any high elevation data. (Elevation >= 10,000ft)
% (2) Removes any data with zeros ice thickness (Thick == 0m)
% (3) Keeps only specified years of data.
% (4) Splits FRAME for ArcGIS Shapefile Error
% (5) Sync NASA ATM Lidar Data
%
% Input Param Structure
%   .input = csv file from atm_sync_to_geocsv.m
%   .remove_high_elev = logical true or false
%   .remove_negative_thick = logical true or false
%   .remove_derived = logical true or false
%   .split_frame = logical true of false
%   .keep_years = cell of seasons to keep. (Empty keeps all)
%   .sync_atm= logical true or false
%   .output = output path for new csv/txt file (With filename and extension)
%
% Output
%   CSV/TXT file ready for ArcGIS Gridding
%
% Additional Information:
% An example season is '2011_Greenland_P3' the keep year would be '2011'
% This format must be used exactly, or the script will not work.
%
% Author: Kyle Purdon
%
% Example:
% ---------------------------------------
% Read in a file with ALL the data over the NWCoast of Greenland, remove
% the high elevation data and zero thickness and keep only data from 2010
% and 2011. Also remove derived quality data, split frame, and sync ATM.
% ---------------------------------------
% param.input = 'C:\Users\kpurdon\projects\greenland\Greenland_NWCoast_AllYears.csv';
% param.remove_high_elev = true;
% param.remove_negative_thick = true;
% param.remove_derived = true;
% param.split_frame = true;
% param.keep_years = {'2011','2010'};
% param.sync_atm = true;
% param.output = 'C:\Users\kpurdon\projects\greenland\Greenland_NWCoast_2010_2011.csv'
% gisdataprep(param);
% ---------------------------------------

%% Print Function Header

fprintf('\n-------------------------------\n');
fprintf('      GIS DATA PREPERATION        \n');
fprintf('-------------------------------\n\n');


%% Load/Read/Close GeoSearch GUI CSV

% Open CSV
fprintf('Reading Input CSV file ... ');
fid = fopen(param.input,'r');
if fid == -1
  error('Cannot Open CSV Input File. Check Path.');
end

% Read CSV
csvdata = textscan(fid,'%f%f%f%f%f%s%f%f%d%s','delimiter',',','headerlines',1);

% Close CSV
cid = fclose(fid);
if cid == -1
  fprintf('CSV File did not close.\n')
end
fprintf('%s\n',datestr(now,'HH:MM:SS'));

%% Get Good Indexes
fprintf('Finding Good Indexes ... ');
if ~isempty(param.keep_years)
  datayears = cell(1,length(csvdata{10}));
  for season_idx = 1:length(csvdata{10})
    datayears{1,season_idx} = csvdata{10}{season_idx}(1:4);
  end
  good_year_idxs = find(ismember(datayears,param.keep_years)==1);
else
  good_year_idxs = 1:length(csvdata{10});
end
if param.remove_high_elev
  good_elev_idxs = find(csvdata{7} < 3048)';
else
  good_elev_idxs = 1:length(csvdata{7});
end
if param.remove_negative_thick
  good_thick_idxs = find(csvdata{4} >= 0)';
else
  good_thick_idxs = 1:length(csvdata{5});
end
if param.remove_derived
  good_quality_idxs = find(csvdata{9} < 3)';
else
  good_quality_idxs = 1:length(csvdata{5});
end

good_idxs = intersect(intersect(intersect(good_elev_idxs,good_thick_idxs),good_quality_idxs),good_year_idxs);
fprintf('%s\n',datestr(now,'HH:MM:SS'));

%% Remove Bad Data
fprintf('Removing Data ... ');
for data_idx = 1:length(csvdata);
  csvdata{data_idx} = csvdata{data_idx}(good_idxs);
end
fprintf('%s\n',datestr(now,'HH:MM:SS'));

%% SYNC Atm Data
if param.sync_atm  
  fprintf('\nSyncing ATM Data\n');
  fprintf('---------------------------------\n');  
  newdata = sync_atm_data(csvdata,'P:\metadata\ATM_smooth_nadir\',500);
  fprintf('---------------------------------\n');
else
  newdata = csvdata;
end

%% Split Frame
if param.split_frame
  fprintf('Splitting Frame ... ');
  [yyyymmdd,segment,frame] = split_frame(strrep(newdata{6},'_',''));
  fprintf('%s\n',datestr(now,'HH:MM:SS'));
end

%% Write New File (STANDARD OR SPLIT FRAME)(ATM OR NO ATM)

fprintf('\nWriting New File ... ');
fid_out = fopen(param.output,'w');
if fid_out == -1
  error('Cannot Open Output File. Check Path.');
end

if ~param.sync_atm % NO ATM Sync
  if ~param.split_frame % Standard FRAME
    % Print Header
    fprintf(fid_out,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
      'LAT','LON','UTCTime','THICK','ELEVATION','FRAME','SURFACE','BOTTOM','QUALITY','SEASON');
    % Write Data
    for data_idx = 1:length(newdata{1})
      fprintf(fid_out,'%2.6f,%2.6f,%5.4f,%6.2f,%4.4f,%s,%6.2f,%6.2f,%01d,%s,\n',...
        newdata{1}(data_idx),newdata{2}(data_idx),newdata{3}(data_idx),newdata{4}(data_idx),newdata{5}(data_idx),...
        newdata{6}{data_idx},newdata{7}(data_idx),newdata{8}(data_idx),newdata{9}(data_idx),newdata{10}{data_idx});
    end
  else % Split FRAME
    % Print Header
    fprintf(fid_out,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
      'LAT','LON','UTCTime','THICK','ELEVATION','YYYYMMDD','SEGMENT','FRAME','SURFACE',...
      'BOTTOM','QUALITY','SEASON');
    % Write Data
    for data_idx = 1:length(newdata{1})
      fprintf(fid_out,'%2.6f,%2.6f,%5.4f,%6.2f,%4.4f,%s,%s,%s,%6.2f,%6.2f,%01d,%s,\n',...
        newdata{1}(data_idx),newdata{2}(data_idx),newdata{3}(data_idx),newdata{4}(data_idx),newdata{5}(data_idx),...
        yyyymmdd{data_idx},segment{data_idx},frame{data_idx},newdata{7}(data_idx),newdata{8}(data_idx),newdata{9}(data_idx),newdata{10}{data_idx});
    end
  end
else % ATM Sync
  if ~param.split_frame % Standard FRAME
    % Print Header
    fprintf(fid_out,'%s,%s,%s,%s,%s,%s,%s,%s,%s%s,%s,%s,%s\n',...
      'LAT','LON','UTCTime','THICK','ELEVATION','FRAME','SURFACE',...
      'BOTTOM','QUALITY','SEASON','A_SURF','A_BED','DATATYPE');
    % Write Data
    for data_idx = 1:length(newdata{1})
      fprintf(fid_out,'%2.6f,%2.6f,%5.4f,%6.2f,%4.4f,%s,%6.2f,%6.2f,%01d,%s,%6.2f,%6.2f,%s\n',...
        newdata{1}(data_idx),...  %LAT 2.6f
        newdata{2}(data_idx),...  %LON 2.6f
        newdata{3}(data_idx),...  %Time 5.4g
        newdata{4}(data_idx),...  %Thick 6.2g
        newdata{5}(data_idx),...  %Elevation 4.4f
        newdata{6}{data_idx},...  %Frame s
        newdata{7}(data_idx),...  %Surface 6.2f
        newdata{8}(data_idx),...  %Bottom 6.2f
        newdata{9}(data_idx),...  %Quality 01d
        newdata{10}{data_idx},... %Season s
        newdata{11}(data_idx),... %A_SURF 6.2f
        newdata{12}(data_idx),... %A_BED 6.2f
        newdata{13}{data_idx});   %Datatype s;
    end
  else % Split FRAME
    % Print Header
    fprintf(fid_out,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
      'LAT','LON','UTCTime','THICK','ELEVATION','YYYYMMDD','SEGMENT','FRAME','SURFACE'...
      ,'BOTTOM','QUALITY','SEASON','A_SURF','A_BED','DATATYPE');
    % Write Data
    for data_idx = 1:length(newdata{1})
      fprintf(fid_out,'%2.6f,%2.6f,%5.4f,%6.2f,%4.4f,%s,%s,%s,%6.2f,%6.2f,%01d,%s,%6.2f,%6.2f,%s\n',...
        newdata{1}(data_idx),...  %LAT 2.6f
        newdata{2}(data_idx),...  %LON 2.6f
        newdata{3}(data_idx),...  %Time 5.4g
        newdata{4}(data_idx),...  %Thick 6.2g
        newdata{5}(data_idx),...  %Elevation 4.4f
        yyyymmdd{data_idx},...    %YYYYMMDD s
        segment{data_idx},...     %Segment s
        frame{data_idx},...       %Frame s
        newdata{7}(data_idx),...  %Surface 6.2f
        newdata{8}(data_idx),...  %Bottom 6.2f
        newdata{9}(data_idx),...  %Quality 01d
        newdata{10}{data_idx},... %Season s
        newdata{11}(data_idx),... %A_SURF 6.2f
        newdata{12}(data_idx),... %A_BED 6.2f
        newdata{13}{data_idx});   %Datatype s
    end
  end
end

% Close File
cidf = fclose(fid_out);
if cidf == -1
  fprintf('Output File release error.\n');
end
% Complete
fprintf('%s\n',datestr(now,'HH:MM:SS'));

fprintf('Done ... %s\n',datestr(now,'HH:MM:SS'));

clear cidf fid_out
end