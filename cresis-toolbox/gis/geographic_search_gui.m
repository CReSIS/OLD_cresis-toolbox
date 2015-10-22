% Script Information/Details
% ==================================================
% Interactive Geographic Search/Download
% Author: Kyle Purdon
% Contributors: John Paden, Aric Beaver
% Version: V3.7 01/05/2012
% Contact Information: cresis_data@cresis.ku.edu
% Contact Subject: "Geographic Search Tool"
% Toolboxes Used: Standard, Mapping, Stats, Images
% Known Bugs: None
% Planned Updates: None
% Additional Information:
% see also run_geographic_search_gui.m
% ==================================================

function geographic_search_gui(param)
try
  %% User Input and Error Checking
  
  if param.runType == 2
    %Interactive Input Mode
    
    clear param.location_en;
    clear param.radarID;
    clear param.download_en;
    clear param.download_dir;
    clear param.file_name;
    clear param.extID;
    clear param.download_tid;
    clear param.down_sample;
    clear param.dsamp_meters;
    clear param.image_en;
    clear param.geotiff_fn;
    
    fprintf('=======================================================\n');
    fprintf('Interactive Geographic Search\n');
    fprintf('=======================================================\n');
    fprintf('Please use the following prompts to enter the information needed.\n');
    fprintf('Pay attention to the example format for how to enter values.\n\n')
    fprintf('-----------------------------------------------------------------\n');
    
    % User enters location of interest.
    fprintf('Enter your location of interest.\n');
    fprintf('(1)=Greenland, (2)=Antarctica or (3)=Canada\n');
    while true
      param.location_en = input('What is your location of interest? (1,2,3) : ');
      if param.location_en == 1
        fprintf('Your location of interest is Greenland.\n');
        break
      elseif param.location_en == 2
        fprintf('Your location of interest is Antarctica.\n');
        break
      elseif param.location_en == 3
        fprintf('Your location of interest is Canada.\n');
        break
      else
        fprintf('\n------ Invalid Location Entered ------\n')
      end
    end
    fprintf('-----------------------------------------------------------------\n');
    
    % User enters radar of interest.
    fprintf('Enter the radar you want to use (Default: RDS)\n');
    fprintf('(1)=RDS, (2)=Snow or (3)=Accum (4)=KUBand\n');
    while true
      param.radarID = input('What radar data do you want? (1,2,3,4) : ');
      if param.radarID == 1
        fprintf('You selected the RDS radar.\n');
        break
      elseif param.radarID == 2
        fprintf('You selected the Snow radar.\n');
        break
      elseif param.radarID == 3
        fprintf('You selected the Accum radar.\n');
        break
      elseif param.radarID == 4
        fprintf('You selected the KUBand radar.\n');
        break
      else
        fprintf('\n------ Invalid Value Entered ------\n')
      end
    end
    
    if param.radarID == 1
      fprintf('-----------------------------------------------------------------\n');
      %Get the TEMP Base path
      while true
        if param.dataGetType == 2
          fprintf('Enter a base path for the TEMP directory. (No Spaces)\n');
          fprintf('Ex. C:\\Users\\Username\\AppData\\Local\\Temp\\ \n');
          fprintf('Ex. /Users/TMP/ \n');
          param.tPath = input('Enter the Path: ','s');
          param.tPath = fullfile(param.tPath,filesep);
          if ~exist(param.tPath,'dir')
            mkdir(param.tPath);
            break
          else
            break
          end
        else
          break;
        end
      end
      
      fprintf('-----------------------------------------------------------------\n');
      
      % Does the user want to download or print a list ?
      fprintf('This script can automatically download/merge all the files in your study area.\n');
      fprintf('(1)=Yes (2)=No\n');
      while true
        param.download = input('Download the files? (1,2) : ');
        if param.download == 1
          fprintf('The files WILL be downloaded.');
          param.download_en = true;
          break
        elseif param.download == 2
          fprintf('The files WILL NOT be downloaded.');
          param.download_en = false;
          break
        else
          fprintf('\n------Invalid Decision.(Must enter 1 or 2) ------\n');
        end
      end
      
      %If the user wants to download the files get the input directory.
      % Also get the download_type and check the input.
      if param.download_en
        fprintf('\n-----------------------------------------------------------------\n');
        % Get directory, check if exists.
        fprintf('Enter the directory for your output files below.\n\n');
        fprintf('The directory should not contain any spaces.\n')
        fprintf('Ex.(C:\\Users\\Projects\\OutputFolder\\)\n\n');
        while true
          param.download_dir = input('Enter the Output directory : ','s');
          param.dPath_check = exist(param.download_dir,'dir');
          if param.dPath_check == 0
            mkdir(param.download_dir);
            break
          else
            break
          end
        end
        param.download_dir = fullfile(param.download_dir,filesep);
        fprintf('-----------------------------------------------------------------\n');
        % Get filename
        fprintf('Enter the filename for your output file.\n');
        while true
          param.file_name = input('Enter the filename: ','s');
          if ischar(param.file_name)
            break;
          else
            fprintf('\n------ Invalid filename. Must be a string.(Word) ------\n');
          end
        end
        fprintf('-----------------------------------------------------------------\n');
        %Get file extension
        fprintf('Enter the extension for your output file.\n')
        fprintf('(1) = csv  (2) = mat  (3) = both\n');
        while true
          param.extID = input('Enter the extension (1,2,3): ');
          if param.extID == 1 || param.extID == 2 || param.extID == 3
            break
          else
            fprintf('\n------ Invalid extension. Must be 1,2,or 3. ------\n');
          end
        end
        fprintf('-----------------------------------------------------------------\n');
        % Get csv_good or csv, check input
        fprintf('Specify (1)="csv_good" or (2)="csv" for the following prompt.\n\n');
        fprintf('csv includes ALL data.\n');
        fprintf('csv_good includes only data with surface and bottom picks\n\n');
        while true
          param.dlID = input('Enter choice (1,2): ');
          if param.dlID == 1
            param.download_type = 'csv_good';
            break
          elseif param.dlID == 2
            param.download_type = 'csv';
            break
          else
            fprintf('\nInvalid Decision. (Must enter 1 or 2)\n');
          end
        end
        fprintf('-----------------------------------------------------------------\n');
        % User decides if they want to downsample the output file.
        fprintf('Would you like to downsample the output file? (Change Spacing)\n');
        fprintf('(1)=Yes (2)=No\n');
        while true
          param.down_sample = input('\nEnter choice (1,2): ');
          if param.down_sample == 1
            param.down_sample = true;
            break
          elseif param.down_sample == 2
            param.down_sample = false;
            break
          else
            fprintf('\n------- Invalid Choice -------\n')
          end
        end
        fprintf('-----------------------------------------------------------------\n');
        if param.down_sample
          fprintf('What sample spacing would you like? (Meters)\n');
          while true
            param.dsamp_meters = input('\nEnter desired spacing: ');
            if isnumeric(param.dsamp_meters)
              break
            else
              fprintf('\n ------- Invalid number -------\n')
            end
          end
        end
      end
      
    end
    
    %Specify path to GeoTIFF or go with the default.
    fprintf('\n-----------------------------------------------------------------\n');
    fprintf('This script allows you to use any raster image for the background.\n');
    fprintf('You may specify your own local or web path, OR, use the default images provided.\n\n');
    fprintf('------- Default Images -------\n\n');
    fprintf('\tGreenland: Landsat-7 Greenland_natural_150m.tif\n');
    fprintf('\tCanada: Landsat-7 Canada_250m.tif\n');
    fprintf('\tAntarctica: Landsat-7 Antarctica_LIMA_480m.tif\n');
    fprintf('\n------- Example User Path -------\n\n');
    fprintf('\tC:\\Users\\Images\\someimage.tif\n');
    fprintf('\nFor default enter 1, for custom enter 2.\n')
    while true
      param.image_en = input('Enter choice (1,2): ');
      if param.image_en == 1
        if param.dataGetType == 1
          %Local paths to the GeoTIFFs
          if ispc
            switch param.location_en
              case 1;
                param.geotiff_fn = '\\titan\projects\GIS_data\greenland\Landsat-7\Greenland_natural_150m.tif';
                break
              case 2;
                param.geotiff_fn = '\\titan\projects\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif';
                break
              case 3;
                param.geotiff_fn = '\\titan\projects\GIS_data\canada\Landsat-7\Canada_250m.tif';
                break
            end
          else
            switch param.location_en
              case 1;
                param.geotiff_fn = '/cresis/projects/GIS_data/greenland/Landsat-7/Greenland_natural_150m.tif';
                break
              case 2;
                param.geotiff_fn = '/cresis/projects/GIS_data/antarctica/Landsat-7/Antarctica_LIMA_480m.tif';
                break
              case 3;
                param.geotiff_fn = '/cresis/projects/GIS_data/canada/Landsat-7/Canada_250m.tif';
                break
            end
          end
        else
          %FTP paths to the GeoTIFFs
          switch param.location_en
            case 1;
              param.geotiff_fn = '/data/picker/GIS_data/greenland/Landsat-7/Greenland_natural_150m.tif';
            case 2;
              param.geotiff_fn = '/data/picker/GIS_data/antarctica/Landsat-7/Antarctica_LIMA_480m.tif';
            case 3;
              param.geotiff_fn = '/data/picker/GIS_data/canada/Landsat-7/Canada_250m.tif';
          end
        end
      elseif param.image_en == 2
        fprintf('Enter the path/filename to your custom image below.\n');
        param.geotiff_fn = input('Enter Path with filename and extension: ','s');
        if exist(param.geotiff_fn,'file');
          break
        else
          fprintf('\nPath does not exist. Please enter a valid path\n');
        end
      else
        fprintf('\n------Invalid choice. Must be "1" or "2"------\n');
      end
      
      break
    end
    
  elseif param.runType == 1
    %Non-Interactive Input Mode (Use param)
    
    fprintf('=======================================================\n');
    fprintf('Interactive Geographic Search\n');
    fprintf('=======================================================\n');
    % Check params for errors.
    % Check location_en
    if param.location_en == 1
      fprintf('Your location of interest is Greenland.\n');
    elseif param.location_en == 2
      fprintf('Your location of interest is Antarctica.\n');
    elseif param.location_en == 3
      fprintf('Your location of interest is Canada.\n');
    else
      error('Invalid Location Entered. Check your input values.')
    end
    
    %Check RADARID Params
    if param.radarID == 1
      fprintf('Download for RDS Radars Supported\n');
    elseif param.radarID == 2
      fprintf('Download for Snow Radars NOT Supported\n');
    elseif param.radarID == 3
      fprintf('Download for Accum Radars NOT Supported\n');
    elseif param.radarID == 4
      fprintf('Download for KUBand Radars NOT Supported\n');
    else
      error('Invalid radarID Entered. Check your input values.');
    end
    
    
    if param.radarID == 1
      
      % Check download params
      if param.download_en == 1 || param.download_en == 0
      else
        error('Please enter "true" or "false" for download_en');
      end
      if param.download_en
        %Check download directory.
        param.dPath_check = exist(param.download_dir,'dir');
        if param.dPath_check == 0
          mkdir(param.download_dir);
        else
          param.download_dir = fullfile(param.download_dir,filesep);
        end
        %Check the output type
        if param.extID == 1 || param.extID == 2 || param.extID == 3
        else
          error('Invalid output type. Must be 1(CSV), 2(MAT), or 3(BOTH)');
        end
        %Check the download_type
        if param.download_tid == 1
          param.download_type = 'csv_good';
        elseif param.download_tid == 2
          param.download_type = 'csv';
        else
          error('Invalid download type. Must be (1) "csv_good" or  (2) "csv"');
        end
        %Set download filename
        param.download_out_fn = strcat(param.download_dir,param.file_name);
        % Check the downsampling input
        if param.down_sample == 1 || param.down_sample == 0
        else
          error('Please enter "true" or "false" for down_sample');
        end
        if ~isnumeric(param.dsamp_meters)
          error('dsamp_factor (N) must be a number.');
        end
      end
    end
    % Check GeoTiff params
    if param.image_en == 2
      if exist(param.geotiff_fn,'file');
      else
        error('Invalid path to GeoTiff, enter a valid path.');
      end
    elseif param.image_en == 1
      if param.dataGetType == 1
        if ispc
          switch param.location_en
            case 1;
              param.geotiff_fn = '\\titan\data3\GIS_data\greenland\Landsat-7\Greenland_natural_150m.tif';
            case 2;
              param.geotiff_fn = '\\titan\data3\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_peninsula.tif';
            case 3;
              param.geotiff_fn = '\\titan\data3\GIS_data\canada\Landsat-7\Canada_250m.tif';
          end
        else
          switch param.location_en
            case 1;
              param.geotiff_fn = '/cresis/data3/GIS_data/greenland/Landsat-7/Greenland_natural_150m.tif';
            case 2;
              param.geotiff_fn = '/cresis/data3/GIS_data/antarctica/Landsat-7/Antarctica_LIMA_peninsula.tif';
            case 3;
              param.geotiff_fn = '/cresis/data3/GIS_data/canada/Landsat-7/Canada_250m.tif';
          end
        end
      else
        %FTP paths to the GeoTIFFs
        switch param.location_en
          case 1;
            param.geotiff_fn = '/data/picker/GIS_data/greenland/Landsat-7/Greenland_natural_150m.tif';
          case 2;
            param.geotiff_fn = '/data/picker/GIS_data/antarctica/Landsat-7/Antarctica_LIMA_peninsula.tif';
          case 3;
            param.geotiff_fn = '/data/picker/GIS_data/canada/Landsat-7/Canada_250m.tif';
        end
      end
      
    else
      fprintf('param.runType is invalid. Should be (1 or 2)\n');
    end
    
  end
  
  fprintf('-----------------------------------------------------------------\n');
  fprintf('User Settings Complete. Please wait while the files load.\n');
  fprintf('-----------------------------------------------------------------\n');
  
  
  %% Set Data Paths
  
  % SPECIFY DATA LOCATIONS FOR LOCAL OR FTP DATA
  % MUST SEPERATE PC AND UNIX PATHS FOR LOCAL
  if param.radarID == 1
    if param.dataGetType == 1
      %Local Data Paths
      if ispc
        fns = {};
        seasons = {};
        if param.location_en == 1 || param.location_en == 3;
          % Greenland and NE Canada
          fns{end+1} = '\\titan\scratch2\mdce\icards\1993_Greenland_P3\CSARP_post\csv\Browse_1993_Greenland_P3.csv';
          seasons{end+1,1} = '1993_Greenland_P3';
          fns{end+1} = '\\titan\scratch2\mdce\icards\1995_Greenland_P3\CSARP_post\csv\Browse_1995_Greenland_P3.csv';
          seasons{end+1,1} = '1995_Greenland_P3';
          fns{end+1} = '\\titan\scratch2\mdce\icards\1996_Greenland_P3\CSARP_post\csv\Browse_1996_Greenland_P3.csv';
          seasons{end+1,1} = '1996_Greenland_P3';
          fns{end+1} = '\\titan\scratch2\mdce\icards\1997_Greenland_P3\CSARP_post\csv\Browse_1997_Greenland_P3.csv';
          seasons{end+1,1} = '1997_Greenland_P3';
          fns{end+1} = '\\titan\scratch2\mdce\icards\1998_Greenland_P3\CSARP_post\csv\Browse_1998_Greenland_P3.csv';
          seasons{end+1,1} = '1998_Greenland_P3';
          fns{end+1} = '\\titan\scratch2\mdce\icards\1999_Greenland_P3\CSARP_post\csv\Browse_1999_Greenland_P3.csv';
          seasons{end+1,1} = '1999_Greenland_P3';
          fns{end+1} = '\\titan\scratch2\mdce\icards\2001_Greenland_P3\CSARP_post\csv\Browse_2001_Greenland_P3.csv';
          seasons{end+1,1} = '2001_Greenland_P3';
          fns{end+1} = '\\titan\scratch2\mdce\icards\2002_Greenland_P3\CSARP_post\csv\Browse_2002_Greenland_P3.csv';
          seasons{end+1,1} = '2002_Greenland_P3';
          fns{end+1} = '\\titan\scratch2\mdce\acords\2003_Greenland_P3\CSARP_post\csv\Browse_2003_Greenland_P3.csv';
          seasons{end+1,1} = '2003_Greenland_P3';
          fns{end+1} = '\\titan\scratch2\mdce\acords\2005_Greenland_TO\CSARP_post\csv\Browse_2005_Greenland_TO.csv';
          seasons{end+1,1} = '2005_Greenland_TO';
          fns{end+1} = '\\titan\scratch2\mdce\mcrds\2006_Greenland_TO\CSARP_post\csv\Browse_2006_Greenland_TO.csv';
          seasons{end+1,1} = '2006_Greenland_TO';
          fns{end+1} = '\\titan\scratch2\mdce\mcrds\2007_Greenland_P3\CSARP_post\csv\Browse_2007_Greenland_P3.csv';
          seasons{end+1,1} = '2007_Greenland_P3';
          fns{end+1} = '\\titan\scratch2\mdce\mcrds\2008_Greenland_TO\CSARP_post\csv\Browse_2008_Greenland_TO.csv';
          seasons{end+1,1} = '2008_Greenland_TO';
          %fns{end+1} = '\\titan\scratch2\mdce\mcrds\2008_Greenland_Ground\CSARP_post\csv\Browse_2008_Greenland_Ground.csv';
          %seasons{end+1,1} = '2008_Greenland_Ground';
          fns{end+1} = '\\titan\scratch2\mdce\mcrds\2009_Greenland_TO\CSARP_post\csv\Browse_2009_Greenland_TO.csv';
          seasons{end+1,1} = '2009_Greenland_TO';
          fns{end+1} = '\\titan\scratch2\mdce\mcords\2010_Greenland_DC8\CSARP_post\csv\Browse_2010_Greenland_DC8.csv';
          seasons{end+1,1} = '2010_Greenland_DC8';
          fns{end+1} = '\\titan\scratch2\mdce\mcords\2010_Greenland_P3\CSARP_post\csv\Browse_2010_Greenland_P3.csv';
          seasons{end+1,1} = '2010_Greenland_P3';
          fns{end+1} = '\\titan\scratch2\mdce\mcords\2011_Greenland_TO\CSARP_post\csv\Browse_2011_Greenland_TO.csv';
          seasons{end+1,1} = '2011_Greenland_TO';
          fns{end+1} = '\\titan\scratch2\mdce\mcords2\2011_Greenland_P3\CSARP_post\csv\Browse_2011_Greenland_P3.csv';
          seasons{end+1,1} = '2011_Greenland_P3';
          fns{end+1} = '\\titan\scratch2\mdce\mcords2\2012_Greenland_P3\CSARP_post\csv\Browse_2012_Greenland_P3.csv';
          seasons{end+1,1} = '2012_Greenland_P3';
        else
          % Antarctica
          fns{end+1} = '\\titan\scratch2\mdce\icards\2002_Antarctica_P3chile\CSARP_post\csv\Browse_2002_Antarctica_P3chile.csv';
          seasons{end+1,1} = '2002_Antarctica_P3chile';
          fns{end+1} = '\\titan\scratch2\mdce\acords\2004_Antarctica_P3chile\CSARP_post\csv\Browse_2004_Antarctica_P3chile.csv';
          seasons{end+1,1} = '2004_Antarctica_P3chile';
          fns{end+1} = '\\titan\scratch2\mdce\mcords\2009_Antarctica_DC8\CSARP_post\csv\Browse_2009_Antarctica_DC8.csv';
          seasons{end+1,1} = '2009_Antarctica_DC8';
          fns{end+1} = '\\titan\scratch2\mdce\mcords\2009_Antarctica_TO\CSARP_post\csv\Browse_2009_Antarctica_TO.csv';
          seasons{end+1,1} = '2009_Antarctica_TO';
          fns{end+1} = '\\titan\scratch2\mdce\mcords\2010_Antarctica_DC8\CSARP_post\csv\Browse_2010_Antarctica_DC8.csv';
          seasons{end+1,1} = '2010_Antarctica_DC8';
          fns{end+1} = '\\titan\scratch2\mdce\mcords2\2011_Antarctica_TO\CSARP_post\csv\Browse_2011_Antarctica_TO.csv';
          seasons{end+1,1} = '2011_Antarctica_TO';
          fns{end+1} = '\\titan\scratch2\mdce\mcords\2011_Antarctica_DC8\CSARP_post\csv\Browse_2011_Antarctica_DC8.csv';
          seasons{end+1,1} = '2011_Antarctica_DC8';
        end
      else
        fns = {};
        seasons = {};
        if param.location_en == 1 || param.location_en == 3;
          % Greenland and NE Canada
          fns{end+1} = '/cresis/scratch2/mdce/icards/1993_Greenland_P3/CSARP_post/csv/Browse_1993_Greenland_P3.csv';
          seasons{end+1,1} = '1993_Greenland_P3';
          fns{end+1} = '/cresis/scratch2/mdce/icards/1995_Greenland_P3/CSARP_post/csv/Browse_1995_Greenland_P3.csv';
          seasons{end+1,1} = '1995_Greenland_P3';
          fns{end+1} = '/cresis/scratch2/mdce/icards/1996_Greenland_P3/CSARP_post/csv/Browse_1996_Greenland_P3.csv';
          seasons{end+1,1} = '1996_Greenland_P3';
          fns{end+1} = '/cresis/scratch2/mdce/icards/1997_Greenland_P3/CSARP_post/csv/Browse_1997_Greenland_P3.csv';
          seasons{end+1,1} = '1997_Greenland_P3';
          fns{end+1} = '/cresis/scratch2/mdce/icards/1998_Greenland_P3/CSARP_post/csv/Browse_1998_Greenland_P3.csv';
          seasons{end+1,1} = '1998_Greenland_P3';
          fns{end+1} = '/cresis/scratch2/mdce/icards/1999_Greenland_P3/CSARP_post/csv/Browse_1999_Greenland_P3.csv';
          seasons{end+1,1} = '1999_Greenland_P3';
          fns{end+1} = '/cresis/scratch2/mdce/icards/2001_Greenland_P3/CSARP_post/csv/Browse_2001_Greenland_P3.csv';
          seasons{end+1,1} = '2001_Greenland_P3';
          fns{end+1} = '/cresis/scratch2/mdce/icards/2002_Greenland_P3/CSARP_post/csv/Browse_2002_Greenland_P3.csv';
          seasons{end+1,1} = '2002_Greenland_P3';
          fns{end+1} = '/cresis/scratch2/mdce/acords/2003_Greenland_P3/CSARP_post/csv/Browse_2003_Greenland_P3.csv';
          seasons{end+1,1} = '2003_Greenland_P3';
          fns{end+1} = '/cresis/scratch2/mdce/acords/2005_Greenland_TO/CSARP_post/csv/Browse_2005_Greenland_TO.csv';
          seasons{end+1,1} = '2005_Greenland_TO';
          fns{end+1} = '/cresis/scratch2/mdce/mcrds/2006_Greenland_TO/CSARP_post/csv/Browse_2006_Greenland_TO.csv';
          seasons{end+1,1} = '2006_Greenland_TO';
          fns{end+1} = '/cresis/scratch2/mdce/mcrds/2007_Greenland_P3/CSARP_post/csv/Browse_2007_Greenland_P3.csv';
          seasons{end+1,1} = '2007_Greenland_P3';
          fns{end+1} = '/cresis/scratch2/mdce/mcrds/2008_Greenland_TO/CSARP_post/csv/Browse_2008_Greenland_TO.csv';
          seasons{end+1,1} = '2008_Greenland_TO';
          %fns{end+1} = '/cresis/scratch2/mdce/mcrds/2008_Greenland_Ground/CSARP_post/csv/Browse_2008_Greenland_Ground.csv';
          %seasons{end+1,1} = '2008_Greenland_Ground';
          fns{end+1} = '/cresis/scratch2/mdce/mcrds/2009_Greenland_TO/CSARP_post/csv/Browse_2009_Greenland_TO.csv';
          seasons{end+1,1} = '2009_Greenland_TO';
          fns{end+1} = '/cresis/scratch2/mdce/mcords/2010_Greenland_DC8/CSARP_post/csv/Browse_2010_Greenland_DC8.csv';
          seasons{end+1,1} = '2010_Greenland_DC8';
          fns{end+1} = '/cresis/scratch2/mdce/mcords/2010_Greenland_P3/CSARP_post/csv/Browse_2010_Greenland_P3.csv';
          seasons{end+1,1} = '2010_Greenland_P3';
          fns{end+1} = '/cresis/scratch2/mdce/mcords/2011_Greenland_TO/CSARP_post/csv/Browse_2011_Greenland_TO.csv';
          seasons{end+1,1} = '2011_Greenland_TO';
          fns{end+1} = '/cresis/scratch2/mdce/mcords2/2011_Greenland_P3/CSARP_post/csv/Browse_2011_Greenland_P3.csv';
          seasons{end+1,1} = '2011_Greenland_P3';
          fns{end+1} = '/cresis/scratch2/mdce/mcords2/2012_Greenland_P3/CSARP_post/csv/Browse_2012_Greenland_P3.csv';
          seasons{end+1,1} = '2012_Greenland_P3';
        else
          % Antarctica
          fns{end+1} = '/cresis/scratch2/mdce/icards/2002_Antarctica_P3chile/CSARP_post/csv/Browse_2002_Antarctica_P3chile.csv';
          seasons{end+1,1} = '2002_Antarctica_P3chile';
          fns{end+1} = '/cresis/scratch2/mdce/acords/2004_Antarctica_P3chile/CSARP_post/csv/Browse_2004_Antarctica_P3chile.csv';
          seasons{end+1,1} = '2004_Antarctica_P3chile';
          fns{end+1} = '/cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_post/csv/Browse_2009_Antarctica_DC8.csv';
          seasons{end+1,1} = '2009_Antarctica_DC8';
          fns{end+1} = '/cresis/scratch2/mdce/mcords/2009_Antarctica_TO/CSARP_post/csv/Browse_2009_Antarctica_TO.csv';
          seasons{end+1,1} = '2009_Antarctica_TO';
          fns{end+1} = '/cresis/scratch2/mdce/mcords/2010_Antarctica_DC8/CSARP_post/csv/Browse_2010_Antarctica_DC8.csv';
          seasons{end+1,1} = '2010_Antarctica_DC8';
          fns{end+1} = '/cresis/scratch2/mdce/mcords2/2011_Antarctica_TO/CSARP_post/csv/Browse_2011_Antarctica_TO.csv';
          seasons{end+1,1} = '2011_Antarctica_TO';
          fns{end+1} = '/cresis/scratch2/mdce/mcords/2011_Antarctica_DC8/CSARP_post/csv/Browse_2011_Antarctica_DC8.csv';
          seasons{end+1,1} = '2011_Antarctica_DC8';
        end
      end
      
    else
      %FTP Data Paths
      
      fns = {};
      seasons = {};
      if param.location_en == 1 || param.location_en == 3;
        % Greenland and NE Canada
        fns{end+1} = '/data/rds/1993_Greenland_P3/csv/Browse_1993_Greenland_P3.csv';
        seasons{end+1,1} = '1993_Greenland_P3';
        fns{end+1} = '/data/rds/1995_Greenland_P3/csv/Browse_1995_Greenland_P3.csv';
        seasons{end+1,1} = '1995_Greenland_P3';
        fns{end+1} = '/data/rds/1996_Greenland_P3/csv/Browse_1996_Greenland_P3.csv';
        seasons{end+1,1} = '1996_Greenland_P3';
        fns{end+1} = '/data/rds/1997_Greenland_P3/csv/Browse_1997_Greenland_P3.csv';
        seasons{end+1,1} = '1997_Greenland_P3';
        fns{end+1} = '/data/rds/1998_Greenland_P3/csv/Browse_1998_Greenland_P3.csv';
        seasons{end+1,1} = '1998_Greenland_P3';
        fns{end+1} = '/data/rds/1999_Greenland_P3/csv/Browse_1999_Greenland_P3.csv';
        seasons{end+1,1} = '1999_Greenland_P3';
        fns{end+1} = '/data/rds/2001_Greenland_P3/csv/Browse_2001_Greenland_P3.csv';
        seasons{end+1,1} = '2001_Greenland_P3';
        fns{end+1} = '/data/rds/2002_Greenland_P3/csv/Browse_2002_Greenland_P3.csv';
        seasons{end+1,1} = '2002_Greenland_P3';
        fns{end+1} = '/data/rds/2003_Greenland_P3/csv/Browse_2003_Greenland_P3.csv';
        seasons{end+1,1} = '2003_Greenland_P3';
        fns{end+1} = '/data/rds/2005_Greenland_TO/csv/Browse_2005_Greenland_TO.csv';
        seasons{end+1,1} = '2005_Greenland_TO';
        fns{end+1} = '/data/rds/2006_Greenland_TO/csv/Browse_2006_Greenland_TO.csv';
        seasons{end+1,1} = '2006_Greenland_TO';
        fns{end+1} = '/data/rds/2007_Greenland_P3/csv/Browse_2007_Greenland_P3.csv';
        seasons{end+1,1} = '2007_Greenland_P3';
        fns{end+1} = '/data/rds/2008_Greenland_TO/csv/Browse_2008_Greenland_TO.csv';
        seasons{end+1,1} = '2008_Greenland_TO';
        % fns{end+1} = '/data/rds/2008_Greenland_Ground/csv/Browse_2008_Greenland_Ground.csv';
        % seasons{end+1,1} = '2008_Greenland_Ground';
        fns{end+1} = '/data/rds/2009_Greenland_TO/csv/Browse_2009_Greenland_TO.csv';
        seasons{end+1,1} = '2009_Greenland_TO';
        fns{end+1} = '/data/rds/2010_Greenland_DC8/csv/Browse_2010_Greenland_DC8.csv';
        seasons{end+1,1} = '2010_Greenland_DC8';
        fns{end+1} = '/data/rds/2010_Greenland_P3/csv/Browse_2010_Greenland_P3.csv';
        seasons{end+1,1} = '2010_Greenland_P3';
        fns{end+1} = '/data/rds/2011_Greenland_P3/csv/Browse_2011_Greenland_P3.csv';
        seasons{end+1,1} = '2011_Greenland_P3';
        fns{end+1} = '/data/rds/2011_Greenland_TO/csv/Browse_2011_Greenland_TO.csv';
        seasons{end+1,1} = '2011_Greenland_TO';
      else
        % Antarctica
        fns{end+1} = '/data/rds/2002_Antarctica_P3chile/csv/Browse_2002_Antarctica_P3chile.csv';
        seasons{end+1,1} = '2002_Antarctica_P3chile';
        fns{end+1} = '/data/rds/2004_Antarctica_P3chile/csv/Browse_2004_Antarctica_P3chile.csv';
        seasons{end+1,1} = '2004_Antarctica_P3chile';
        fns{end+1} = '/data/rds/2009_Antarctica_DC8/csv/Browse_2009_Antarctica_DC8.csv';
        seasons{end+1,1} = '2009_Antarctica_DC8';
        fns{end+1} = '/data/rds/2009_Antarctica_TO/csv/Browse_2009_Antarctica_TO.csv';
        seasons{end+1,1} = '2009_Antarctica_TO';
        fns{end+1} = '/data/rds/2010_Antarctica_DC8/csv/Browse_2010_Antarctica_DC8.csv';
        seasons{end+1,1} = '2010_Antarctica_DC8';
        fns{end+1} = '/data/rds/2011_Antarctica_DC8/csv/Browse_2011_Antarctica_DC8.csv';
        seasons{end+1,1} = '2011_Antarctica_DC8';
      end
      
    end
  elseif param.radarID == 2 %Snow Radar
    if param.dataGetType == 1
      %Local Data
      if ispc
        %PC
        fns = {};
        seasons = {};
        if param.location_en == 1 || param.location_en == 3;
          %Greenland/Canada
          fns{end+1} = '\\titan\scratch2\mdce\snow\2011_Greenland_P3\CSARP_post\csv\Browse_2011_Greenland_P3.csv';
          seasons{end+1,1} = '2011_Greenland_P3';
        else
          %Antarctica
          %                 fns{end+1} = '';
          %                 seasons{end+1,1} = '';
        end
      else
        %UNIX
        fns = {};
        seasons = {};
        if param.location_en == 1 || param.location_en == 3;
          %Greenland/Canada
          fns{end+1} = '/cresis/scratch2/mdce/anow/2011_Greenland_P3/CSARP_post/csv/Browse_2011_Greenland_P3.csv';
          seasons{end+1,1} = '2011_Greenland_P3';
        else
          %Antarctica
          %                 fns{end+1} = '';
          %                 seasons{end+1,1} = '';
        end
      end
    else
      %FTP Data
      if param.location_en == 1 || param.location_en == 3;
        %Greenland/Canada
        fns{end+1} = '/data/snow/2011_Greenland_P3/csv/Browse_2011_Greenland_P3.csv';
        seasons{end+1,1} = '2011_Greenland_P3';
      else
        %Antarctica
        %             fns{end+1} = '';
        %             seasons{end+1,1} = '';
      end
    end
  elseif param.radarID == 3 %Accum Radar
    if param.dataGetType == 1
      %Local Data
      if ispc
        %PC
        fns = {};
        seasons = {};
        if param.location_en == 1 || param.location_en == 3;
          %Greenland/Canada
          fns{end+1} = '\\titan\scratch2\mdce\accum\2011_Greenland_P3\CSARP_post\csv\Browse_2011_Greenland_P3.csv';
          seasons{end+1,1} = '2011_Greenland_P3';
        else
          %Antarctica
          %                 fns{end+1} = '';
          %                 seasons{end+1,1} = '';
        end
      else
        %UNIX
        fns = {};
        seasons = {};
        if param.location_en == 1 || param.location_en == 3;
          %Greenland/Canada
          fns{end+1} = '/cresis/scratch2/mdce/accum/2011_Greenland_P3/CSARP_post/csv/Browse_2011_Greenland_P3.csv';
          seasons{end+1,1} = '2011_Greenland_P3';
        else
          %Antarctica
          %                 fns{end+1} = '';
          %                 seasons{end+1,1} = '';
        end
      end
    else
      %FTP Data
      if param.location_en == 1 || param.location_en == 3;
        %Greenland/Canada
        fns{end+1} = '/data/accum/2011_Greenland_P3/csv/Browse_2011_Greenland_P3.csv';
        seasons{end+1,1} = '2011_Greenland_P3';
      else
        %Antarctica
        %             fns{end+1} = '';
        %             seasons{end+1,1} = '';
      end
    end
  elseif param.radarID == 4 %KUBand Radar
    if param.dataGetType == 1
      %Local Data
      if ispc
        %PC
        fns = {};
        seasons = {};
        if param.location_en == 1 || param.location_en == 3;
          %Greenland/Canada
          fns{end+1} = '\\titan\scratch2\mdce\kuband\2011_Greenland_P3\CSARP_post\csv\Browse_2011_Greenland_P3.csv';
          seasons{end+1,1} = '2011_Greenland_P3';
        else
          %Antarctica
          %                 fns{end+1} = '';
          %                 seasons{end+1,1} = '';
        end
      else
        %UNIX
        fns = {};
        seasons = {};
        if param.location_en == 1 || param.location_en == 3;
          %Greenland/Canada
          fns{end+1} = '/cresis/scratch2/mdce/kuband/2011_Greenland_P3/CSARP_post/csv/Browse_2011_Greenland_P3.csv';
          seasons{end+1,1} = '2011_Greenland_P3';
        else
          %Antarctica
          %                 fns{end+1} = '';
          %                 seasons{end+1,1} = '';
        end
      end
    else
      %FTP Data
      if param.location_en == 1 || param.location_en == 3;
        %Greenland/Canada
        fns{end+1} = '/data/kuband/2011_Greenland_P3/csv/Browse_2011_Greenland_P3.csv';
        seasons{end+1,1} = '2011_Greenland_P3';
      else
        %Antarctica
        %             fns{end+1} = '';
        %             seasons{end+1,1} = '';
      end
    end
  end
  
  
  %% LOAD ALL OF THE BROWSE FILES (LOCAL OR FTP)
  
  if param.dataGetType == 1
    %Open Local Data
    
    %--------------------------------------------
    %      Open all of the BROWSE Files
    %--------------------------------------------
    fprintf('Opening sources\n');
    %Check if any sources exist for specified params.
    if isempty(fns)
      fprintf('========================================================\n');
      fprintf('                  NO DATA ERROR                         \n');
      fprintf('========================================================\n');
      fprintf('There is no data available for Geographic Search to Use.\n');
      fprintf('Contact cresis_data@cresis.ku.edu for more information. \n');
      fprintf('========================================================\n');
      error('No Available Data.');
    end
    lat_list = [];
    lon_list = [];
    frm_list = [];
    seasonID = [];
    for fn_idx = 1:length(fns)
      fn = fns{fn_idx};
      fprintf('%s: %s\n', datestr(now,'HH:MM:SS'), fn);
      %Open and parse the Brows CSV file based on radar type.
      data_ptr = fopen(fn,'r');
      
      if param.radarID == 1
        % Radar depth sounder format:
        C = textscan(data_ptr, '%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
        [LAT,LON,TIME,THICK,ELEVATION,FRAME,SURFACE,BOTTOM,QUALITY] = deal(C{:});
      else
        %Other radar format:
        fprintf('Data download not supported for Non-RDS Radars.\n');
        C = textscan(data_ptr, '%f%f%f%f%s','headerlines',1,'delimiter',',');
        [LAT,LON,ELEVATION,TIME,FRAME] = deal(C{:});
      end
      frm_list = cat(1,frm_list,FRAME);
      seasonID = cat(1,seasonID,ones(length(FRAME),1).*fn_idx);
      
      % Close data pointer
      if isa(data_ptr,'double')
        cid_sourceload = fclose(data_ptr);
        if cid_sourceload ~= 0
          fprint('Source Load Data Pointer failed to close.\n');
        end
      end
      lat_list = cat(1,lat_list,LAT);
      lon_list = cat(1,lon_list,LON);
    end
    
  else
    %Create an FTP Connection
    %Use MGET to download data
    
    %Create the temporary directory in the system temp folder.
    param.tPath = tempname(param.tPath);
    param.tmpPath = fullfile(param.tPath,filesep);
    mkdir(param.tmpPath);
    
    %Open the FTP Connection
    FTP = ftp('data.cresis.ku.edu');
    
    %--------------------------------------------
    %      Open all of the BROWSE Files
    %--------------------------------------------
    fprintf('Opening sources\n');
    lat_list = [];
    lon_list = [];
    frm_list = [];
    seasonID = [];
    for fn_idx = 1:length(fns)
      fn = fns{fn_idx};
      fprintf('%s: %s\n', datestr(now,'HH:MM:SS'), fn);
      
      %Get the file and save it to the temp path
      mget(FTP,fn,param.tmpPath);
      
      %Figure out the system type, and change path accordingly.
      if ~ispc
        %UNIX Path, Do Nothing.
      else
        %Windows Path, Switch from / to \
        fnr = findstr(fn,'/');
        fn(fnr) = '\';
      end
      
      %Open and parse the Brows CSV file based on radar type.
      fn = strcat(param.tmpPath,fns{fn_idx});
      data_ptr = fopen(fn,'r');
      
      if param.radarID == 1
        % Radar depth sounder format:
        C = textscan(data_ptr, '%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
        [LAT,LON,TIME,THICK,ELEVATION,FRAME,SURFACE,BOTTOM,QUALITY] = deal(C{:});
      else
        %Other radar format:
        fprintf('Data download not supported for Non-RDS Radars.\n');
        C = textscan(data_ptr, '%f%f%f%f%s','headerlines',1,'delimiter',',');
        [LAT,LON,ELEVATION,TIME,FRAME] = deal(C{:});
      end
      frm_list = cat(1,frm_list,FRAME);
      seasonID = cat(1,seasonID,ones(length(FRAME),1).*fn_idx);
      
      % Close data pointer
      if isa(data_ptr,'double')
        cid_sourceload = fclose(data_ptr);
        if cid_sourceload ~= 0
          fprint('Source Load Data Pointer failed to close.\n');
        end
      end
      lat_list = cat(1,lat_list,LAT);
      lon_list = cat(1,lon_list,LON);
    end
    
  end
  
  %% Plot Figures
  
  % PLOT THE FIRGURE (GeoTIFF and Flights)
  % (LOCAL OR FTP)
  
  if param.dataGetType == 1
    %Load local GeoTIFF and plot
    
    %--------------------------
    % Plotting the GeoTiff
    %--------------------------
    fprintf('-----------------------------------------------------------------\n');
    fprintf('Loading and Plotting GeoTIFF. Please wait ... ');
    proj = geotiffinfo(param.geotiff_fn);
    [RGB, R, tmp] = geotiffread(param.geotiff_fn);
    R = R/1e3;
    if ~exist('fig_h','var') || isempty(fig_h)
      fig_h = figure('name','Interactive Geographic Search Plot'); clf;
    else
      fig_h = figure('name','Interactive Geographic Search Plot');
      figure(fig_h); clf;
    end
    mapshow(RGB, R);
    if exist('lat','var') && ~isempty(lat)
      [X,Y] = projfwd(proj,lat,lon);
      X = X/1e3;
      Y = Y/1e3;
      hold on;
      plot(X,Y,varargin{:});
      plot(X(1),Y(1),'bo');
    end
    fprintf('%s', datestr(now,'HH:MM:SS\n'));
    
    %--------------------------
    % Plotting the flights
    %--------------------------
    fprintf('Loading Flightlines ... ');
    if exist('LAT','var') && ~isempty(LAT)
      [xF,yF] = projfwd(proj,lat_list,lon_list);
      xF = xF/1e3;
      yF = yF/1e3;
      hold on
      plot(xF,yF,'-r');
    end
    
    %--------------------------
    % Set Plot Params
    %--------------------------
    title('------ SEE COMMAND WINDOW FOR INSTRUCTIONS ------');
    xlabel('X (Km)');
    ylabel('Y (Km)');
    hold off;
    fprintf('%s', datestr(now,'HH:MM:SS\n'));
    
    %--------------------------
    % Print plot instructions
    %--------------------------
    fprintf('\n-------------------------------------------------------\n');
    fprintf('INSTRUCTIONS FOR PLOT: DRAWING YOUR STUDY AREA\n')
    fprintf('-------------------------------------------------------\n');
    
    fprintf('(1) Left click to draw a polygon around the area you want.\n');
    fprintf('(2) Double click once to end your polygon (Does not have to be closed)\n');
    fprintf('(3) You can then stretch, shrink or move the polygon.\n');
    fprintf('(4) Double click INSIDE the polygon to complete input.\n');
    
    
  else
    %MGET FTP GeoTIFF and plot
    
    fprintf('-----------------------------------------------------------------\n');
    %Download Geotiff
    fprintf('Downloading GeoTIFF ... ');
    mget(FTP,param.geotiff_fn,param.tmpPath);
    
    %Add tmpPath name to paramgeotiff_fn
    param.geotiff_fn = strcat(param.tmpPath,param.geotiff_fn);
    
    %Figure out the system type, and change path accordingly.
    if ~ispc
      %UNIX Path, Do Nothing.
    else
      %Windows Path, Switch from / to \
      fnr = findstr(param.geotiff_fn,'/');
      param.geotiff_fn(fnr) = '\';
    end
    fprintf('%s', datestr(now,'HH:MM:SS\n'));
    
    %Plot GeoTIFF
    fprintf('Plotting GeoTIFF ... ');
    proj = geotiffinfo(param.geotiff_fn);
    [RGB, R, tmp] = geotiffread(param.geotiff_fn);
    R = R/1e3;
    if ~exist('fig_h','var') || isempty(fig_h)
      fig_h = figure('name','Interactive Geographic Search Plot'); clf;
    else
      fig_h = figure('name','Interactive Geographic Search Plot');
      figure(fig_h); clf;
    end
    mapshow(RGB, R);
    if exist('lat','var') && ~isempty(lat)
      [X,Y] = projfwd(proj,lat,lon);
      X = X/1e3;
      Y = Y/1e3;
      hold on;
      plot(X,Y,varargin{:});
      plot(X(1),Y(1),'bo');
    end
    fprintf('%s', datestr(now,'HH:MM:SS\n'));
    
    %Plot Flightlines
    fprintf('Plotting Flights ... ');
    if exist('LAT','var') && ~isempty(LAT)
      [xF,yF] = projfwd(proj,lat_list,lon_list);
      xF = xF/1e3;
      yF = yF/1e3;
      hold on
      plot(xF,yF,'-r');
    end
    fprintf('%s', datestr(now,'HH:MM:SS\n'));
    
    %--------------------------
    % Set Plot Params
    %--------------------------
    title('------ SEE COMMAND WINDOW FOR INSTRUCTIONS ------');
    xlabel('X (Km)');
    ylabel('Y (Km)');
    hold off;
    
    %--------------------------
    % Print plot instructions
    %--------------------------
    fprintf('\n-------------------------------------------------------\n');
    fprintf('INSTRUCTIONS FOR PLOT: DRAWING YOUR STUDY AREA\n')
    fprintf('-------------------------------------------------------\n');
    
    fprintf('(1) Left click to draw a polygon around the area you want.\n');
    fprintf('(2) Double click once to end your polygon (Does not have to be closed)\n');
    fprintf('(3) You can then stretch, shrink or move the polygon.\n');
    fprintf('(4) Double click INSIDE the polygon to complete input.\n');
    
  end
  
  %% Get User Figure Input (Polygon)
  
  
  poly_handle = impoly;
  position = wait(poly_handle);
  [polyPts] = getPosition(poly_handle);
  close('Interactive Geographic Search Plot');
  pause(0.1);
  xPoly = polyPts(:,1);
  yPoly = polyPts(:,2);
  fprintf('-----------------------------------------------------------------\n');
  fprintf('User input to figure complete ... %s\n',datestr(now,'HH:MM:SS'));
  
  
  %% Create List of Good Segments
  
  
  fprintf('-----------------------------------------------------------------\n');
  fprintf('Creating files for output or download ... ');
  
  %Get a logical matrix of values within the polygon.
  IN = inpolygon(xF,yF,xPoly,yPoly);
  
  %Get only the good values for frame and seasonID based on 'IN'
  frm_list = frm_list(IN);
  seasonID = seasonID(IN);
  
  %Get the frames that are unique and the seasonID for those frames
  [frm_list gFrmIDX] = unique(frm_list);
  seasonID = seasonID(gFrmIDX);
  
  %Create a list of segments based on the good, unique frame list
  seg_list = {};
  for frm_idx = 1:length(frm_list)
    seg_list{end+1,1} = sprintf('%s_%s', frm_list{frm_idx}(1:8), ...
      frm_list{frm_idx}(9:10));
  end
  
  %Get on the good values for segment and the seasonID for those segments
  [seg_list gSegIDX] = unique(seg_list);
  seasonID = seasonID(gSegIDX);
  
  %Print complete time for 'File output creation'
  fprintf('%s', datestr(now,'HH:MM:SS\n'));
  
  
  %% Load/Download and Decimate Data
  if param.radarID == 1
    if param.dataGetType == 1
      %Load local data
      
      fprintf('-----------------------------------------------------------------\n');
      fprintf('Beginning download process. Please wait.\n')
      fprintf('-----------------------------------------------------------------\n');
      
      %Begin download and save process.
      if param.download_en
        Latitude = [];
        Longitude = [];
        UTC_time_sod = [];
        Thickness = [];
        Elevation = [];
        Frame = [];
        Surface = [];
        Bottom = [];
        Quality = [];
        Season = [];
        for seg_idx = 1:length(seg_list)
          fprintf('Downloading: %s ',seg_list{seg_idx});
          %Get the data path and create the fileID
          
          if ispc
            localDataPath = fns{seasonID(seg_idx)}(1:findstr(fns{seasonID(seg_idx)},'csv\')-1);
            fileID = sprintf('%s%s\\Data_%s.csv', ...
              localDataPath, param.download_type, seg_list{seg_idx});
            data_ptr = fopen(fileID,'r');
          else
            localDataPath = fns{seasonID(seg_idx)}(1:findstr(fns{seasonID(seg_idx)},'csv/')-1);
            fileID = sprintf('%s%s//Data_%s.csv', ...
              localDataPath, param.download_type, seg_list{seg_idx});
            data_ptr = fopen(fileID,'r');
          end
          if isempty(data_ptr)
            continue;
          end
          %Read in the data
          C = textscan(data_ptr,'%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
          [LAT,LON,TIME,THICK,ELEVATION,FRAME,SURFACE,BOTTOM,QUALITY] = deal(C{:});
          
          %Create list of seasons in downloaded dataset.
          SEASON = ones(length(LAT),1);
          SEASON = seasons(SEASON.*seasonID(seg_idx));
          
          %If no csv_good data exists at the last segment. End the script and
          %error out.
          no_dat_id = 1;
          if isempty(LAT)
            fprintf('---- Segment contains no data.\n')
            no_dat_id = 0;
          end
          if seg_idx == length(seg_list) && isempty(Latitude)
            if isa(data_ptr,'double')
              % This is a file id so we need to close it
              cidd = fclose(data_ptr);
            end
            error('No good data in your StudyArea.');
            break
          end
          
          
          %If data is good, continue.
          if isa(data_ptr,'double')
            % This is a file id so we need to close it
            cidd = fclose(data_ptr);
          end
          
          if ~param.down_sample
            fprintf('... %s\n',datestr(now,'HH:MM:SS'));
          end
          
          %Create logical of good indexes (Inside the study polygon)
          [xF,yF] = projfwd(proj,LAT,LON);
          xF = xF/1e3;   yF = yF/1e3;
          IN = inpolygon(xF,yF,xPoly,yPoly);
          
          % Keep only data in the study area (Use the IN from inpolygon)
          LAT=LAT(IN);      LON=LON(IN);         ELEVATION=ELEVATION(IN);
          TIME=TIME(IN);    FRAME=FRAME(IN);     QUALITY=QUALITY(IN);
          THICK=THICK(IN);  BOTTOM=BOTTOM(IN);   SURFACE=SURFACE(IN);
          SEASON=SEASON(IN);
          
          if isempty(LAT) && no_dat_id
            fprintf('---- No data inside study area.\n')
            no_dat_id = 0;
          end
          
          %Decimate data by dsamp_meters with low-pass filter
          if no_dat_id
            if param.down_sample
              fprintf('---- Downsampling: %s ... ',seg_list{seg_idx});
              
              along_track = geodetic_to_along_track(LAT,LON,ELEVATION-SURFACE);
              decim_idxs = get_equal_alongtrack_spacing_idxs(along_track,param.dsamp_meters);
              
              % Put the decimated, IN, data into the final variables.
              Latitude = cat(1,Latitude,LAT(decim_idxs));
              Longitude = cat(1,Longitude,LON(decim_idxs));
              UTC_time_sod = cat(1,UTC_time_sod,TIME(decim_idxs));
              Elevation = cat(1,Elevation,ELEVATION(decim_idxs));
              Frame = cat(1,Frame,FRAME(decim_idxs));
              Quality = cat(1,Quality,QUALITY(decim_idxs));
              Surface = cat(1,Surface,SURFACE(decim_idxs));
              Thickness = cat(1,Thickness,THICK(decim_idxs));
              Bottom = cat(1,Bottom,BOTTOM(decim_idxs));
              Season = cat(1,Season,SEASON(decim_idxs));
              fprintf('%s', datestr(now,'HH:MM:SS\n'));
            else
              % Put the un-decimated data into the final variables.
              Latitude = cat(1,Latitude,LAT);
              Longitude = cat(1,Longitude,LON);
              UTC_time_sod = cat(1,UTC_time_sod,TIME);
              Elevation = cat(1,Elevation,ELEVATION);
              Frame = cat(1,Frame,FRAME);
              Quality = cat(1,Quality,QUALITY);
              Surface = cat(1,Surface,SURFACE);
              Thickness = cat(1,Thickness,THICK);
              Bottom = cat(1,Bottom,BOTTOM);
              Season = cat(1,Season,SEASON);
            end
          end
          
        end
        fprintf('Download Complete... %s\n',datestr(now,'HH:MM:SS'));
        fprintf('-----------------------------------------------------------------\n');
      end
    else
      %Download FTP data
      
      fprintf('-----------------------------------------------------------------\n');
      fprintf('Beginning download process. Please wait.\n')
      fprintf('-----------------------------------------------------------------\n');
      
      if param.download_en
        Latitude = [];
        Longitude = [];
        UTC_time_sod = [];
        Thickness = [];
        Elevation = [];
        Frame = [];
        Surface = [];
        Bottom = [];
        Quality = [];
        Season = [];
        
        for seg_idx = 1:length(seg_list)
          fprintf('Downloading: %s ',seg_list{seg_idx});
          param.fileID = sprintf('/data/rds/%s/%s/Data_%s.csv', ...
            seasons{seasonID(seg_idx)}, param.download_type, seg_list{seg_idx});
          mget(FTP,param.fileID,param.tmpPath);
          %Figure out the system type, and change path accordingly.
          if ~ispc
            %UNIX Path, Do Nothing.
          else
            %Windows Path, Switch from / to \
            fnr = findstr(param.fileID,'/');
            param.fileID(fnr) = '\';
          end
          data_ptr = fopen(strcat(param.tmpPath,param.fileID),'r');
          
          if isempty(data_ptr)
            continue;
          end
          
          C = textscan(data_ptr,'%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
          [LAT,LON,TIME,THICK,ELEVATION,FRAME,SURFACE,BOTTOM,QUALITY] = deal(C{:});
          
          SEASON = ones(length(LAT),1);
          SEASON = seasons(SEASON.*seasonID(seg_idx));
          
          %If no csv_good data exists at the last segment. End the script and
          %error out.
          if isempty(LAT)
            fprintf('Segment contains no data.\n')
            if seg_idx == length(seg_list) && isempty(Latitude)
              if isa(data_ptr,'double')
                % This is a file id so we need to close it
                cidd = fclose(data_ptr);
              end
              error('No good data in your StudyArea.');
            end
            break
          end
          
          %If data is good, continue.
          if isa(data_ptr,'double')
            % This is a file id so we need to close it
            cidd = fclose(data_ptr);
          end
          
          if ~param.down_sample
            fprintf('... %s\n',datestr(now,'HH:MM:SS'));
          end
          
          %Create logical of good indexes (Inside the study polygon)
          [xF,yF] = projfwd(proj,LAT,LON);
          xF = xF/1e3;   yF = yF/1e3;
          IN = inpolygon(xF,yF,xPoly,yPoly);
          
          % Keep only data in the study area (Use the IN from inpolygon)
          LAT=LAT(IN);      LON=LON(IN);         ELEVATION=ELEVATION(IN);
          TIME=TIME(IN);    FRAME=FRAME(IN);     QUALITY=QUALITY(IN);
          THICK=THICK(IN);  BOTTOM=BOTTOM(IN);   SURFACE=SURFACE(IN);
          SEASON=SEASON(IN);
          
          %Decimate data by dsamp_meters with low-pass filter
          if param.down_sample
            fprintf('---- Downsampling: %s ... ',seg_list{seg_idx});
            
            % Get the indexes of along_track spaced points.
            along_track = geodetic_to_along_track(LAT,LON,ELEVATION-SURFACE);
            decim_idxs = get_equal_alongtrack_spacing_idxs(along_track,param.dsamp_meters);
            
            % Put the decimated, IN, data into the final variables.
            Latitude = cat(1,Latitude,LAT(decim_idxs));
            Longitude = cat(1,Longitude,LON(decim_idxs));
            UTC_time_sod = cat(1,UTC_time_sod,TIME(decim_idxs));
            Elevation = cat(1,Elevation,ELEVATION(decim_idxs));
            Frame = cat(1,Frame,FRAME(decim_idxs));
            Quality = cat(1,Quality,QUALITY(decim_idxs));
            Surface = cat(1,Surface,SURFACE(decim_idxs));
            Thickness = cat(1,Thickness,THICK(decim_idxs));
            Bottom = cat(1,Bottom,BOTTOM(decim_idxs));
            Season = cat(1,Season,SEASON(decim_idxs));
            fprintf('%s', datestr(now,'HH:MM:SS\n'));
          else
            % Put the un-decimated data into the final variables.
            Latitude = cat(1,Latitude,LAT);
            Longitude = cat(1,Longitude,LON);
            UTC_time_sod = cat(1,UTC_time_sod,TIME);
            Elevation = cat(1,Elevation,ELEVATION);
            Frame = cat(1,Frame,FRAME);
            Quality = cat(1,Quality,QUALITY);
            Surface = cat(1,Surface,SURFACE);
            Thickness = cat(1,Thickness,THICK);
            Bottom = cat(1,Bottom,BOTTOM);
            Season = cat(1,Season,SEASON);
          end
        end
        fprintf('Download Complete... %s\n',datestr(now,'HH:MM:SS'));
        fprintf('-----------------------------------------------------------------\n');
      end
    end
  end
  
  %% Save the output file/s
  
  if param.radarID == 1
    if param.download_en
      % Save the data into a file depending on the file extension supplied
      % in the user input.
      fprintf('Saving file %s ...',param.file_name);
      if param.extID == 2
        % Save to MAT -v6 file
        param.download_out_fn = strcat(param.download_dir,param.file_name,'.mat');
        save(param.download_out_fn,'-v6','Latitude','Longitude','UTC_time_sod','Thickness','Elevation','Frame','Surface','Bottom','Quality','Season');
      elseif param.extID == 1
        % Save to CSV file
        param.download_out_fn = strcat(param.download_dir,param.file_name,'.csv');
        fid_csv = fopen(param.download_out_fn,'w');
        fprintf(fid_csv,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
          'LAT','LON','UTCTime','THICK','ELEVATION','FRAME','SURFACE','BOTTOM','QUALITY','SEASON');
        for txt_idx = 1:length(Latitude)
          fprintf(fid_csv,'%2.6f,%2.6f,%5.4f,%6.2f,%4.4f,%s,%6.2f,%6.2f,%01d,%s,\n',...
            Latitude(txt_idx),Longitude(txt_idx),...
            UTC_time_sod(txt_idx),Thickness(txt_idx),...
            Elevation(txt_idx),Frame{txt_idx},Surface(txt_idx),Bottom(txt_idx),...
            Quality(txt_idx),Season{txt_idx});
        end
        fprintf('%s', datestr(now,'HH:MM:SS\n'));
        cidf = fclose(fid_csv);
        if cidf == 0
        else
          fprintf('CSV file release error.\n');
        end
      elseif param.extID == 3 %Save both a CSV and a MAT file
        %Save to CSV file
        param.download_out_fn = strcat(param.download_dir,param.file_name,'.csv');
        fid_csv = fopen(param.download_out_fn,'w');
        fprintf(fid_csv,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
          'LAT','LON','UTCTime','THICK','ELEVATION','FRAME','SURFACE','BOTTOM','QUALITY','SEASON');
        for txt_idx = 1:length(Latitude)
          fprintf(fid_csv,'%2.6f,%2.6f,%5.4f,%6.2f,%4.4f,%s,%6.2f,%6.2f,%01d,%s,\n',...
            Latitude(txt_idx),Longitude(txt_idx),...
            UTC_time_sod(txt_idx),Thickness(txt_idx),...
            Elevation(txt_idx),Frame{txt_idx},Surface(txt_idx),Bottom(txt_idx),...
            Quality(txt_idx),Season{txt_idx});
        end
        fprintf('%s', datestr(now,'HH:MM:SS\n'));
        cidf = fclose(fid_csv);
        if cidf == 0
        else
          fprintf('CSV file release error.\n');
        end
        %Save to MAT file
        clear param.download_out_fn;
        param.download_out_fn = strcat(param.download_dir,param.file_name,'.mat');
        save(param.download_out_fn,'-v6','Latitude','Longitude','UTC_time_sod','Thickness','Elevation','Frame','Surface','Bottom','Quality','Season');
      end
    end
  end
  
  
  %Print Segments and Frames
  fprintf('-----------------------------------------------------------------\n');
  fprintf('The Segments in your Study Area:\n');
  fprintf('-----------------------------------------------------------------\n');
  for seg_idx = 1:5:length(seg_list)
    if seg_idx+5 < length(seg_list)
      fprintf('%s | %s | %s | %s | %s |\n',seg_list{seg_idx},seg_list{seg_idx+1},seg_list{seg_idx+2},seg_list{seg_idx+3},seg_list{seg_idx+4});
    else
      if (length(seg_list) - seg_idx) == 1
        fprintf('%s |\n',seg_list{seg_idx});
      end
      if (length(seg_list) - seg_idx) == 2
        fprintf('%s | %s |\n',seg_list{seg_idx},seg_list{seg_idx+1});
      end
      if (length(seg_list) - seg_idx) == 3
        fprintf('%s | %s | %s |\n',seg_list{seg_idx},seg_list{seg_idx+1},seg_list{seg_idx+2});
      end
      if (length(seg_list) - seg_idx) == 4
        fprintf('%s | %s | %s | %s |\n',seg_list{seg_idx},seg_list{seg_idx+1},seg_list{seg_idx+2},seg_list{seg_idx+3});
      end
    end
  end
  
  fprintf('-----------------------------------------------------------------\n');
  fprintf('The Frames in your Study Area:\n');
  fprintf('-----------------------------------------------------------------\n');
  for frm_idx = 1:5:length(frm_list)
    if frm_idx+5 < length(frm_list)
      fprintf('%s | %s | %s | %s | %s |\n',frm_list{frm_idx},frm_list{frm_idx+1},frm_list{frm_idx+2},frm_list{frm_idx+3},frm_list{frm_idx+4});
    else
      if (length(frm_list) - frm_idx) == 1
        fprintf('%s |\n',frm_list{frm_idx});
      end
      if (length(frm_list) - frm_idx) == 2
        fprintf('%s | %s |\n',frm_list{frm_idx},frm_list{frm_idx+1});
      end
      if (length(frm_list) - frm_idx) == 3
        fprintf('%s | %s | %s |\n',frm_list{frm_idx},frm_list{frm_idx+1},frm_list{frm_idx+2});
      end
      if (length(frm_list) - frm_idx) == 4
        fprintf('%s | %s | %s | %s |\n',frm_list{frm_idx},frm_list{frm_idx+1},frm_list{frm_idx+2},frm_list{frm_idx+3});
      end
    end
  end
  
  fprintf('-----------------------------------------------------------------\n');
  disp('Get the files here: <a href="https://data.cresis.ku.edu/">CReSIS Data</a>');
  fprintf('-----------------------------------------------------------------\n');
  
  
  %% Final clean up
  
  if param.dataGetType == 2
    fprintf('Closing FTP Connection ...');
    close(FTP);
    fprintf('%s', datestr(now,'HH:MM:SS\n'));
    fprintf('Clearing Temporary Files ...');
    rmdir(param.tPath,'s');
    fprintf('%s', datestr(now,'HH:MM:SS\n'));
  end
  fprintf('Interactive Geographic Search Complete ... %s\n',datestr(now,'HH:MM:SS\n'));
  
catch ME
  if param.dataGetType == 2
    fprintf('----- MATLAB EXCEPTION -----\n');
    if isfield(param,'tpath')
      try
        fprintf('Clearing Temporary Files ...');
        rmdir(param.tPath,'s');
        fprintf('%s', datestr(now,'HH:MM:SS\n'));
      catch ME
        fprintf('Unable to clear temp files. Clear them manually.\n')
        rethrow(ME);
      end
    end
    rethrow(ME)
  else
    rethrow(ME)
  end
end
end
