% make_coverage_maps.m
%
% Used to create individual maps of each flight season's radar depth 
%  sounder (RDS) coverage. This script produces: 
%    1. individual JPEG maps (with Landsat GeoTIFF as background) 
%    2. text file with statistics on:
%       i.  percent time the radar was turned on
%       ii. percent of good radar data exists for the season.
%    3. a combined MATLAB figure (.fig) of all of the seasons.
%
% This script calls the individual scripts:
%  coverage_maps_old_without_gps
%  coverage_maps_old_with_gps
%  coverage_maps_with_param_file
%
% Authors: Steven Foga, John Paden
%
% See also: coverage_maps_old_without_gps, coverage_maps_old_with_gps,
%   coverage_maps_with_param_file, coverage_maps_master_fig_files

fprintf('===========================================================\n');
fprintf('coverage_maps\n\n');

global gRadar;

% ===================================================================
% User Settings
% ===================================================================

% location = 'Greenland';
% location = 'Canada';
location = 'Antarctica';

out_dir = '/cresis/scratch1/petan/coverage_maps/';
% out_dir = 'Z:\sfoga\coverage_maps\';

if strcmpi(location,'Greenland') || strcmpi(location,'Canada')
  season_names = {};
  season_names{end+1,1} = 'icards/1993_Greenland_P3';
%   season_names{end+1,1} = 'icards/1995_Greenland_P3';
%   season_names{end+1,1} = 'icards/1996_Greenland_P3';
%   season_names{end+1,1} = 'icards/1997_Greenland_P3';
%   season_names{end+1,1} = 'icards/1998_Greenland_P3';
%   season_names{end+1,1} = 'icards/1999_Greenland_P3';
%   season_names{end+1,1} = 'icards/2001_Greenland_P3';
%   season_names{end+1,1} = 'icards/2002_Greenland_P3';
%   season_names{end+1,1} = 'acords/2003_Greenland_P3';
%   season_names{end+1,1} = 'acords/2005_Greenland_TO';
%   season_names{end+1,1} = 'mcrds/2006_Greenland_TO';
%   season_names{end+1,1} = 'mcrds/2007_Greenland_P3';
%   season_names{end+1,1} = 'mcrds/2008_Greenland_Ground';
%   season_names{end+1,1} = 'mcrds/2008_Greenland_TO';
%   season_names{end+1,1} = 'mcrds/2009_Greenland_TO';
%   season_names{end+1,1} = 'mcords/2010_Greenland_DC8';
%   season_names{end+1,1} = 'mcords/2010_Greenland_P3';
%   season_names{end+1,1} = 'mcords/2011_Greenland_TO';
%   season_names{end+1,1} = 'mcords2/2011_Greenland_P3';
% season_names{end+1,1} = 'mcords2/2012_Greenland_P3';
elseif strcmpi(location,'Antarctica')
  season_names = {};
  season_names{end+1,1} = 'icards/2002_Antarctica_P3chile';
  season_names{end+1,1} = 'acords/2004_Antarctica_P3chile';
  season_names{end+1,1} = 'mcords/2009_Antarctica_DC8';
  season_names{end+1,1} = 'mcords/2009_Antarctica_TO';
  season_names{end+1,1} = 'mcords/2010_Antarctica_DC8';
  season_names{end+1,1} = 'mcords/2011_Antarctica_DC8';
  season_names{end+1,1} = 'mcords2/2011_Antarctica_TO';
  season_names{end+1,1} = 'mcords2/2012_Antarctica_DC8';
else
  return
end

% ===================================================================
% Automated Section
% ===================================================================
global gRadar

if strcmpi(location,'Greenland')
  
  geotiff_fns = {};
  geotiff_fns{1} = ct_filename_gis(gRadar,'greenland/Landsat-7/Greenland_natural_250m.tif');
  

  fig_size = [50 50 600 800];
  fig_paper_position = [0.25 2.5 6 8];
  fig_legend_size = [0.5742 0.1341 0.3051 0.0992];

elseif strcmpi(location,'Canada')

  geotiff_fns = {};
  geotiff_fns{1} = ct_filename_gis(gRadar,'canada/Landsat-7/Canada_250m.tif');
  
  fig_size = [50 50 560 420];
  fig_paper_position = [0.25 2.5 6 8];
  fig_legend_size = [100 25 0 0];
  
elseif strcmpi(location,'Antarctica')
  
  geotiff_fns = {};
  geotiff_fns{1} = ct_filename_gis(gRadar,'antarctica/Landsat-7/Antarctica_LIMA_480m.tif');  

  fig_size = [50 50 600 800];
  fig_legend_size = [235 0 0 0];
else
  return
end
 
    
    for season_idx = 1:length(season_names)
      fclose all;
      fprintf('Processing season %s\n', season_names{season_idx});
      
      filesep_idx = find(season_names{season_idx}=='/',1);
      clear param;
      param.radar_name = season_names{season_idx}(1:filesep_idx-1);
      param.season_name = season_names{season_idx}(filesep_idx+1:end);
      param.rdr_season_type = ct_filename_out(param,'','CSARP_post/csv/',true);
      param.rdr_season_good = ct_filename_out(param,'','CSARP_post/csv_good/',true);
      param.all_season_type = strcat(param.rdr_season_type,'Browse_',param.season_name,'.csv');
      
      % Replace stat_good_gre with param.stat_good
      param.stat_good = strcat(param.rdr_season_good,'Browse_',param.season_name,'.csv');
      
      % Replace stat_all_gre with param.stat_all
      param.stat_all = strcat(param.rdr_season_type,'Browse_',param.season_name,'.csv');
      
      
      % Navigation to 'params' toolbox (only needed to find param
      % spreadsheets)
      param.param_path = strrep(gRadar.path,'cresis-toolbox/','');
      param.param_path = strcat(param.param_path,'params/');
      
      param.gps_all = [];
      
      % GPS or Param file paths 
      if strcmp(param.season_name,'2004_Antarctica_P3chile')
        param.gps_all = '/cresis/scratch1/mdce/csarp_support/gps/ACORDS_2004_Antarctica_POS_GPS';
      elseif strcmp(param.season_name,'2003_Greenland_P3') 
        param.gps_all = '/cresis/scratch1/mdce/csarp_support/gps/ACORDS_2003_Greenland_POS_GPS';
      elseif strcmp(param.season_name,'2005_Greenland_TO') 
        param.gps_all = '/cresis/scratch1/mdce/csarp_support/gps/ACORDS_2005_Greenland_POS_GPS';
      elseif strcmp(param.season_name,'2006_Greenland_TO') 
        param.gps_all = '/cresis/scratch1/mdce/csarp_support/gps/MCRDS_2006_Greenland_POS_GPS';
      elseif strcmp(param.season_name,'2007_Greenland_P3')
        param.gps_all = '/cresis/scratch1/mdce/csarp_support/gps/MCRDS_2007_Greenland_POS_GPS';
      elseif strcmp(param.season_name,'2008_Greenland_TO')
        param.gps_all = '/cresis/scratch1/mdce/csarp_support/gps/MCRDS_2008_Greenland_POS_DGPS';
      elseif strcmp(param.season_name,'2008_Greenland_Ground')
        param.gps_all = '/cresis/scratch1/mdce/csarp_support/gps/MCRDS_2008_Greenland_Ground_GPS';
      elseif strcmp(param.season_name,'2009_Greenland_TO')
        param.gps_all = '/cresis/scratch1/mdce/csarp_support/gps/MCRDS_2009_Greenland_POS_GPS';
      elseif strcmp(param.season_name,'2009_Antarctica_DC8')
        param.param_path = strcat(param.param_path,'mcords_param_',param.season_name,'.xls');  
        param.segs = make_segment_list(param.param_path);
      elseif strcmp(param.season_name,'2009_Antarctica_TO')
        param.param_path = strcat(param.param_path,'mcords_param_',param.season_name,'.xls');  
        param.segs = make_segment_list(param.param_path);
      elseif strcmp(param.season_name,'2010_Antarctica_DC8')
        param.param_path = strcat(param.param_path,'mcords_param_',param.season_name,'.xls');  
        param.segs = make_segment_list(param.param_path); 
      elseif strcmp(param.season_name,'2011_Antarctica_DC8')
        param.param_path = strcat(param.param_path,'mcords_param_',param.season_name,'.xls');  
        param.segs = make_segment_list(param.param_path); 
      elseif strcmp(param.season_name,'2011_Antarctica_TO')  
        param.param_path = strcat(param.param_path,'mcords_param_',param.season_name,'.xls');  
        param.segs = make_segment_list(param.param_path); 
      elseif strcmp(param.season_name,'2010_Greenland_P3')
        param.param_path = strcat(param.param_path,'mcords_param_',param.season_name,'.xls');  
        param.segs = make_segment_list(param.param_path);
      elseif strcmp(param.season_name,'2010_Greenland_DC8')
        param.param_path = strcat(param.param_path,'mcords_param_',param.season_name,'.xls');  
        param.segs = make_segment_list(param.param_path);
      elseif strcmp(param.season_name,'2011_Greenland_TO')
        param.param_path = strcat(param.param_path,'mcords_param_',param.season_name,'.xls');  
        param.segs = make_segment_list(param.param_path);
      elseif strcmp(param.season_name,'2011_Greenland_P3')
        param.param_path = strcat(param.param_path,'mcords_param_',param.season_name,'.xls');  
        param.segs = make_segment_list(param.param_path);
      elseif strcmp(param.season_name,'2012_Greenland_P3')
        param.param_path = strcat(param.param_path,'mcords_param_','_param_',param.season_name,'.xls');  
        param.segs = make_segment_list(param.param_path);
      elseif strcmp(param.season_name,'2012_Antarctica_DC8')
        param.param_path = strcat(param.param_path,'mcords_param_','_param_',param.season_name,'.xls');  
        param.segs = make_segment_list(param.param_path);
      else
          fprintf('ERROR: No GPS files exist (2002 or before), or no matching season found.\n\n');
      end
      
      % No ct_filename exists for scratch1 (GPS file location), so
      % filepaths conversion must be performed manually.
      if ispc && ~isempty(param.gps_all)
          param.gps_all = strrep(param.gps_all,'/cresis/scratch1/','Y:\');
          param.gps_all = strrep(param.gps_all,'/','\');
      else
      end
        
      if strcmp(param.season_name,'2002_Antarctica_P3chile') || strcmp(param.season_name,'1993_Greenland_P3') || ...
          strcmp(param.season_name,'1995_Greenland_P3') || strcmp(param.season_name,'1996_Greenland_P3') || ...
          strcmp(param.season_name,'1997_Greenland_P3') || strcmp(param.season_name,'1998_Greenland_P3') || ...
          strcmp(param.season_name,'1999_Greenland_P3') || strcmp(param.season_name,'2001_Greenland_P3') || ...
          strcmp(param.season_name,'2002_Greenland_P3') || strcmp(param.season_name,'2002 Antarctica P3chile')
        
        coverage_maps_old_without_gps;
      
      elseif strcmp(param.season_name,'2003_Greenland_P3') || strcmp(param.season_name,'2005_Greenland_TO') ...
          || strcmp(param.season_name,'2006_Greenland_TO') || strcmp(param.season_name,'2007_Greenland_P3') ...
          || strcmp(param.season_name,'2008_Greenland_Ground') || strcmp(param.season_name,'2008_Greenland_TO') ...
          || strcmp(param.season_name,'2009_Greenland_TO') || strcmp(param.season_name,'2004_Antarctica_P3chile')
        
        coverage_maps_old_with_gps;
       
      else % All missions flown after 2009_Greenland_TO
        coverage_maps_with_param_file;
      end
    end
    return;