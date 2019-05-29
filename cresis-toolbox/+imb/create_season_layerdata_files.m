% Creates a season layerdata file containing the lat, lon
% The picker loads this file when plotting flightlines 

function create_season_layerdata_files

  global gRadar;

  %% Select Seasons
  season_names = {};
  %   season_names{end+1} = 'rds_param_1993_Greenland_P3';
  %   season_names{end+1} = 'rds_param_1995_Greenland_P3';
  %   season_names{end+1} = 'rds_param_1996_Greenland_P3';
  %   season_names{end+1} = 'rds_param_1997_Greenland_P3';
  %   season_names{end+1} = 'rds_param_1998_Greenland_P3';
  %   season_names{end+1} = 'rds_param_1999_Greenland_P3';
  %   season_names{end+1} = 'rds_param_2001_Greenland_P3';
  %   season_names{end+1} = 'rds_param_2002_Antarctica_P3chile';
  %   season_names{end+1} = 'rds_param_2002_Greenland_P3';
  %   season_names{end+1} = 'rds_param_2003_Greenland_P3';
  %   season_names{end+1} = 'rds_param_2004_Antarctica_P3chile';
  %   season_names{end+1} = 'rds_param_2005_Greenland_TO';
  %   season_names{end+1} = 'rds_param_2006_Greenland_TO';
  %   season_names{end+1} = 'rds_param_2008_Greenland_Ground';
  %   season_names{end+1} = 'rds_param_2008_Greenland_TO';
  %   season_names{end+1} = 'rds_param_2009_Antarctica_DC8';
  %   season_names{end+1} = 'rds_param_2009_Antarctica_TO';
  %   season_names{end+1} = 'rds_param_2009_Greenland_TO';
  %   season_names{end+1} = 'rds_param_2009_Greenland_TO_wise';
  %   season_names{end+1} = 'rds_param_2010_Antarctica_DC8';
  %   season_names{end+1} = 'rds_param_2010_Greenland_DC8';
  %   season_names{end+1} = 'rds_param_2010_Greenland_P3';
  %   season_names{end+1} = 'rds_param_2010_Greenland_TO_wise';
  %   season_names{end+1} = 'rds_param_2011_Antarctica_DC8';
  %   season_names{end+1} = 'rds_param_2011_Antarctica_TO';
  %   season_names{end+1} = 'rds_param_2011_Greenland_P3';
  %   season_names{end+1} = 'rds_param_2011_Greenland_TO';
  %   season_names{end+1} = 'rds_param_2012_Antarctica_DC8';
  %   season_names{end+1} = 'rds_param_2012_Greenland_P3';
  %   season_names{end+1} = 'rds_param_2013_Antarctica_Basler';
  %   season_names{end+1} = 'rds_param_2013_Antarctica_P3';
  %   season_names{end+1} = 'rds_param_2013_Greenland_P3';
  %   season_names{end+1} = 'rds_param_2014_Antarctica_DC8';
  %   season_names{end+1} = 'rds_param_2014_Greenland_P3';
  %   season_names{end+1} = 'rds_param_2015_Greenland_C130';
  %   season_names{end+1} = 'rds_param_2015_Greenland_Polar6';
  %   season_names{end+1} = 'rds_param_2016_Antarctica_DC8';
  %   season_names{end+1} = 'rds_param_2016_Greenland_G1XB';
  %   season_names{end+1} = 'rds_param_2016_Greenland_P3';
  %   season_names{end+1} = 'rds_param_2016_Greenland_Polar6';
  %   season_names{end+1} = 'rds_param_2016_Greenland_TOdtu';
  %   season_names{end+1} = 'rds_param_2017_Antarctica_Basler';
  %   season_names{end+1} = 'rds_param_2017_Antarctica_P3';
  %   season_names{end+1} = 'rds_param_2017_Antarctica_Polar6';
    season_names{end+1} = 'rds_param_2017_Greenland_P3';
  %   season_names{end+1} = 'rds_param_2018_Antarctica_DC8';
  %   season_names{end+1} = 'rds_param_2018_Antarctica_Ground';
  %   season_names{end+1} = 'rds_param_2018_Greenland_P3';
  %   season_names{end+1} = 'rds_param_2018_Greenland_Polar6';

  %% Select the destinaion where the files are to be saved
  out_dir = 'Y:\rohan\coverage_maps_loader\';

  if ~exist(out_dir,'dir')
    mkdir(out_dir);
  end

  %% Setting up params
  layer_params = [];
  idx = 1;
  layer_params(idx).name = 'bottom';
  layer_params(idx).source = 'layerData';
  % layer_params(idx).layerdata_source = 'CSARP_post/layerData';
  layer_params(idx).layerdata_source = 'layerData';

  %% Reading params, layerdata files, and storing lat, lon in the season layerdata file
  for season_idx = 1:length(season_names)

    lat = [];
    lon = [];  
    
    %Reading params from spreadsheet
    disp('\n\n')
    params = read_param_xls(ct_filename_param(strcat(season_names{season_idx},'.xls')),'','post');
    params = ct_set_params(params,'cmd.generic',1);
    params = ct_set_params(params,'cmd.generic',0,'cmd.notes','do not process');
    disp(season_names{season_idx})

    for param_idx = 1:length(params)
      param = params(param_idx);
      param = merge_structs(param,gRadar);

      %Reading layerData
      try
        layer = opsLoadLayers(param,layer_params);
      catch
        continue;
      end

      if isempty(layer)
        continue;
      end;

      %Looping through the layer data
      for layer_idx = 1:length(layer)

        %Checking for bad file
        if length(layer(layer_idx).lat) ~= length(layer(layer_idx).type)
          fprintf('    Bad file\n');
          continue;
        end
        lats = layer(layer_idx).lat;
        lons = layer(layer_idx).lon;
        lat = [lat lats];
        lon = [lon lons];
        
      end
    end
    save(char(strcat('./',season_names{season_idx},'.mat')), 'lat', 'lon');

  end
end