function create_season_layerdata_files
% create_season_layerdata_files
%
% Creates a season layer file containing the lat, lon. The imb.picker loads
% this file when plotting flightlines without OPS. Loading all the
% individual CSARP_layerData files would be slow so this helps speed up the
% loading.
%
% Output filenames are of the form:
% /cresis/snfs1/dataproducts/csarp_support/layer/layer_accum_2018_Antarctica_TObas.mat
%
% Output file variables are all 1 by Nx vectors of the same length,
% segments are terminated with NaN:
%   lat: latitude in degrees
%   lon: longitude in degrees
%   frm: full frame ID 2019020401123
%   elev: elevation in meters
%   surf: surface two way travel time in seconds
%   bottom: bottom two way travel time in seconds
%   quality: integer enumeration of bottom quality, 1=good, 2=moderate,
%     3=poor or derived from another source
%
% Author: John Paden, Rohan Choudhari

%% Select season parameter files
param_fns = {};
% param_fns{end+1} = 'accum_param_2018_Antarctica_TObas.xls';
% param_fns{end+1} = 'rds_param_1993_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_1995_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_1996_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_1997_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_1998_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_1999_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2001_Greenland_P3.xls';
% % % param_fns{end+1} = 'rds_param_2002_Antarctica_P3chile.xls'; % May not work
% param_fns{end+1} = 'rds_param_2002_Greenland_P3.xls';
% % param_fns{end+1} = 'rds_param_2003_Greenland_P3.xls'; % May not work
% % param_fns{end+1} = 'rds_param_2004_Antarctica_P3chile.xls'; % May not work
% % param_fns{end+1} = 'rds_param_2005_Greenland_TO.xls'; % May not work
% param_fns{end+1} = 'rds_param_2006_Greenland_TO.xls';
% % param_fns{end+1} = 'rds_param_2007_Greenland_P3.xls'; % May not work
% param_fns{end+1} = 'rds_param_2008_Greenland_Ground.xls';
% param_fns{end+1} = 'rds_param_2008_Greenland_TO.xls';
% param_fns{end+1} = 'rds_param_2009_Antarctica_DC8.xls';
% param_fns{end+1} = 'rds_param_2009_Antarctica_TO.xls';
% param_fns{end+1} = 'rds_param_2009_Greenland_TO.xls';
% param_fns{end+1} = 'rds_param_2009_Greenland_TO_wise.xls';
% param_fns{end+1} = 'rds_param_2010_Antarctica_DC8.xls';
% param_fns{end+1} = 'rds_param_2010_Greenland_DC8.xls';
% param_fns{end+1} = 'rds_param_2010_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2010_Greenland_TO_wise.xls';
% param_fns{end+1} = 'rds_param_2011_Antarctica_DC8.xls';
% param_fns{end+1} = 'rds_param_2011_Antarctica_TO.xls';
% param_fns{end+1} = 'rds_param_2011_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2011_Greenland_TO.xls';
% param_fns{end+1} = 'rds_param_2012_Antarctica_DC8.xls';
% param_fns{end+1} = 'rds_param_2012_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2013_Antarctica_Basler.xls';
% param_fns{end+1} = 'rds_param_2013_Antarctica_P3.xls';
% param_fns{end+1} = 'rds_param_2013_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2014_Antarctica_DC8.xls';
% param_fns{end+1} = 'rds_param_2014_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2015_Greenland_C130.xls';
% param_fns{end+1} = 'rds_param_2015_Greenland_Polar6.xls';
% param_fns{end+1} = 'rds_param_2016_Antarctica_DC8.xls';
% param_fns{end+1} = 'rds_param_2016_Greenland_G1XB.xls';
% param_fns{end+1} = 'rds_param_2016_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2016_Greenland_Polar6.xls';
% param_fns{end+1} = 'rds_param_2016_Greenland_TOdtu.xls';
% param_fns{end+1} = 'rds_param_2017_Antarctica_Basler.xls';
% param_fns{end+1} = 'rds_param_2017_Antarctica_P3.xls';
% param_fns{end+1} = 'rds_param_2017_Antarctica_Polar6.xls';
% param_fns{end+1} = 'rds_param_2017_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2018_Antarctica_DC8.xls';
% param_fns{end+1} = 'rds_param_2018_Antarctica_Ground.xls';
% param_fns{end+1} = 'rds_param_2018_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2018_Greenland_Polar6.xls';
% param_fns{end+1} = 'rds_param_2019_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2019_Antarctica_GV.xls';
param_fns{end+1} = 'snow_param_2012_Greenland_P3.xls';

%% Setup layer load parameters
if 0
  layer_params = struct('name','surface');
  layer_params.source = 'layerdata';
  % layer_params(2).layerdata_source = 'CSARP_post/layerData';
  layer_params(2).name = 'bottom';
  layer_params(2).source = 'layerdata';
  % layer_params(2).layerdata_source = 'CSARP_post/layerData';
  layer_params(2).existence_check = false;
else
  % HACK!!!
  keyboard
  layer_params = struct('name','surface');
  layer_params.source = 'layerdata';
  layer_params(1).layerdata_source = 'layerData_koenig';
  layer_params(2).name = 'bottom';
  layer_params(2).source = 'layerdata';
  layer_params(2).layerdata_source = 'layerData_koenig';
  layer_params(2).existence_check = false;
  for lay_idx=2:30
    layer_params(lay_idx+1).name = sprintf('Koenig_%d',lay_idx);
    layer_params(lay_idx+1).source = 'layerdata';
    layer_params(lay_idx+1).layerdata_source = 'layerData_koenig';
    layer_params(lay_idx+1).existence_check = false;
  end
end

%% Loop to load each season
global gRadar;
for param_idx = 1:length(param_fns)
  
  % Initialize variables to be extracted from layers
  lat = [];
  lon = [];
  frm = [];
  surf = [];
  bottom = [];
  elev = [];
  quality = [];
  
  % Read in parameter spreadsheet
  param_fn = ct_filename_param(param_fns{param_idx});
  params = read_param_xls(param_fn,'');
  if 0
  params = ct_set_params(params,'cmd.generic',1);
  params = ct_set_params(params,'cmd.generic',0,'cmd.notes','do not process');
  else
    % HACK!!!
    keyboard
  params = ct_set_params(params,'cmd.generic',0);
  params = ct_set_params(params,'cmd.generic',1,'day_seg','20120330_04');
  params = ct_set_params(params,'debug',1);
  end
  disp(param_fns{param_idx})
  
  %% Load each segment
  for param_idx = 1:length(params)
    param = params(param_idx);
    param = merge_structs(param,gRadar);
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
      fprintf('%s\tdo not process\n', param.day_seg);
      continue;
    end
    
    % Reading layerData
    try
      layer = opsLoadLayers(param,layer_params);
    catch ME
      fprintf('%s\terror!!!\t%s\n', param.day_seg, ME.getReport);
      continue;
    end
    
    if isempty(layer)
      fprintf('%s\tno_layers_returned!!!\n', param.day_seg);
      continue;
    end;
    
    if isempty(layer(1).lat)
      fprintf('%s\tempty_layer!!!\n', param.day_seg);
      continue;
    end;
    
    if any(isnan(layer(1).twtt))
      fprintf('%s\tsurface_NaN!!!\n', param.day_seg);
    end;
    
    % Checking for inconsistent field lengths
    if length(layer(1).lat) ~= length(layer(1).lon) ...
      || length(layer(1).lat) ~= length(layer(1).frm) ...
      || length(layer(1).lat) ~= length(layer(1).elev) ...
      || length(layer(1).lat) ~= length(layer(1).twtt) ...
      || length(layer(1).lat) ~= length(layer(2).twtt) ...
      || length(layer(1).lat) ~= length(layer(2).quality)
      fprintf('%s\tmismatch_lengths\n');
      keyboard;
      continue;
    end
    
    % Concatenate data
    fprintf('%s\tloaded\t%d\t%d\n', param.day_seg, length(layer(1).lat), sum(~isnan(layer(2).twtt)));
    lat = [lat layer(1).lat NaN];
    lon = [lon layer(1).lon NaN];
    % Store full frame ID number 20190204_01_003 --> 2019020401003
    frm = [frm str2num(param.day_seg([1:8,10:11]))*1000+layer(1).frm NaN];
    elev = [elev layer(1).elev NaN];
    surf = [surf layer(1).twtt NaN];
    bottom = [bottom layer(2).twtt NaN];
    quality = [quality layer(2).quality NaN];
  end
  
  %% Save output
  out_fn_dir = ct_filename_support(param,'layer','');
  out_fn_name = sprintf('layer_%s_%s_%s.mat', param.post.ops.location, ct_output_dir(param.radar_name), param.season_name);
  out_fn = fullfile(out_fn_dir,out_fn_name);
  fprintf('  Saving %s\n\n', out_fn);
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  save(out_fn,'lat','lon','frm','elev','surf','bottom','quality');
end
