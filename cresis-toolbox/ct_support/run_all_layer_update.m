% script run_all_layer_update
% run_all_layer_update
%
% Run layer_update on all seasons.
%   
% Author: John Paden
%
% See also: run_all_layer_update, layer_update

%% User Settings
% =========================================================================

%% User Settings: Select seasons
% -------------------------------------------------------------------------
param_fns = {};
% param_fns{end+1} = 'accum_param_2015_Antarctica_Ground.xls';
% param_fns{end+1} = 'accum_param_2018_Antarctica_TObas.xls';
% param_fns{end+1} = 'accum_param_2019_Antarctica_TObas.xls';
% param_fns{end+1} = 'rds_param_1993_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_1995_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_1996_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_1997_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_1998_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_1999_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2001_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2002_Antarctica_P3chile.xls'; % May not work
% param_fns{end+1} = 'rds_param_2002_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2003_Greenland_P3.xls'; % May not work
% param_fns{end+1} = 'rds_param_2004_Antarctica_P3chile.xls'; % May not work
% param_fns{end+1} = 'rds_param_2005_Greenland_TO.xls'; % May not work
% param_fns{end+1} = 'rds_param_2006_Greenland_TO.xls';
% param_fns{end+1} = 'rds_param_2007_Greenland_P3.xls'; % May not work
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
param_fns{end+1} = 'rds_param_2012_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2013_Antarctica_Basler.xls';
% param_fns{end+1} = 'rds_param_2013_Antarctica_P3.xls';
% param_fns{end+1} = 'rds_param_2013_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2014_Antarctica_DC8.xls';
% param_fns{end+1} = 'rds_param_2014_Greenland_P3.xls';
% param_fns{end+1} = 'rds_param_2015_Greenland_C130.xls';
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
% param_fns{end+1} = 'rds_param_2019_Antarctica_Ground.xls';
% param_fns{end+1} = 'rds_param_2019_Antarctica_GV.xls';
% param_fns{end+1} = 'snow_param_2012_Greenland_P3.xls';
% param_fns{end+1} = 'snow_param_2019_SouthDakota_CESSNA.xls';

param_override = [];
param_override.layer_update.in_path = 'layerData';
param_override.layer_update.out_path = 'layer';
param_override.layer_update.frames_records_en = true;

%% Automated Section
% =========================================================================
% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%% Loop to process each season
for param_idx = 1:length(param_fns)
  
  % Read in parameter spreadsheet
  param_fn = ct_filename_param(param_fns{param_idx});
  fprintf('Reading %s\n', param_fn);
  params = read_param_xls(param_fn,'');
  
  if isempty(params)
    continue;
  end
  
  % Run all segments (except "do not process")
  if 1
    params = ct_set_params(params,'cmd.generic',1);
    params = ct_set_params(params,'cmd.generic',0,'cmd.notes','do not process');
  else
    % HACK!!!
    keyboard
    params = ct_set_params(params,'cmd.generic',0);
    params = ct_set_params(params,'cmd.generic',1,'day_seg','20120330_04');
  end
  
  %% Update each segment
  for param_idx = 1:length(params)
    param = params(param_idx);
    
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
      fprintf('%s\tdo not process\n', param.day_seg);
      continue;
    end
    
    try
      % layer_update(param,param_override);
      layer_update
    catch ME
      fprintf('%s\terror!!!\t%s\n', param.day_seg, ME.getReport);
      continue;
    end
  end
end
