% script run_check_surface
%
% Runs check_surface.m
%
% cat /N/dcwan/projects/cresis/output/ct_tmp/check_surface/snow/2017_Greenland_P3/*.txt
%
% Author: John Paden

%% User Settings
param_override = [];

% params = read_param_xls(ct_filename_param('snow_param_2009_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2009_Antarctica_DC8.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_DC8.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2010_Antarctica_DC8.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2011_Antarctica_DC8.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2012_Antarctica_DC8.xls'));
% params = ct_set_params(params,'radar.wfs(1).nz_valid',[0 1]);

% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2015_Greenland_Polar6.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2016_Antarctica_DC8.xls'));
params = read_param_xls(ct_filename_param('snow_param_2017_Arctic_Polar5.xls'));
params = ct_set_params(params,'radar.nz_valid',[0 1 2 3]);

% params.cmd.generic=1;
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20170410_01');
% params = ct_set_params(params,'cmd.generic',1,'cmd.mission_names','^sea.*');
% params = ct_set_params(params,'cmd.generic',1,'cmd.mission_names','(?(?!^sea.*)^.*)');
% params = ct_set_params(params,'cmd.generic',1);
% params = ct_set_params(params,'cmd.generic',0,'cmd.notes','Do not process');

param_override.check_surface.debug_plots = {'visible','twtt','gps','nz'};

param_override.check_surface.save_records_en = false;

%% Automated Section
% =====================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Process each of the segments
ctrl_chain = {};
for param_idx = 1:length(params)
  param = params(param_idx);
  if isfield(param.cmd,'generic') && ~iscell(param.cmd.generic) && ~ischar(param.cmd.generic) && param.cmd.generic
    check_surface(param,param_override);
  end
end
