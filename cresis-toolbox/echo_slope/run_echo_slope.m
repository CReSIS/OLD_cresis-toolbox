% script run_echo_slope
%
% Script for running echo_slope
%
% Authors: Kevin Moore
%
% See also: echo_slope.m

%/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'));

% Syntax for running a specific segment and frame by overriding parameter spreadsheet values
%params = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO.xls'),'20091228_01');
%params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20140508_01');
params = ct_set_params(params,'cmd.frms',[57]);

%Set echo slope parameters

params = ct_set_params(params,'echo_slope.out_path','echo_slope');
params = ct_set_params(params,'echo_slope.in_path','standard');

%set rows and columns of the slop tiles
param_override.echo_slope.rows = 100;
param_override.echo_slope.cols = 100;


%set the maximum and minimum slope of the tiles in degrees
param_override.echo_slope.min_slope = -40;
param_override.echo_slope.max_slope = 40;


%set the number of tiles to create 
param_override.echo_slope.n = 11;

%set the sigma factor
param_override.echo_slope.sigma_factor = 2;


%param_override.sched.type = 'no scheduler'; % Example to override default cluster settings

%% Automated Section
% =====================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%echo_slope(params, param_override);

%Process each of the segments
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  echo_slope(param,param_override);
end

% Post process code (only include if necessary)