% Script run_post
%
% Loads the "post" worksheet from the parameter spreadsheet and then calls
% post with this information.
%
% Authors: Theresa Stumpf, John Paden, Reece Mathews
%
% See also: post.m

%% User Settings
param_override = [];

% params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'));
% params = read_param_xls(ct_filename_param('rds_param_2018_Antarctica_Ground.xls'));
params = read_param_xls(ct_filename_param('snow_param_2018_Greenland_P3.xls'));

% Example to run a specific segment and frame by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20180320_01');
params = ct_set_params(params,'cmd.frms',[3]);

% param_override.post.layers = struct('name', 'surface', 'source', 'layerData');

param_override.post.echo.elev_comp = 2;
param_override.post.echo.depth = '[min(Surface_Elev)-6 max(Surface_Elev)+10 ]';
param_override.post.echo.position = [1 1 1947 426];

if 0
  param_override.post.data_dirs = {'qlook_noise'};
  param_override.post.out_path = 'post_noise_final';
else
  param_override.post.data_dirs = {'standard'};
  param_override.post.out_path = 'post_standard_final';
end

params = ct_set_params(params,'post.maps_en',1);
params = ct_set_params(params,'post.echo_en',1);
params = ct_set_params(params,'post.layers_en',1);
params = ct_set_params(params,'post.data_en',0);
params = ct_set_params(params,'post.csv_en',0);
params = ct_set_params(params,'post.concat_en',0);
params = ct_set_params(params,'post.pdf_en',0);

%% Automated Section
% =====================================================================
% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%% Process each of the segments
% =====================================================================
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  post(param,param_override);
end

%% Create season-wide concatenated and browse files (CSV and KML)
% =====================================================================
concatenate_csv_kml = false;
for param_idx = 1:length(params)
  param = params(param_idx);
  cmd = param.cmd;
  if ct_generic_en(param)
    continue;
  end
  if param.post.concat_en
    concatenate_csv_kml = true;
    param = params(param_idx);
    break;
  end
end
if concatenate_csv_kml
  fprintf('Creating season concatenated and browse files (%s)\n', datestr(now));
  if ~isempty(param.post.out_path)
    post_path = ct_filename_out(param,param.post.out_path,'',1);
  else
    post_path = ct_filename_out(param,'post','',1);
  end
  % Create concatenated and browse files for all data
  csv_base_dir = fullfile(post_path,'csv');
  kml_base_dir = fullfile(post_path,'kml');
  if ~exist(kml_base_dir,'dir')
    mkdir(kml_base_dir)
  end
  run_concatenate_thickness_files(csv_base_dir,kml_base_dir,param);
  
  % Create concatenated and browse files for all data with thickness
  csv_base_dir = fullfile(post_path,'csv_good');
  kml_base_dir = fullfile(post_path,'kml_good');
  if ~exist(kml_base_dir,'dir')
    mkdir(kml_base_dir)
  end
  run_concatenate_thickness_files(csv_base_dir,kml_base_dir,param);
end
