% script tomo.run_surface_extractor.m
%
% Description: Tracks surfaces of data slices. Calls tomo.collate.
%
% See also: tomo.collate
%
% Author: John Paden, Jordan Sprick, and Mingze Xu


% fn_dir: Directory where files are at

% params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'20140401_03','post');
% params.cmd.generic = 1;
% params.cmd.frms = 44;
params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'','post');

surf_extract = [];
surf_extract.out_dir = 'CSA_music';

surf_extract.ice_mask_fn = ct_filename_gis([],fullfile('canada','ice_mask','03_rgi50_ArcticCanadaNorth/03_rgi50_ArcticCanadaNorth.mat'));

surf_extract.geotiff_fn = ct_filename_gis([],fullfile('canada','DEM','SPI_All.tif'));

% img_01: right looking (positive theta)
% img_02: nadir looking (theta == 0)
% img_03: left looking (negative theta)
surf_extract.nadir_img = 2;

surf_extract.add_layers_flag = 1;
surf_extract.ice_twtt_flag = 1;
surf_extract.extract_flag = 1;
surf_extract.theta_calibrated = 1;

%% Automated loading section 
% =========================================================================

global gRadar;

clear('param_override');

% Input checking
if ~exist('params','var')
  error('Use run_master: A struct array of parameters must be passed in\n');
end
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  param.surf_extract = surf_extract;
  tomo.collate(param,param_override);
end
