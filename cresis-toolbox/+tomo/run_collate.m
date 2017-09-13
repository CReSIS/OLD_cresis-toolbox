% script tomo.run_collate
%
% Description: Run script for tomo.collate
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfData,
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

clear('param_override');
param_override = [];  
param_override.sched.type = 'no scheduler';
param_override.sched.submit_arguments = '-l nodes=1:ppn=1,pmem=8000mb,walltime=120:00';

%% 2014_Greenland_P3 User Settings

% parameter spreadsheet of day or day_seg to be run
params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'20140325_07','post');
params.cmd.frms = [1 2];
params.cmd.generic = 1;

tomo_collate = [];

% directory from which radar slices are acquired
% fused images will also be saved in this location
tomo_collate.in_dir = 'paden_music';

% directory to which output layer data will be exported
tomo_collate.out_dir = 'paden_surfData';
% DEM used to extract surface layer information
tomo_collate.geotiff_fn = ct_filename_gis([],fullfile('arctic','ArcticDEM','2014_Greenland_P3_20140325.tif'));
tomo_collate.geotiff_bad_value = -32767;
% Ocean mask used to help fill DEM
% tomo_collate.ocean_mask_fn = ct_filename_tmp([],fullfile('ocean_mask','2014_Greenland_P3_20140506_01.tif'));

% calibrated steering vector
tomo_collate.sv_cal_fn = ct_filename_ct_tmp(rmfield(params(1),'day_seg'),'','sv_calibration','theta_cal.mat');

% ice mask fn
tomo_collate.ice_mask_fn = ct_filename_gis([],'canada/ice_mask/03_rgi50_ArcticCanadaNorth/03_rgi50_ArcticCanadaNorth.mat');

% img_01: right looking (positive theta)
% img_02: nadir looking (theta == 0)
% img_03: left looking (negative theta)
tomo_collate.imgs = [1 2 3];
tomo_collate.master_img_idx = 2;

% pixel intensity threshold, pixels above the data_threshold will be set to
% the value of data_threshold
tomo_collate.data_threshold = 13.5;

% middle of DOA axis
tomo_collate.mid = -1;   % default
% weight of smoothness enforcement
tomo_collate.smooth_weight = -1;    % default
% variance of smoothness enforcement
tomo_collate.smooth_var = -1;   % default
% expected twtt offset between DOA bins
tomo_collate.smooth_slope = zeros(1,63);
% DOA bins and along-track slices to trim off from each edge [top bottom left right]
tomo_collate.bounds_relative = [3 2 0 0];
% Number of refine/extract loops to run
tomo_collate.max_loops = 50;

tomo_collate.layer_source = 'layerData';

%% Collate Script Configuration

% when set to true, runs fuse_images.m
tomo_collate.fuse_images_flag = false;

  % when set to true, images will be fused vertically in fuse_images.m
  % rather than horizontally
  tomo_collate.vertical_fuse = false;
  
  % Holds pair parameters of image pair fusings
  % for each pair:
  %   1st element: Time after surface to begin fuse
  %   2nd element: Time after surface to end fuse
  %  tomo_collate.img_comb = [0 -inf 4e-6 -inf]; % only for vertical_fuse
  tomo_collate.img_comb = [0 3.5e-5]; % only for vertical_fuse

% when set to true, runs add_icemask_surfacedem.m
tomo_collate.add_icemask_surfacedem_flag = false;

% when set to true, runs create_surfData.m
tomo_collate.create_surfData_flag = true;

%% Automated loading section
% =========================================================================

global gRadar;

% Input checking
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
  
  param.tomo_collate = tomo_collate;
  tomo.collate(param,param_override);
end

