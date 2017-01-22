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

%% 2014_Greenland_P3 20140401_03 User Settings

params = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'),'20161014_05','post');
params.cmd.generic = 1;
params.cmd.frms = 21;


tomo_collate = [];
% tomo_collate.in_dir = 'music_sravya';
tomo_collate.in_dir = 'music3D_jordan';
tomo_collate.out_dir = 'surfData';
% tomo_collate.ice_mask_fn = ct_filename_gis([],fullfile('canada','ice_mask','03_rgi50_ArcticCanadaNorth/03_rgi50_ArcticCanadaNorth.mat'));
% tomo_collate.geotiff_fn = ct_filename_gis([],fullfile('canada','DEM','SPI_All.tif'));
tomo_collate.geotiff_fn = ct_filename_gis([],fullfile('antarctica','DEM','BEDMAP2','original_data','bedmap2_tiff','bedmap2_surface.tif'));
tomo_collate.sv_cal_fn = ct_filename_ct_tmp(rmfield(params(1),'day_seg'),'','sv_calibration','theta_cal.mat');

% img_01: right looking (positive theta)
% img_02: nadir looking (theta == 0)
% img_03: left looking (negative theta)
tomo_collate.imgs = [1 2 3];
tomo_collate.master_img_idx = 2;

tomo_collate.fuse_images_flag = false;
  tomo_collate.vertical_fuse = false;
  % tomo_collate.vertical_fuse = true;
  tomo_collate.img_comb = [0 -inf 4e-6 -inf]; % only for vertical_fuse
tomo_collate.add_icemask_surfacedem_flag = false;
tomo_collate.create_surfData_flag = true;
tomo_collate.data_threshold = 13.5;

% middle of DOA axis
tomo_collate.mid = -1;   % default
% weight of smoothness enforcement
tomo_collate.smooth_weight = -1;    % default
% variance of smoothness enforcement
tomo_collate.smooth_var = -1;   % default
% expected twtt offset between DOA bins
tomo_collate.smooth_slope = zeros(1,63);

tomo_collate.layer_source = 'layerData';

%% Automated loading section 
% =========================================================================

global gRadar;

% Input checking
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

param_override.sched.submit_arguments = '-l nodes=1:ppn=1,pmem=8000mb,walltime=120:00';

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  param.tomo_collate = tomo_collate;
  tomo.collate(param,param_override);
end

