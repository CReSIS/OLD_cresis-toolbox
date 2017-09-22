% script tomo.run_collate
%
% Description: Run script for tomo.collate
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfData,
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

clear('param_override');

%% User Settings
param_override = [];  
param_override.sched.type = 'no scheduler';
param_override.sched.submit_arguments = '-l nodes=1:ppn=1,pmem=8000mb,walltime=120:00';

% example_setup = 'horizontal';
example_setup = 'vertical';

tomo_collate = [];
if strcmpi(example_setup,'vertical')
  %% Vertical multiwaveform fuse example
  params = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO.xls'),'20091224_01','post');
  params.cmd.frms = 16;
  params.cmd.generic = 1;
  
  % .in_dir: ct_filename_out directory to use at input, fused image will be stored here.
  tomo_collate.in_dir = 'nick_music';
  
  % .out_dir: ct_filename_out directory to which output surfData will be exported
  tomo_collate.out_dir = 'paden_surfData';
  
  % .geotiff_fn: DEM used to extract surface layer information
  tomo_collate.geotiff_fn = ct_filename_gis([],'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_surface.tif');
  tomo_collate.geotiff_bad_value = 32767;
  
  % .sv_cal_fn: filename containing steering vector calibration, leave empty to not use
  tomo_collate.sv_cal_fn = '';
  
  % .ice_mask_fn: filename of ice mask, leave empty to not use
  tomo_collate.ice_mask_fn = '';

  % .imgs: list of images II to use from .in_dir (Data_img_II*.mat). These
  %   should be listed from left most beam to right most beam when
  %   horizontal fuse is used. They should be listed from top to bottom
  %   when vertical fuse is used.
  tomo_collate.imgs = [1 2];
  tomo_collate.master_img_idx = 1;

  % .dem_guard: additional region around flight line to search for DEM points
  %   Setting too high slows the process down, setting too low will miss
  %   DEM points needed to properly represent the surface.
  tomo_collate.dem_guard = 16e3;
  
  % .dem_per_slice_guard: additional region around each slice to search for DEM points
  %   Setting too high slows the process down, setting too low will miss
  %   DEM points needed to properly represent the surface.
  tomo_collate.dem_per_slice_guard = 2500;
  
  % .data_threshold: pixel values above this will be clipped (set to this value)
  tomo_collate.data_threshold = 13.5;
  
  % middle of DOA axis
  tomo_collate.mid = -1;   % default
  % weight of smoothness enforcement
  tomo_collate.smooth_weight = -1;    % default
  % weight for 
  tomo_collate.smooth_var = -1;
  % expected twtt offset between DOA bins
  tomo_collate.smooth_slope = zeros(1,63);
  
  % .bounds_relative: DOA bins and along-track slices to trim off from each edge [top bottom left right]
  tomo_collate.bounds_relative = [3 2 0 0];
  % Number of refine/extract loops to run
  tomo_collate.max_loops = 50;
  
  tomo_collate.layer_source = 'ops';
  
  % .img_comb: Same as get_heights and combine worksheets. This is
  %   used for vertical using only. For N images,
  %   there will be (N-1) combines performed. Each combine is described by
  %   three numbers so that there should be (N-1)*3 elements in this
  %   vector. Each set of three numbers indicates the following:
  %     1st element: Minimum time to begin combine
  %     2nd element: Minimum time after ice surface to begin combine
  %     3rd element: Time at end of the preceeding waveform to not use
  tomo_collate.img_comb = [0 3.5e-5 1e-6];
  
  % .fuse: 'vertical' or 'horizontal' fuse
  tomo_collate.fuse_method = 'vertical';
  
  % .fuse_images_flag: runs fuse_images.m when true
  tomo_collate.fuse_images_flag = true;
  
  % .add_icemask_surfacedem_flag: runs add_icemask_surfacedem.m when true
  tomo_collate.add_icemask_surfacedem_flag = true;
  
  % create_surfData_flag: runs create_surfData.m when true
  tomo_collate.create_surfData_flag = true;

elseif strcmpi(example_setup,'horizontal')
  %% Horizontal multibeam fuse example
  params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'20140401_03','post');
  params.cmd.frms = 37;
  params.cmd.generic = 1;
  
  % .in_dir: ct_filename_out directory to use at input, fused image will be stored here.
  tomo_collate.in_dir = 'paden_music';
  
  % .out_dir: ct_filename_out directory to which output surfData will be exported
  tomo_collate.out_dir = 'paden_surfData';
  
  % .geotiff_fn: DEM used to extract surface layer information
  tomo_collate.geotiff_fn = ct_filename_gis([],fullfile('arctic','ArcticDEM','2014_Greenland_P3_20140401_03.tif'));
  tomo_collate.geotiff_bad_value = -32767;
  
  % .sv_cal_fn: filename containing steering vector calibration, leave empty to not use
  tomo_collate.sv_cal_fn = ct_filename_ct_tmp(rmfield(params(1),'day_seg'),'','sv_calibration','theta_cal.mat');
  
  % .ice_mask_fn: filename of ice mask, leave empty to not use
  tomo_collate.ice_mask_fn = ct_filename_gis([],'canada/ice_mask/03_rgi50_ArcticCanadaNorth/03_rgi50_ArcticCanadaNorth.mat');
  
  % .ocean_mask_fn: filename of ocean mask
  tomo_collate.ocean_mask_fn = ct_filename_gis([],fullfile('world','land_mask','Land_Mask_IDL_jharbeck','GSHHS_f_L1.shp'));

  % .imgs: list of images II to use from .in_dir (Data_img_II*.mat). These
  %   should be listed from left most beam to right most beam when
  %   horizontal fuse is used. They should be listed from top to bottom
  %   when vertical fuse is used.
  tomo_collate.imgs = [1 2 3];

  % .dem_guard: additional region in meters around flight line to search for DEM points
  %   Setting too high slows the process down, setting too low will miss
  %   DEM points needed to properly represent the surface.
  tomo_collate.dem_guard = 8e3;
  
  % .dem_per_slice_guard: additional region in meters around each slice to search for DEM points
  %   Setting too high slows the process down, setting too low will miss
  %   DEM points needed to properly represent the surface.
  tomo_collate.dem_per_slice_guard = 120;
  
  % .data_threshold: pixel values above this will be clipped (set to this value)
  tomo_collate.data_threshold = 13.5;
  
  % middle of DOA axis
  tomo_collate.mid = -1;   % default
  % weight of smoothness enforcement
  tomo_collate.smooth_weight = -1;    % default
  % weight for 
  tomo_collate.smooth_var = -1;
  % expected twtt offset between DOA bins
  tomo_collate.smooth_slope = zeros(1,63);
  
  % .bounds_relative: DOA bins and along-track slices to trim off from each edge [top bottom left right]
  tomo_collate.bounds_relative = [3 2 0 0];
  % Number of refine/extract loops to run
  tomo_collate.max_loops = 50;
  
  tomo_collate.layer_source = 'ops';
  
  % .img_comb: Same as get_heights and combine worksheets. This is
  %   used for vertical using only. For N images,
  %   there will be (N-1) combines performed. Each combine is described by
  %   three numbers so that there should be (N-1)*3 elements in this
  %   vector. Each set of three numbers indicates the following:
  %     1st element: Minimum time to begin combine
  %     2nd element: Minimum time after ice surface to begin combine
  %     3rd element: Time at end of the preceeding waveform to not use
  tomo_collate.img_comb = [0 3.5e-5 1e-6];
  
  % .fuse: 'vertical' or 'horizontal' fuse
  tomo_collate.fuse_method = 'horizontal';
  
  % .fuse_images_flag: runs fuse_images.m when true
  tomo_collate.fuse_images_flag = true;
  
  % .add_icemask_surfacedem_flag: runs add_icemask_surfacedem.m when true
  tomo_collate.add_icemask_surfacedem_flag = true;
  
  % create_surfData_flag: runs create_surfData.m when true
  tomo_collate.create_surfData_flag = true;
  
else
  error('%s is not a valid example_setup.', example_setup);
end

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

