% script tomo.run_collate
%
% Description: Run script for tomo.collate
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfData
%
% Author: John Paden, Jordan Sprick, Mingze Xu, and Victor Berger

%% User Setup
% =====================================================================
param_override = [];

% example_setup = 'horizontal';
% example_setup = 'vertical';
example_setup = 'grid';

tomo_collate = [];
if strcmpi(example_setup,'vertical')
  %% Vertical multiwaveform fuse example
  params = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO.xls'),'','post');
  params = ct_set_params(params,'cmd.generic',0);
  params = ct_set_params(params,'cmd.generic',1,'day_seg','20091224_01|20091228|20091229|20100101_02|20100103|20100104_01|20100112');
  params = ct_set_params(params,'cmd.frms',[]);
  
  % .in_path: ct_filename_out directory to use at input, fused image will be stored here.
  tomo_collate.in_path = 'music3D';
  
  % .out_path: ct_filename_out directory to which output surfData will be exported
  tomo_collate.out_path = 'surfData';

  % .imgs: list of images II to use from .in_path (Data_img_II*.mat). These
  %   should be listed from left most beam to right most beam when
  %   horizontal fuse is used. They should be listed from top to bottom
  %   when vertical fuse is used.
  tomo_collate.imgs = {1,2};
  tomo_collate.master_img_idx = 1;
  
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
  
  % .geotiff_fn: DEM used to extract surface layer information
  tomo_collate.geotiff_fn = ct_filename_gis([],'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_surface.tif');
  tomo_collate.geotiff_bad_value = 32767;
  
  % .sv_cal_fn: filename containing steering vector calibration, leave empty to not use
  tomo_collate.sv_cal_fn = '';
  
  % .ice_mask_fn: filename of ice mask, leave empty to not use
  tomo_collate.ice_mask_fn = '';

  % .dem_guard: additional region around flight line to search for DEM points
  %   Setting too high slows the process down, setting too low will miss
  %   DEM points needed to properly represent the surface.
  tomo_collate.dem_guard = 16e3;
  
  % .dem_per_slice_guard: additional region around each slice to search for DEM points
  %   Setting too high slows the process down, setting too low will miss
  %   DEM points needed to properly represent the surface.
  tomo_collate.dem_per_slice_guard = 2500;
  
  % .bounds_relative: DOA bins and along-track slices to trim off from each edge [top bottom left right]
  tomo_collate.bounds_relative = [3 2 0 0];
  
  % .layer_params: parameter structure for opsLoadLayer (first layer should
  %   be ice top and second layer should be ice bottom)
  tomo_collate.layer_params = struct('name','surface','source','layerdata');
  tomo_collate.layer_params(2).name = 'bottom';
  tomo_collate.layer_params(2).source = 'layerdata';
  
  % surfData_mode: surfData mode ('overwrite','fillgaps', or 'append', note that append with the
  %   same surface name as an existing surface will overwrite that surface whereas fillgaps
  %   will leave the surface untouched if it already exists)
  tomo_collate.surfData_mode = 'append';
  
  % surfdata_cmds: surfdata commands to run
  tomo_collate.surfdata_cmds = [];
  tomo_collate.surfdata_cmds(end+1).cmd = 'trws';
  tomo_collate.surfdata_cmds(end).surf_names = {'bottom trws','bottom'};
  tomo_collate.surfdata_cmds(end).visible = true;
  
  % .fuse_images_flag: runs fuse_images.m when true
  tomo_collate.fuse_images_flag = true;
  
  % .add_icemask_surfacedem_flag: runs add_icemask_surfacedem.m when true
  tomo_collate.add_icemask_surfacedem_flag = true;
  
  % create_surfData_flag: runs create_surfData.m when true
  tomo_collate.create_surfData_flag = true;
  
  param_override.tomo_collate = tomo_collate;

elseif strcmpi(example_setup,'horizontal')
  %% Horizontal multibeam fuse example
  params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'','post');
  params = ct_set_params(params,'cmd.generic',0);
%   params = ct_set_params(params,'cmd.generic',1,'day_seg','20140325_05|20140325_06|20140325_07|20140401_03|20140506_01');
  params = ct_set_params(params,'cmd.generic',1,'day_seg','20140325_05');
  params = ct_set_params(params,'cmd.frms',[2]);

  % .in_path: ct_filename_out directory to use at input, fused image will be stored here.
  tomo_collate.in_path = 'music3D';
  
  % .out_path: ct_filename_out directory to which output surfData will be exported
  tomo_collate.out_path = 'surfData_no_MC';

  % .imgs: list of images II to use from .in_path (Data_img_II*.mat). These
  %   should be listed from left most beam to right most beam when
  %   horizontal fuse is used. They should be listed from top to bottom
  %   when vertical fuse is used.
  tomo_collate.imgs = {[1 2 3]};

  % .geotiff_fn: DEM used to extract surface layer information
  tomo_collate.geotiff_fn = ct_filename_gis([],fullfile('arctic','ArcticDEM','2014_Greenland_P3_20140401_03.tif'));
  tomo_collate.geotiff_bad_value = -32767;
  
  % .sv_cal_fn: filename containing steering vector calibration, leave empty to not use
  tomo_collate.sv_cal_fn = ct_filename_ct_tmp(rmfield(params(1),'day_seg'),'','sv_calibration','theta_cal.mat');
  
  % .ice_mask_fn: filename of ice mask, leave empty to not use
  tomo_collate.ice_mask_fn = ct_filename_gis([],'canada/ice_mask/03_rgi50_ArcticCanadaNorth/03_rgi50_ArcticCanadaNorth.mat');
  
  % .ocean_mask_fn: filename of ocean mask
  tomo_collate.ocean_mask_fn = ct_filename_gis([],fullfile('world','land_mask','Land_Mask_IDL_jharbeck','GSHHS_f_L1.shp'));

  % .dem_guard: additional region in meters around flight line to search for DEM points
  %   Setting too high slows the process down, setting too low will miss
  %   DEM points needed to properly represent the surface.
  tomo_collate.dem_guard = 12e3;
  
  % .dem_per_slice_guard: additional region in meters around each slice to search for DEM points
  %   Setting too high slows the process down, setting too low will miss
  %   DEM points needed to properly represent the surface.
  tomo_collate.dem_per_slice_guard = 240;
  
  % .bounds_relative: DOA bins and along-track slices to trim off from each edge [top bottom left right]
  tomo_collate.bounds_relative = [3 2 0 0];

  % .layer_params: parameter structure for opsLoadLayer (first layer should
  %   be ice top and second layer should be ice bottom)
  tomo_collate.layer_params = struct('name','surface','source','layerdata');
  tomo_collate.layer_params(2).name = 'bottom';
  tomo_collate.layer_params(2).source = 'layerdata';

  % surfData_mode: surfData mode ('overwrite','fillgaps', or 'append', note that append with the
  %   same surface name as an existing surface will overwrite that surface whereas fillgaps
  %   will leave the surface untouched if it already exists)
  tomo_collate.surfData_mode = 'append';
  
  % surfdata_cmds: surfdata commands to run
  tomo_collate.surfdata_cmds = [];
  
  tomo_collate.surfdata_cmds(end+1).cmd = 'detect';
  tomo_collate.surfdata_cmds(end).surf_names = 'bottom detect';
  tomo_collate.surfdata_cmds(end).visible = false;
  tomo_collate.surfdata_cmds(end).data_threshold = 13.5;
  
  tomo_collate.surfdata_cmds(end+1).cmd = 'extract';
  tomo_collate.surfdata_cmds(end).surf_names = 'bottom extract';
  tomo_collate.surfdata_cmds(end).visible = false;
  tomo_collate.surfdata_cmds(end).data_threshold = 13.5;
  
  tomo_collate.surfdata_cmds(end+1).cmd = 'viterbi';
  tomo_collate.surfdata_cmds(end).surf_names = 'bottom viterbi';
  tomo_collate.surfdata_cmds(end).visible = false;
  tomo_collate.surfdata_cmds(end).smooth_weight = 55; % schu
  tomo_collate.surfdata_cmds(end).smooth_var = -1; 
  tomo_collate.surfdata_cmds(end).repulsion = 150; % schu
  tomo_collate.surfdata_cmds(end).egt_weight = 10; 
  tomo_collate.surfdata_cmds(end).ice_bin_thr = 3;
  tomo_collate.surfdata_cmds(end).CF_sensory_distance = 50;
  tomo_collate.surfdata_cmds(end).CF_max_cost = 200;
  tomo_collate.surfdata_cmds(end).CF_lambda = 0.075;
  
  tomo_collate.surfdata_cmds(end+1).cmd = 'trws';
  tomo_collate.surfdata_cmds(end).surf_names = {'bottom trws','bottom'};
  tomo_collate.surfdata_cmds(end).visible = true;
  tomo_collate.surfdata_cmds(end).smooth_weight = [22 22];
  tomo_collate.surfdata_cmds(end).smooth_var = 32;
  tomo_collate.surfdata_cmds(end).max_loops = 50;
  tomo_collate.surfdata_cmds(end).CF_sensory_distance = 50;
  tomo_collate.surfdata_cmds(end).CF_max_cost = 200;
  tomo_collate.surfdata_cmds(end).CF_lambda = 0.075;
  
%   tomo_collate.surfdata_cmds(end+1).cmd = 'dem';
%   tomo_collate.surfdata_cmds(end).surf_names = 'bottom dem';
%   tomo_collate.surfdata_cmds(end).visible = false;
%   tomo_collate.surfdata_cmds(end).plot_name_values = {'color','red','marker','^'};
%   tomo_collate.surfdata_cmds(end).dem_fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_DEM_nofilt/20140401_03/20140401_03_007_bottom.tif';
%   tomo_collate.surfdata_cmds(end).dem_bad_value = -32767;
%   
%   tomo_collate.surfdata_cmds(end+1).cmd = 'dem';
%   tomo_collate.surfdata_cmds(end).surf_names = 'bottom dem2';
%   tomo_collate.surfdata_cmds(end).visible = false;
%   tomo_collate.surfdata_cmds(end).plot_name_values = {'color','black','marker','^'};
%   tomo_collate.surfdata_cmds(end).dem_fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_DEM_nofilt/20140401_03/20140401_03_033_bottom.tif';
%   tomo_collate.surfdata_cmds(end).dem_bad_value = -32767;
  
  % .fuse_images_flag: runs fuse_images.m when true
  tomo_collate.fuse_images_flag = true;
  
  % .add_icemask_surfacedem_flag: runs add_icemask_surfacedem.m when true
  tomo_collate.add_icemask_surfacedem_flag = true;
  
  % create_surfData_flag: runs create_surfData.m when true
  tomo_collate.create_surfData_flag = true;

  param_override.tomo_collate = tomo_collate;
  
elseif strcmpi(example_setup,'grid')
  %% Grid multiwaveform fuse example
  params = read_param_xls(ct_filename_param('rds_param_2018_Greenland_P3.xls'),'','post');
  
%   params = ct_set_params(params,'cmd.generic',0);
%   params = ct_set_params(params,'cmd.generic',1,'day_seg','20180406_01');
%   params = ct_set_params(params,'cmd.frms',[]);

  % 20180418_06_013 --> 3D object in ice
  % 20180405_01_056 --> 3D object in ice
  params = ct_set_params(params,'cmd.generic',0);
%   params = ct_set_params(params,'cmd.generic',1,'day_seg','20180404_02');
%   params = ct_set_params(params,'cmd.frms',,'day_seg','20180404_02');
%   params = ct_set_params(params,'cmd.generic',1,'day_seg','20180406_01');
%   params = ct_set_params(params,'cmd.frms',[2],'day_seg','20180406_01');
%   params = ct_set_params(params,'cmd.generic',1,'day_seg','20180418_06');
%   params = ct_set_params(params,'cmd.frms',[13],'day_seg','20180418_06');
%   params = ct_set_params(params,'cmd.generic',1,'day_seg','20180405_01');
%   params = ct_set_params(params,'cmd.frms',[56],'day_seg','20180405_01');
  params = ct_set_params(params,'cmd.generic',1,'day_seg','20180418_05');
  params = ct_set_params(params,'cmd.frms',[],'day_seg','20180418_05');

  params = ct_set_params(params,'array.method','music');
  params = ct_set_params(params,'array.three_dim.en',1);
  for param_idx = 1:length(params)
    if params(param_idx).cmd.generic
      if length(params(param_idx).radar.wfs) == 4
        params(param_idx).array.imgs = {[ones([7 1]),[6:12].'],[2*ones([7 1]),[6:12].'],[3*ones([7 1]),[6:12].'],[4*ones([7 1]),[6:12].']};
        tomo_collate.imgs = {[1 2],[3 4]};
        tomo_collate.img_comb = [0 3e-6 1e-6];
        params = ct_set_params(params,'array.Nsv',64);
        params = ct_set_params(params,'array.Nsig',2);
      elseif length(params(param_idx).radar.wfs) == 6
        params(param_idx).array.imgs = {[ones([7 1]),[6:12].'],[2*ones([7 1]),[6:12].'],[3*ones([7 1]),[6:12].'],[4*ones([7 1]),[6:12].'],[5*ones([7 1]),[6:12].'],[6*ones([7 1]),[6:12].']};
        tomo_collate.imgs = {[1 2],[3 4],[5 6]};
        tomo_collate.img_comb = [0 3e-6 1e-6 0 10e-6 3e-6];
        params = ct_set_params(params,'array.Nsv',64);
        params = ct_set_params(params,'array.Nsig',2);
      elseif length(params(param_idx).radar.wfs) == 2
        % imgs: No change required
        tomo_collate.imgs = {1};
        tomo_collate.img_comb = [];
        params = ct_set_params(params,'array.Nsv',128);
        params = ct_set_params(params,'array.Nsig',4);
      else
        keyboard
      end
    end
  end
  
  % .in_path: ct_filename_out directory to use at input, fused image will be stored here.
  tomo_collate.in_path = 'music_imgs4_Nsig2';
  
  % .surf_out_path: ct_filename_out directory to use at output for surfData
  tomo_collate.surf_out_path = 'surfData';
%   tomo_collate.surf_out_path = 'surfData_englacial';

  % .imgs: list of images II to use from .in_path (Data_img_II*.mat). These
  %   should be listed from left most beam to right most beam when
  %   horizontal fuse is used. They should be listed from top to bottom
  %   when vertical fuse is used.
  tomo_collate.master_img_idx = 1;
  
  % .img_comb: Same as get_heights and combine worksheets. This is
  %   used for vertical using only. For N images,
  %   there will be (N-1) combines performed. Each combine is described by
  %   three numbers so that there should be (N-1)*3 elements in this
  %   vector. Each set of three numbers indicates the following:
  %     1st element: Minimum time to begin combine
  %     2nd element: Minimum time after ice surface to begin combine
  %     3rd element: Time at end of the preceeding waveform to not use
  
  % .fuse_columns: aligns with .imgs, each entry should contain 2xN-1
  % entries where N is the length of the corresponding cell in .imgs. Each
  % column of 2 numbers represents the start/stop columns to blend with.
  % There should be one column per interface between two images that needs
  % to be blended in the horizontal dimension. E.g. If blending 3 images in
  % the horizontal dimension, the entry should be 2x2. If blending 2
  % images, the entry should be 2x1. If there is only one image for a
  % particular vertical index, then fuse_columns should be empty.
  tomo_collate.fuse_columns = {[30 36].',[30 36].'};
  
  % .sv_cal_fn: filename containing steering vector calibration, leave empty to not use
  tomo_collate.sv_cal_fn = '';
  
  % .ice_mask_fn: filename of ice mask, leave empty to not use
  tomo_collate.ice_mask_fn = 'greenland/IceMask/GimpIceMask_90m_v1.1.tif';

  % .dem_guard: additional region around flight line to search for DEM points
  %   Setting too high slows the process down, setting too low will miss
  %   DEM points needed to properly represent the surface.
  tomo_collate.dem_guard = 16e3;
  
  % .dem_per_slice_guard: additional region around each slice to search for DEM points
  %   Setting too high slows the process down, setting too low will miss
  %   DEM points needed to properly represent the surface.
  tomo_collate.dem_per_slice_guard = 2500;
  
  % .bounds_relative: DOA bins and along-track slices to trim off from each edge [top bottom left right]
  tomo_collate.bounds_relative = [3 2 0 0];
  
  % .layer_params: parameter structure for opsLoadLayer (first layer should
  %   be ice top and second layer should be ice bottom)
  tomo_collate.layer_params = struct('name','surface','source','layerdata');
  tomo_collate.layer_params(2).name = 'bottom';
  tomo_collate.layer_params(2).source = 'layerdata';
%   tomo_collate.layer_params(2).name = 'float';
%   tomo_collate.layer_params(2).source = 'ops';
  
  % surfData_mode: surfData mode ('overwrite','fillgaps', or 'append', note that append with the
  %   same surface name as an existing surface will overwrite that surface whereas fillgaps
  %   will leave the surface untouched if it already exists)
  tomo_collate.surfData_mode = 'overwrite';
  
  % surfdata_cmds: surfdata commands to run
  tomo_collate.surfdata_cmds = [];
  tomo_collate.surfdata_cmds(end+1).cmd = 'trws';
  tomo_collate.surfdata_cmds(end).surf_names = {'bottom trws','bottom'};
  tomo_collate.surfdata_cmds(end).visible = true;
  
  % .fuse_images_flag: runs fuse_images.m when true
  tomo_collate.fuse_images_flag = true;
  
  % .add_icemask_surfacedem_flag: runs add_icemask_surfacedem.m when true
  tomo_collate.add_icemask_surfacedem_flag = true;
  
  % create_surfData_flag: runs create_surfData.m when true
  tomo_collate.create_surfData_flag = true;
  
  param_override.tomo_collate = tomo_collate;

else
  error('%s is not a valid example_setup.', example_setup);
end

dbstop if error;
param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
% param_override.cluster.rerun_only = true;
% param_override.cluster.desired_time_per_job  = 240*60;
% param_override.cluster.cpu_time_mult  = 2;
% param_override.cluster.mem_mult  = 2;

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
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  ctrl_chain{end+1} = tomo.collate(param,param_override);
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);
