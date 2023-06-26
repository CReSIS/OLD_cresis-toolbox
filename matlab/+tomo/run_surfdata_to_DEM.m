% script run_surfdata_to_DEM
%
% Run script for tomo.surfdata_to_DEM
%
% Authors: Jordan Sprick, Nick Holschuh, John Paden, Victor Berger
%
% See also: tomo.surfdata_to_DEM

%% User Settings
params = [];
dem = [];
if 1
  dem.doa_method_flag = true;
   params = read_param_xls(ct_filename_param('rds_param_2011_Greenland_P3.xls'));
%   params.cmd.generic = 1;
%   params.cmd.frms = [];
  params = ct_set_params(params,'cmd.generic',0);

  params = ct_set_params(params,'cmd.generic',1,'day_seg','20110317_03');
  params = ct_set_params(params,'cmd.frms',[1]);
  
  % surfdata_source: input surfData directory (ct_filename_out)
  dem.surfdata_source = 'test_surfData_doa';
%   dem.surfdata_source = 'surfData_englacial';
  
  % input_dir_name: input radar 3D image directory (ct_filename_out)
  dem.input_dir_name = 'test_music3D_mle';
  
  % output_dir_name: string containing output directory (ct_filename_out)
  dem.output_dir_name = 'test_DEM_doa';
%   dem.output_dir_name = 'DEM_englacial';
  
  % geotiff_fn: the projection information is taken from this file and this
  %   file is used for creating the maps
  dem.geotiff_fn = ct_filename_gis(params(1),fullfile('greenland','Landsat-7','Greenland_natural_90m.tif'));

  % DOA_trim: remove this many direction of arrival bins on each side of the
  %   image
%   dem.DOA_trim = 5;

  % DOA limits for slice browsing (outside this limits will not be displayed)
  dem.doa_limits = [-25 25]; 
  
  % med_filt: medfilt2 arguments for spatial filtering (leave blank for no filtering)
  dem.med_filt = [];
  
  % figure_dots_per_km: scalar representing the number of pixels per km in
  %   the figures
  dem.figure_dots_per_km = 20;
  
  % theta_cal_fn: string containing Theta calibration filename (leave blank
  %   if no calibration used)
%   dem.theta_cal_fn = ct_filename_ct_tmp(rmfield(params(1),'day_seg'),'','sv_calibration','theta_cal.mat');
  
  % ice_mask_ref: string containing the reference citation for the ice mask
  dem.ice_mask_ref = 'http://www.glims.org/RGI/rgi50_dl.html';
  
  % geotiff_ref: string containing the reference citation for the geotiff
  dem.geotiff_ref = 'https://landsat.usgs.gov/';
  
  % DEM_ref: string containing the reference citation for the surface DEM
  dem.DEM_ref = 'DEMs provided by the Polar Geospatial Center under NSF OPP awards 1043681, 1559691 and 1542736.';
  
  % ice_top: string containing surface name for the top of the ice (used for
  % refraction)
  dem.ice_top = 'top';
  
  % surface_names: cell array of strings containing surface names that output
  % data products will be generated for
  dem.surface_names = {'top','bottom'};
%   dem.surface_names = {'top'};
  
  % quality_surface_names: cell array of strings containing surface names
  % for the quality surface that corresponds to each entry in
  % surface_names.
  dem.quality_surface_names = {'top quality','bottom quality'};
%   dem.quality_surface_names = {'bottom quality'};
  
  dem.grid_spacing = 25;
  dem.bad_geotiff_value = 32767;
  
  param_override = struct('dem',dem);
elseif 0
  params = read_param_xls(ct_filename_param('rds_param_2018_Greenland_P3.xls'),'','post');
%   params.cmd.generic = 1;
%   params.cmd.frms = [];
  params = ct_set_params(params,'cmd.generic',0);
%   params = ct_set_params(params,'cmd.generic',1,'day_seg','20180404_02');
%   params = ct_set_params(params,'cmd.frms',[1]);
%   params = ct_set_params(params,'cmd.generic',1,'day_seg','20180418_06');
%   params = ct_set_params(params,'cmd.frms',[13]);
%   params = ct_set_params(params,'cmd.generic',1,'day_seg','20180405_01');
%   params = ct_set_params(params,'cmd.frms',[56]);
  params = ct_set_params(params,'cmd.generic',1,'day_seg','20180406_01');
  params = ct_set_params(params,'cmd.frms',[2]);
  
  % surfdata_source: input surfData directory (ct_filename_out)
  dem.surfdata_source = 'surfData';
%   dem.surfdata_source = 'surfData_englacial';
  
  % input_dir_name: input radar 3D image directory (ct_filename_out)
  dem.input_dir_name = 'music_imgs4_Nsig2';
  
  % output_dir_name: string containing output directory (ct_filename_out)
  dem.output_dir_name = 'DEM';
%   dem.output_dir_name = 'DEM_englacial';
  
  % geotiff_fn: the projection information is taken from this file and this
  %   file is used for creating the maps
  dem.geotiff_fn = ct_filename_gis(params(1),fullfile('greenland','Landsat-7','Greenland_natural_90m.tif'));

  % DOA_trim: remove this many direction of arrival bins on each side of the
  %   image
  dem.DOA_trim = 0%;5; % Set this to 0 for DOA methods
  
  % med_filt: medfilt2 arguments for spatial filtering (leave blank for no filtering)
  dem.med_filt = [];%[5 17];
  
  % figure_dots_per_km: scalar representing the number of pixels per km in
  %   the figures
  dem.figure_dots_per_km = 20;
  
  % theta_cal_fn: string containing Theta calibration filename (leave blank
  %   if no calibration used)
%   dem.theta_cal_fn = ct_filename_ct_tmp(rmfield(params(1),'day_seg'),'','sv_calibration','theta_cal.mat');
  
  % ice_mask_ref: string containing the reference citation for the ice mask
  dem.ice_mask_ref = 'http://www.glims.org/RGI/rgi50_dl.html';
  
  % geotiff_ref: string containing the reference citation for the geotiff
  dem.geotiff_ref = 'https://landsat.usgs.gov/';
  
  % DEM_ref: string containing the reference citation for the surface DEM
  dem.DEM_ref = 'DEMs provided by the Polar Geospatial Center under NSF OPP awards 1043681, 1559691 and 1542736.';
  
  % ice_top: string containing surface name for the top of the ice (used for
  % refraction)
  dem.ice_top = 'top';
  
  % surface_names: cell array of strings containing surface names that output
  % data products will be generated for
  dem.surface_names = {'top','bottom'};
%   dem.surface_names = {'bottom'};
  
  % quality_surface_names: cell array of strings containing surface names
  % for the quality surface that corresponds to each entry in
  % surface_names.
  dem.quality_surface_names = {'top quality','bottom quality'};
%   dem.quality_surface_names = {'bottom quality'};
  
  dem.grid_spacing = 25;
  dem.bad_geotiff_value = 32767;
  
  param_override = struct('dem',dem);
  
elseif 0
  params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'','post');
%   params.cmd.generic = 1;
%   params.cmd.frms = [];
  
  % surfdata_source: input surfData directory (ct_filename_out)
  dem.surfdata_source = 'CSARP_post/surfData';
  
  % input_dir_name: input radar 3D image directory (ct_filename_out)
  dem.input_dir_name = 'CSARP_post/music3D';
  
  % output_dir_name: string containing output directory (ct_filename_out)
  dem.output_dir_name = 'DEM';
  
  % geotiff_fn: the projection information is taken from this file and this
  %   file is used for creating the maps
  dem.geotiff_fn = ct_filename_gis(params(1),fullfile('canada','Landsat-7','Canada_90m.tif'));

  % DOA_trim: remove this many direction of arrival bins on each side of the
  %   image
  dem.DOA_trim = 5;
  
  % med_filt: medfilt2 arguments for spatial filtering (leave blank for no filtering)
  dem.med_filt = [5 17];
  
  % figure_dots_per_km: scalar representing the number of pixels per km in
  %   the figures
  dem.figure_dots_per_km = 20;
  
  % theta_cal_fn: string containing Theta calibration filename (leave blank
  %   if no calibration used)
  dem.theta_cal_fn = ct_filename_ct_tmp(rmfield(params(1),'day_seg'),'','sv_calibration','theta_cal.mat');
  
  % ice_mask_ref: string containing the reference citation for the ice mask
  dem.ice_mask_ref = 'http://www.glims.org/RGI/rgi50_dl.html';
  
  % geotiff_ref: string containing the reference citation for the geotiff
  dem.geotiff_ref = 'https://landsat.usgs.gov/';
  
  % DEM_ref: string containing the reference citation for the surface DEM
  dem.DEM_ref = 'DEMs provided by the Polar Geospatial Center under NSF OPP awards 1043681, 1559691 and 1542736.';
  
  % ice_top: string containing surface name for the top of the ice (used for
  % refraction)
  dem.ice_top = 'top';
  
  % surface_names: cell array of strings containing surface names that output
  % data products will be generated for
%   dem.surface_names = {'top','bottom'};
  dem.surface_names = {'bottom'};
  
  % quality_surface_names: cell array of strings containing surface names
  % for the quality surface that corresponds to each entry in
  % surface_names.
  dem.quality_surface_names = {'top quality','bottom quality'};
%   dem.quality_surface_names = {'bottom quality'};
  
  dem.grid_spacing = 25;
  dem.bad_geotiff_value = 32767;
  
  param_override = struct('dem',dem);
  
else
  params = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO.xls'),'20091224_01','post');
  params.cmd.generic = 1;
  params.cmd.frms = [];
  
  % surfdata_source: input surfData directory (ct_filename_out)
  dem.surfdata_source = 'surfData';
  
  % input_dir_name: input radar 3D image directory (ct_filename_out)
  dem.input_dir_name = 'music3D';
  
  % output_dir_name: string containing output directory (ct_filename_out)
  dem.output_dir_name = 'DEM';
  
  % geotiff_fn: the projection information is taken from this file and this
  %   file is used for creating the maps
  dem.geotiff_fn = ct_filename_gis(params(1),fullfile('antarctica','Landsat-7','Antarctica_LIMA_480m.tif'));
  
  % DOA_trim: remove this many direction of arrival bins on each side of the
  %   image
  dem.DOA_trim = 15;
  
  % med_filt: medfilt2 arguments for spatial filtering (leave blank for no filtering)
  dem.med_filt = [13 9];
  
  % figure_dots_per_km: scalar representing the number of pixels per km in
  %   the figures
  dem.figure_dots_per_km = 20;
  
  % theta_cal_fn: string containing Theta calibration filename (leave blank
  %   if no calibration used)
  dem.theta_cal_fn = '';
  
  % ice_mask_ref: string containing the reference citation for the ice mask
  dem.ice_mask_ref = '';
  
  % geotiff_ref: string containing the reference citation for the geotiff
  dem.geotiff_ref = 'https://landsat.usgs.gov/';
  
  % DEM_ref: string containing the reference citation for the surface DEM
  dem.DEM_ref = '';
  
  % ice_top: string containing surface name for the top of the ice (used for
  % refraction)
  dem.ice_top = 'top';
  
  % surface_names: cell array of strings containing surface names that output
  % data products will be generated for
  dem.surface_names = {'top','bottom'};
  %dem.surface_names = {'bottom'};
  
  % quality_surface_names: cell array of strings containing surface names
  % for the quality surface that corresponds to each entry in
  % surface_names.
  dem.quality_surface_names = {'top quality','bottom quality'};
  %dem.quality_surface_names = {'bottom quality'};
  
  dem.grid_spacing = 25;
  dem.bad_geotiff_value = 32767;
  
  param_override = struct('dem',dem);
  
end

%% Automated data posting section
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
  
  tomo.surfdata_to_DEM(param,param_override);
end
