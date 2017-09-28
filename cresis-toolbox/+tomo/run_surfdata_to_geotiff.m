%% User Settings

if 0
  params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'20140506_01','post');
  params.cmd.generic = 1;
  params.cmd.frms = [3];
  
  % surfdata_source: input surfData directory (ct_filename_out)
  surfdata_source = 'paden_surfData';
  
  % input_dir_name: input radar 3D image directory (ct_filename_out)
  input_dir_name = 'CSA_music';
  
  % output_dir_name: string containing output directory (ct_filename_out)
  output_dir_name = 'paden_DEM';
  
  % geotiff_fn: the projection information is taken from this file and this
  %   file is used for creating the maps
  geotiff_fn = ct_filename_gis(param,fullfile('canada','Landsat-7','Canada_90m.tif'));

  % DOA_trim: remove this many direction of arrival bins on each side of the
  %   image
  DOA_trim = 3;
  
  theta_cal_fn = ct_filename_ct_tmp(rmfield(params(1),'day_seg'),'','sv_calibration','theta_cal.mat');
  
else
  params = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO.xls'),'20091224_01','post');
  params.cmd.generic = 1;
  params.cmd.frms = [16:18];
  
  % surfdata_source: input surfData directory (ct_filename_out)
  surfdata_source = 'paden_surfData';
  
  % input_dir_name: input radar 3D image directory (ct_filename_out)
  input_dir_name = 'nick_music';
  
  % output_dir_name: string containing output directory (ct_filename_out)
  output_dir_name = 'paden_DEM';
  
  % geotiff_fn: the projection information is taken from this file and this
  %   file is used for creating the maps
  geotiff_fn = ct_filename_gis(param,fullfile('antarctica','Landsat-7','Antarctica_LIMA_480m.tif'));

  % DOA_trim: remove this many direction of arrival bins on each side of the
  %   image
  DOA_trim = 15;
  
  theta_cal_fn = '';
  
end

% active_surfs: cell array of surfData surface names
active_surfs = {'ice surface','bottom'};

% ice_mask_ref: string containing the reference citation for the ice mask
ice_mask_ref = 'http://www.glims.org/RGI/rgi50_dl.html';
% ice_mask_ref = [];

% geotiff_ref: string containing the reference citation for the geotiff
% geotiff_ref = [];
geotiff_ref = 'https://landsat.usgs.gov/';

% DEM_ref: string containing the reference citation for the surface DEM
DEM_ref = 'https://theia.cnes.fr/rocket/#/search?page=3&collection=Spirit&view=default';
% DEM_ref = ct_filename_gis(gRadar,fullfile('antarctia','Landsat-7','Antarctica_LIMA.tif'));
% DEM_ref = [];

%% Automated Section

for param_idx = 1:length(params)
  
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  param.dem.geotiff_fn = geotiff_fn;
  param.dem.active_surfs = active_surfs;
  param.dem.DOA_trim = DOA_trim;
  param.dem.ice_mask_ref = ice_mask_ref;
  param.dem.geotiff_ref = geotiff_ref;
  param.dem.DEM_ref = DEM_ref;
  param.dem.surfdata_source = surfdata_source;
  param.dem.input_dir_name = input_dir_name;
  param.dem.output_dir_name = output_dir_name;
  param.dem.theta_cal_fn = theta_cal_fn;
  
  tomo.surfdata_to_geotiff(param);

end