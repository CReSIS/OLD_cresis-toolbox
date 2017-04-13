%% User Settings

% params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'20140401_03','post');
% params = read_param_xls(ct_filename_param('rds_param_2014_Antarctica_DC8.xls'),'20141121_05','post');
params = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'),'20161014_05','post');
params.cmd.generic = 1;
% params.cmd.frms = [37,39,43:48];
% params.cmd.frms = 46;
% params.cmd.frms = 47:48;
% params.cmd.frms = 4;
params.cmd.frms = 21;

out_type = 'music3D_jordan';

% active_surfs = {'ice surface','bottom'};
active_surfs = {'bottom'};
DOA_trim = 3;

% geotiff_fn = ct_filename_gis(gRadar,fullfile('canada','Landsat-7','Canada_90m.tif'));
% geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_surface.tif';
geotiff_fn = ct_filename_gis(gRadar,fullfile('antarctica','Landsat-7','Antarctica_LIMA.tif'));

% ice_mask_ref = 'http://www.glims.org/RGI/rgi50_dl.html';
ice_mask_ref = [];
geotiff_ref = [];
% DEM_ref = 'https://theia.cnes.fr/rocket/#/search?page=3&collection=Spirit&view=default';
DEM_ref = ct_filename_gis(gRadar,fullfile('antarctia','Landsat-7','Antarctica_LIMA.tif'));
% DEM_ref = [];

theta_cal_fn = ct_filename_ct_tmp(rmfield(params(1),'day_seg'),'','sv_calibration','theta_cal.mat');

%% Automated Section

for param_idx = 1:length(params)
  
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  param.geotiff_fn = geotiff_fn;
  param.active_surfs = active_surfs;
  param.DOA_trim = DOA_trim;
  param.ice_mask_ref = ice_mask_ref;
  param.geotiff_ref = geotiff_ref;
  param.DEM_ref = DEM_ref;
  param.out_type = out_type;
  param.theta_cal_fn = theta_cal_fn;
  
  tomo.surfdata_to_geotiff(param);

end