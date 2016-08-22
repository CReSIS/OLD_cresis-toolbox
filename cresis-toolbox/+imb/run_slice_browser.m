% script imb.run_slice_browser
%
% Author: Elijah Paden, John Paden
close all;
% mdata = [];
if ~exist('run_slice_browser_init','var') || ~run_slice_browser_init
  
%   mdata = load('X:/ct_data/rds/2014_Greenland_P3/CSARP_CSA_music/20140401_03/Data_20140401_03_039.mat');
  mdata = load('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_CSA_music/20140401_03/Data_20140401_03_039.mat');
  mdata.ice_mask = logical(mdata.ice_mask);
  
  %% Save surf data into a file
  
  surf_data = load('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_surfData/20140401_03/Data_20140401_03_039');
%   surf_data = load('~/surf_data.mat');
  layer = surf_data.surf;
  
  param = [];
  param.layer_fn = '~/surf_data.mat';
  if ~exist(param.layer_fn,'file')
    save(param.layer_fn,'layer')
  end
  
  theta_cal = load('/cresis/snfs1/dataproducts/ct_data/ct_tmp/sv_calibration/rds/2014_Greenland_P3/theta_cal.mat');
%   theta_cal = load('X:\ct_data\ct_tmp\sv_calibration\rds\2014_Greenland_P3\theta_cal.mat');
  mdata.theta = theta_cal.theta;
  
  mdata.ice_mask = layer(find(strncmp({layer.name},'ice mask',8))).y;
  
%   geotiff_fn = 'X:\GIS_data\canada\Landsat-7\Canada_90m.tif';
  geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/canada/Landsat-7/Canada_90m.tif';
%   ice_mask_fn = 'X:\GIS_data\canada\ice_mask\03_rgi50_ArcticCanadaNorth\03_rgi50_ArcticCanadaNorth.mat';
%   ice_mask_fn = '/cresis/snfs1/dataproducts/GIS_data/canada/ice_mask/03_rgi50_ArcticCanadaNorth/03_rgi50_ArcticCanadaNorth.mat';
  ice_mask_fn = '~/03_rgi50_ArcticCanadaNorth_2.mat';
  fprintf('Loading geotiff and ice mask\n');
  proj = geotiffinfo(geotiff_fn);
  ice_mask = load(ice_mask_fn);
  [DEM, R, tmp] = geotiffread(geotiff_fn);
  fprintf('Done Loading\n');
  
  run_slice_browser_init = true;
  
end

%% Call slice_browser
try; delete(obj); end;
h_control_fig = figure(1); clf;
h_control_axes = axes('Parent',h_control_fig);
h_control_image = imagesc(lp(squeeze(mdata.Topography.img(:,33,:))),'Parent',h_control_axes);
colormap(parula(256))
obj = imb.slice_browser(lp(mdata.Topography.img),h_control_image,param);

try; delete(icemask_tool); end;
icemask_tool = imb.slicetool_icemask();
custom_data.DEM = DEM;
custom_data.R = R;
custom_data.ice_mask = ice_mask;
custom_data.proj = proj;
custom_data.ice_mask_fn = ice_mask_fn;
custom_data.mdata = mdata;
custom_data.sb = obj;
custom_data.reduce_flag = 1;
icemask_tool.set_custom_data(custom_data);
obj.insert_tool(icemask_tool);

try; delete(detect_tool); end;
detect_tool = imb.slicetool_detect();
custom_data.mu = mdata.Topography.mu;
custom_data.sigma = mdata.Topography.sigma;
custom_data.ice_mask = mdata.ice_mask;
custom_data.bottom = interp1(mdata.Time,1:length(mdata.Time),mdata.Bottom);
detect_tool.set_custom_data(custom_data);
obj.insert_tool(detect_tool);

try; delete(threshold_tool); end;
threshold_tool = imb.slicetool_threshold();
custom_data.ice_mask = mdata.ice_mask;
custom_data.mdata = mdata;
custom_data.sb = obj;
threshold_tool.set_custom_data(custom_data);
obj.insert_tool(threshold_tool);

try; delete(extract_tool); end;
extract_tool = imb.slicetool_extract();
custom_data.mu = mdata.Topography.mu;
custom_data.sigma = mdata.Topography.sigma;
custom_data.ice_mask = mdata.ice_mask;
custom_data.bottom = interp1(mdata.Time,1:length(mdata.Time),mdata.Bottom);
extract_tool.set_custom_data(custom_data);
obj.insert_tool(extract_tool);
