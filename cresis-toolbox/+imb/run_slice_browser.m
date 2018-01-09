% script imb.run_slice_browser
%
% Script for running imb.slice_browser.m
%
% Run "clear run_slice_browser_fn;" whenever you change frames or geotiff
%
% Author: Sravya Athinarapu, Elijah Paden, John Paden, Jordan Sprick

%% User Settings
% =========================================================================

if 1
  param.radar_name = 'rds';
  param.season_name = '2014_Greenland_P3';
  out_type = 'music3D';
  surfdata_source = 'surfData';
  param.day_seg = '20140325_05';
  frm = 2;
  geotiff_fn = ct_filename_gis(param,fullfile('canada','Landsat-7','Canada_90m.tif'));
  ice_mask_fn = ct_filename_gis(param,fullfile('canada','ice_mask','03_rgi50_ArcticCanadaNorth','03_rgi50_ArcticCanadaNorth.bin'));
  bounds_relative = [3 2 0 0];
  
elseif 0
  param.radar_name = 'rds';
  param.season_name = '2009_Antarctica_TO';
  out_type = 'music3D';
  surfdata_source = 'surfData';
  param.day_seg = '20091224_01';
  frm = 26;
  geotiff_fn = ct_filename_gis(param,fullfile('antarctica','Landsat-7','Antarctica_LIMA_480m.tif'));
  ice_mask_fn = '';
  bounds_relative = [8 8 0 0];
  
else
  param.radar_name = 'rds';
  param.season_name = '2016_Antarctica_DC8';
  out_type = 'music3D';
  surfdata_source = 'surfData';
  param.day_seg = '20161117_06';
  frm = 1;
  geotiff_fn = ct_filename_gis(param,fullfile('antarctica','Landsat-7','Antarctica_LIMA_480m.tif'));
  ice_mask_fn = '';
  bounds_relative = [8 8 0 0];
end

%% Automated Section
% =========================================================================

fn = fullfile(ct_filename_out(param,out_type,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
if ~exist('run_slice_browser_fn','var') || ~strcmp(run_slice_browser_fn,fn)
  fprintf('Loading data (%s)\n', datestr(now));
  mdata = load(fn);
  
  proj = geotiffinfo(geotiff_fn);
  [DEM, R, tmp] = geotiffread(geotiff_fn);
  
  run_slice_browser_fn = fn;
  fprintf('  Done loading data (%s)\n', datestr(now));
end

sb_param = [];
sb_param.surfdata_fn = fullfile(ct_filename_out(param,surfdata_source,'CSARP_surfData'),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
if ~exist(sb_param.surfdata_fn)
  sb_param.surfdata_fn = '';
end
sb_param.bounds_relative = bounds_relative;

%% Call slice_browser
try; delete(obj); end;
obj = imb.slice_browser(10*log10(mdata.Topography.img),[],sb_param);

try; delete(detect_tool); end;
detect_tool = imb.slicetool_detect();
if isfield(mdata.Topography,'mu')
  % This field is only present after training has been run
  custom_data.mu = mdata.Topography.mu;
  custom_data.sigma = mdata.Topography.sigma;
end
detect_tool.set_custom_data(custom_data);
obj.insert_tool(detect_tool);

try; delete(extract_tool); end;
extract_tool = imb.slicetool_extract();
extract_tool.set_custom_data(custom_data);
obj.insert_tool(extract_tool);

try; delete(max_tool); end;
max_tool = imb.slicetool_max();
obj.insert_tool(max_tool);

try; delete(quality_tool); end;
quality_tool = imb.slicetool_quality();
obj.insert_tool(quality_tool);

try; delete(delete_tool); end;
delete_tool = imb.slicetool_delete();
obj.insert_tool(delete_tool);

try; delete(threshold_tool); end;
threshold_tool = imb.slicetool_threshold();
obj.insert_tool(threshold_tool);

if ~isempty(ice_mask_fn)
  [ice_mask_fn_dir ice_mask_fn_name] = fileparts(ice_mask_fn);
  ice_mask_mat_fn = fullfile(ice_mask_fn_dir,[ice_mask_fn_name '.mat']);
  ice_mask = load(ice_mask_mat_fn,'R','X','Y','proj');
  
  [fid,msg] = fopen(ice_mask_fn,'r');
  if fid < 1
    fprintf('Could not open file %s\n', ice_mask_bin_fn);
    error(msg);
  end
  ice_mask.mask = logical(fread(fid,[length(ice_mask.Y),length(ice_mask.X)],'uint8'));
  fclose(fid);
  
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
  custom_data.ice_mask_layer = 3;
  icemask_tool.set_custom_data(custom_data);
  obj.insert_tool(icemask_tool);
end
