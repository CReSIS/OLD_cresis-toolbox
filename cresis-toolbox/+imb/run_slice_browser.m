% script imb.run_slice_browser
%
% Author: Sravya Athinarapu, Elijah Paden, John Paden, Jordan Sprick

%% User Settings
% =========================================================================

param.radar_name = 'mcords3';
param.season_name = '2014_Greenland_P3';
out_type = 'paden';
surfdata_source = 'surfData2';
% out_type = 'CSA_music';
% surfdata_source = 'surfData';
param.day_seg = '20140401_03';
frm = 37;
% frm = 39;
% frm = 43;
% frm = 44;
% frm = 45;
% frm = 46;
% frm = 47;
% frm = 48;

%% Automated Section
% =========================================================================

fn = fullfile(ct_filename_out(param,out_type,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
if ~exist('run_slice_browser_fn','var') || ~strcmp(run_slice_browser_fn,fn)
  fprintf('Loading data (%s)\n', datestr(now));
  mdata = load(fn);
  
  geotiff_fn = ct_filename_gis(param,fullfile('canada','Landsat-7','Canada_90m.tif'));
  ice_mask_fn = ct_filename_gis(param,fullfile('canada','ice_mask','03_rgi50_ArcticCanadaNorth','03_rgi50_ArcticCanadaNorth.mat'));
  
  proj = geotiffinfo(geotiff_fn);
  ice_mask = load(ice_mask_meta_fn,'R','X','Y','proj');
  fid = fopen(ice_mask_fn,'r');
  ice_mask.mask = logical(fread(fid,[length(ice_mask.Y),length(ice_mask.X)],'uint8'));
  fclose(fid);
  [DEM, R, tmp] = geotiffread(geotiff_fn);
  
  run_slice_browser_fn = fn;
  fprintf('  Done loading data (%s)\n', datestr(now));
end

sb_param = [];
sb_param.layer_fn = fullfile(ct_filename_out(param,surfdata_source,'CSARP_surfData'),sprintf('Data_%s_%03d.mat',param.day_seg,frm));

%% Call slice_browser
try; delete(obj); end;
obj = imb.slice_browser(10*log10(mdata.Topography.img),[],sb_param);

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

try; delete(detect_tool); end;
detect_tool = imb.slicetool_detect();
custom_data.mu = mdata.Topography.mu;
custom_data.sigma = mdata.Topography.sigma;
detect_tool.set_custom_data(custom_data);
obj.insert_tool(detect_tool);

try; delete(threshold_tool); end;
threshold_tool = imb.slicetool_threshold();
custom_data.ice_mask = mdata.ice_mask;
custom_data.theta = mdata.param_combine.array_param.theta;
custom_data.img = mdata.Topography.img;
custom_data.Time = mdata.Time;
custom_data.sb = obj;
threshold_tool.set_custom_data(custom_data);
obj.insert_tool(threshold_tool);

try; delete(extract_tool); end;
extract_tool = imb.slicetool_extract();
custom_data.mu = mdata.Topography.mu;
custom_data.sigma = mdata.Topography.sigma;
extract_tool.set_custom_data(custom_data);
obj.insert_tool(extract_tool);

try; delete(max_tool); end;
max_tool = imb.slicetool_max();
obj.insert_tool(max_tool);
