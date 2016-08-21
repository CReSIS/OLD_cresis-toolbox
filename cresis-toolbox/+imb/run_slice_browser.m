% script imb.run_slice_browser
%
% Author: Elijah Paden, John Paden

% mdata = [];
if ~exist('run_slice_browser_init','var') || ~run_slice_browser_init
    
    mdata = load('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_CSA_music/20140401_03/Data_20140401_03_039.mat');
%     mdata = load('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_CSA_music/20140401_03/Data_20140401_03_044_old.mat');
    mdata.ice_mask = logical(mdata.ice_mask);
    
    %% Save layer data into a file
    layer = [];
    layer(1).x = repmat((1:64).',[1 size(mdata.twtt,2)]);
    layer(1).y = interp1(mdata.Time,1:length(mdata.Time),mdata.twtt);
    layer(1).plot_name_values = {'color','black','marker','x'};
    layer(1).name = 'surface';
    layer(1).surf_layer = 1;
    layer(1).active_layer = 1;
    layer(1).control_layer = 1;
    layer(1).mask_layer = 4;
    
    layer(2).x =  repmat((1:64).',[1 size(mdata.twtt,2)]);
    layer(2).y = mdata.bottom_surface.';
    layer(2).plot_name_values = {'color','blue','marker','^'};
    layer(2).name = 'bottom';
    layer(2).surf_layer = 1;
    layer(2).active_layer = 2;
    layer(2).control_layer = 3;
    layer(2).mask_layer = 4;
    
    layer(3).x = repmat((1:64).',[1 size(mdata.twtt,2)]);
    layer(3).y = NaN * zeros(size(layer(1).y));
    layer(3).y(33,:) = interp1(mdata.Time,1:length(mdata.Time),mdata.Bottom);
    layer(3).plot_name_values = {'color','magenta','marker','+'};
    layer(3).name = 'Bcontrol';
    layer(3).surf_layer = 1;
    layer(3).active_layer = 2;
    layer(3).control_layer = 3;
    layer(3).mask_layer = 4;
    
    layer(4).x = repmat((1:64).',[1 size(mdata.twtt,2)]);
    layer(4).y = mdata.ice_mask;
    %   layer(4).y = layer(1).y;
    %   layer(4).y(~logical(mdata.ice_mask)) = NaN;
    layer(4).plot_name_values = {'color','white','marker','x'};
    layer(4).name = 'ice mask';
    layer(4).surf_layer = 1;
    layer(4).active_layer = 2;
    layer(4).control_layer = 3;
    layer(4).mask_layer = 4;
    
    param = [];
    param.layer_fn = '~/layer_data.mat';
    %if ~exist(param.layer_fn,'file')
    save(param.layer_fn,'layer')
    %end
    run_slice_browser_init = true;
    
    %JORDAN
    theta_cal = load('/cresis/snfs1/dataproducts/ct_data/ct_tmp/sv_calibration/rds/2014_Greenland_P3/theta_cal.mat');
    mdata.theta = theta_cal.theta;
    % geotiff_fn = 'X:\GIS_data\canada\Landsat-7/'
    geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/canada/Landsat-7/Canada_90m.tif';
    % ice_mask_fn = X:\GIS_data\canada\ice_mask\03_rgi50_ArcticCanadaNorth
    ice_mask_fn = '/cresis/snfs1/dataproducts/GIS_data/canada/ice_mask/03_rgi50_ArcticCanadaNorth/03_rgi50_ArcticCanadaNorth.mat';
    fprintf('Loading geotiff and ice mask\n');
    proj = geotiffinfo(geotiff_fn);
    ice_mask = load(ice_mask_fn);
    [DEM, R, tmp] = geotiffread(geotiff_fn);
    fprintf('Done Loading\n');
    
end

%% Call slice_browser
try; delete(obj); end;
h_control_fig = figure();
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
custom_data.ice_mask = mdata.ice_mask;
custom_data.bottom = interp1(mdata.Time,1:length(mdata.Time),mdata.Bottom);
extract_tool.set_custom_data(custom_data);
obj.insert_tool(extract_tool);
