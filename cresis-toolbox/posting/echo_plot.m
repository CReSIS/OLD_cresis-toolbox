function [mdata] = echo_plot(echo_fn,layerdata_source,layer_names)
% [mdata] = echo_plot(echo_fn,layerdata_source,layer_names)
%
% Example of loading data with load_L1B.m and using elevation_compensation.m
%
% Plots data in three different ways:
%   Figure 1: time-delay on y-axis
%   Figure 2: range on y-axis
%   Figure 3: time-delay on y-axis
%
% Step 1: support list of layers using current layerdata.load_layers
% mdata = echo_plot(echo_fn,'layer') % In this case assume third argument is {'surface','bottom'}
% mdata = echo_plot(echo_fn,'layer',{'surface','bottom'}) % Assume first layer is the surface
% Step 2: support layer_params struct (see opsLoadLayers)
% layer_params = struct('name',{'surface','bottom'},'source','layerdata','layerdata_source','layer')
% mdata = echo_plot(echo_fn,layer_params)
% Step 3: allow wild character searches in the layer names
% [mdata,h] = echo_plot(echo_fn,'',{'surface','.*'})
% Step 4: return convenient list of GUI handles
% [mdata,h] = echo_plot(echo_fn,'layer',{'surface','bottom'})
%
% Examples:
%
% fn = 'IRMCR1B_V01_20130408_01_020.nc';
% mdata = echo_plot(fn);
%
% fn = '/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_post/CSARP_qlook/20120412_01/Data_img_01_20120412_01_001.mat';
% [mdata,lay] = echo_plot(fn,'CSARP_post/layerData');
%
% Author: John Paden
%
% See also: load_L1B.m, elevation_compensation.m


layer_names = {'surface','layer_01', 'layer_02','bottom'};
if ispc
  echo_fn = 'X:\ct_data\rds\2018_Antarctica_DC8\CSARP_post\CSARP_standard\20181016_01\Data_20181016_01_021.mat';
else
  echo_fn = '/cresis/snfs1/dataproducts/ct_data/rds/2018_Antarctica_DC8/CSARP_post/CSARP_standard/20181016_01/Data_20181016_01_021.mat';
end
% layerdata_source = 'layer';
mdata = load_L1B(echo_fn);
layerdata_source = echo_param(mdata);
master = [];
master.GPS_time = mdata.GPS_time;
master.Elevation = mdata.Elevation;
master.Latitude = mdata.Latitude;
master.Longitude = mdata.Longitude;
layer_params = [];

for idx = 1:length(layer_names)
  layer_params(idx).name = {'surface','layer_01','layer_02'};
  layer_params(idx).source = 'layerdata';
  layer_params(idx).layerdata_source = 'layer';
end
global gRadar

if isstruct(layerdata_source)
  param = echo_param(mdata);
  
  param = merge_structs(param,gRadar);
  
  fprintf('opsLoadLayers %s\n', param.day_seg);
  [layers] = opsLoadLayers(param,layer_params);
  layers_twtt = {};
  for lay_idx = 1:length(layer_names)
    ops_layer = [];
    ops_layer{1}.gps_time = layers(lay_idx).gps_time;
    ops_layer{1}.type = layers(lay_idx).type;
    ops_layer{1}.quality = layers(lay_idx).quality;
    ops_layer{1}.twtt = layers(lay_idx).twtt;
    ops_layer{1}.type(isnan(ops_layer{1}.type)) = 2;
    ops_layer{1}.quality(isnan(ops_layer{1}.quality)) = 1;
    lay = opsInterpLayersToMasterGPSTime(master,ops_layer,[300 60]);
    layers_twtt{lay_idx} = lay.layerData{1}.value{2}.data;
  end
else
  if ~exist('layer_names','var') || isempty(layer_names)
    layer_names = {'surface','bottom'};
  end
  
  [layers_twtt{1:length(layer_names)}] = layerdata.load_layers(mdata,layerdata_source,layer_names{:});
  
end

% mdata.Data: Nt by Nx matrix, Nt rows, Nx columns, t for fast-time, x for along-track

% nargin, Using {:} to pass arguments

% if exist('layerdata_source','var')
%   % Replace L1B echogram surface and bottom with L2 surface and bottom from layer data file
%   [mdata.Surface,mdata.Bottom] = layerdata.load_layers(mdata,layerdata_source,'surface','bottom');
%   % A = {layerdata.load_layers(mdata,layerdata_source,third_arg{:})}
% end

% if ~exist('defaultlayer','var') || isempty(defaultlayer)
%   defaultlayer = {'surface','bottom'};
% end

% mdata.Data: Nt by Nx matrix, Nt rows, Nx columns, t for fast-time, x for along-track

% nargin, Using {:} to pass arguments

% if exist('layerdata_source','var')
%   % Replace L1B echogram surface and bottom with L2 surface and bottom from layer data file
%   [mdata.Surface,mdata.Bottom] = layerdata.load_layers(mdata,layerdata_source,'surface','bottom');
%   % A = {layerdata.load_layers(mdata,layerdata_source,third_arg{:})}
% end



% if iscell(defaultlayer)
%   % Replace L1B echogram surface and bottom with L2 surface and bottom from layer data file
%   [mdata.Surface,mdata.Bottom] = layerdata.load_layers(mdata,layerdata_source,defaultlayer{1}, defaultlayer{length(defaultlayer)});
%
%
%
%
%   %   if length(defaultlayer) > 2
%   %   for i = 2:(length(defaultlayer)-1)
%   %     mdata.xlayer{i} = layerdata.load_layers(mdata,layerdata_source,defaultlayer{i});
%   %   end
%   %   end
%   % A = {layerdata.load_layers(mdata,layerdata_source,third_arg{:})}
% end
%% Set which bins to plot
param.ylims_bins = [-inf inf];
good_bins = round(max(1,min(mdata.Surface)+param.ylims_bins(1)) : min(max(mdata.Surface)+param.ylims_bins(2),size(mdata.Data,1)));

figure(1); clf;
imagesc([],mdata.Time(good_bins)*1e6,10*log10(mdata.Data(good_bins,:)))
xlabel('Range line');
ylabel('Two way travel time (us)')
colormap(1-gray(256))
hold on
plot(mdata.Surface*1e6);
for layer_idx = 2:length(layers_twtt)
  plot(layers_twtt{layer_idx}*1e6);
end
hold off;

%% Elevation Correction Example
param = [];
param.update_surf = false;
param.filter_surf = true;
if 1
  param.er_ice = 3.15;
  param.elev_comp = 3;
  %param.depth = '[min(Surface_Elev)-20 max(Surface_Elev)+2]';
  [mdata_WGS84,depth_good_idxs] = elevation_compensation(mdata,param,layers_twtt); % PADEN: ADD third argument here (matrix of TWTT, Nlayer by Nx)
  
  %% Plot versus range
  figure(2); clf;
  imagesc([],mdata_WGS84.Elevation(1) - mdata_WGS84.Elevation_Fasttime(depth_good_idxs),10*log10(mdata_WGS84.Data(depth_good_idxs,:)));
  xlabel('Range line');
  ylabel('Range (m)')
  colormap(1-gray(256))
  hold on
  plot(mdata_WGS84.Elevation(1) - mdata_WGS84.Surface_Elev);
  for layer_idx = 2:length(layers_twtt)
    plot(mdata_WGS84.Elevation(1) - mdata_WGS84.layers_elev{layer_idx});
  end
  hold off;
  
  %% Plot versus WGS-84 elevation
  figure(3); clf;
  imagesc([],mdata_WGS84.Elevation_Fasttime(depth_good_idxs),10*log10(mdata_WGS84.Data(depth_good_idxs,:)));
  xlabel('Range line');
  ylabel('WGS-84 (m)')
  set(gca,'YDir','normal')
  colormap(1-gray(256))
  hold on
  plot(mdata_WGS84.Surface_Elev);
  for layer_idx = 2:length(layers_twtt)
    plot(mdata_WGS84.layers_elev{layer_idx});
  end
  hold off;
  
else
  param.er_ice = linspace(1.5,3.15,101);
  param.er_depth = linspace(0,100,101);
  param.elev_comp = 2;
  param.depth = '[-2 20]';
  [mdata_WGS84,depth_axis] = elevation_compensation(mdata,param);
  
  %% Plot versus depth
  figure(2); clf;
  imagesc([],depth_axis,10*log10(mdata_WGS84.Data));
  xlabel('Range line');
  ylabel('Depth (m)')
  colormap(1-gray(256))
end