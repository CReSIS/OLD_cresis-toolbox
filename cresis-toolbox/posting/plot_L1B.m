% script plot_L1B
%
% Example of loading data with load_L1B.m and using elevation_compensation.m
%
% Plots data in three different ways:
%   Figure 1: time-delay on y-axis
%   Figure 2: range on y-axis
%   Figure 3: time-delay on y-axis
%
% Author: John Paden
%
% See also: load_L1B.m, elevation_compensation.m

%  fn = 'IRMCR1B_V01_20130408_01_020.nc';
%  mdata = load_L1B(fn);

fn = '/cresis/snfs1/dataproducts/ct_data/snow/2014_Greenland_P3/CSARP_post/CSARP_qlook/20140324_01/Data_20140324_01_010.mat';
mdata = load_L1B(fn);
if 0
  % Replace L1B echogram surface with L2 surface from layer data file
  fn = '/cresis/snfs1/dataproducts/ct_data/snow/2014_Greenland_P3/CSARP_post/CSARP_layerData/20140324_01/Data_20140324_01_010.mat';
  lay = load(fn);
  mdata.Surface = interp1(lay.GPS_time,lay.layerData{1}.value{2}.data,mdata.GPS_time);
end

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
hold off;

%% Elevation Correction Example
param = [];
param.update_surf = false;
param.filter_surf = true;
if 1
  param.er_ice = 3.15;
  param.elev_comp = 3;
  %param.depth = '[min(Surface_Elev)-20 max(Surface_Elev)+2]';
  [mdata_WGS84,depth_good_idxs] = elevation_compensation(mdata,param);
  
  %% Plot versus range
  figure(2); clf;
  imagesc([],mdata_WGS84.Elevation(1) - mdata_WGS84.Elevation_Fasttime(depth_good_idxs),10*log10(mdata_WGS84.Data(depth_good_idxs,:)));
  xlabel('Range line');
  ylabel('Range (m)')
  colormap(1-gray(256))
  hold on
  plot(mdata_WGS84.Elevation(1) - mdata_WGS84.Surface_Elev);
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

return;
