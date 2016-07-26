% script imb.run_slice_browser
%
% Author: Elijah Paden, John Paden

if ~exist('run_slice_browser_init','var') || ~run_slice_browser_init
  mdata = load('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_CSA_music/20140401_03/Data_20140401_03_044.mat');
  
  %% Save layer data into a file
  layer = [];
  layer(1).x = repmat((1:64).',[1 size(mdata.twtt,2)]);
  layer(1).y = interp1(mdata.Time,1:length(mdata.Time),mdata.twtt);
  layer(1).plot_name_values = {'color','black','marker','x'};
  layer(1).name = 'surface';
  
  layer(2).x =  repmat((1:64).',[1 size(mdata.twtt,2)]);
  layer(2).y = mdata.bottom_surface.';
  layer(2).plot_name_values = {'color','blue','marker','^'};
  layer(2).name = 'bottom';
  
  layer(3).x = repmat((1:64).',[1 size(mdata.twtt,2)]);
  layer(3).y = NaN * zeros(size(layer(1).y));
  layer(3).plot_name_values = {'color','magenta','marker','+'};
  layer(3).name = 'Bcontrol';
  
  param = [];
  param.layer_fn = 'layer_data.mat';
  %if ~exist(param.layer_fn,'file')
    save(param.layer_fn,'layer')
  %end
  run_slice_browser_init = true;
end

%% Call slice_browser
try; delete(obj); end;
h_control_fig = figure(1); clf;
h_control_axes = axes('Parent',h_control_fig);
h_control_image = imagesc(lp(squeeze(mdata.Topography.img(:,33,:))),'Parent',h_control_axes);
colormap(jet(256))
obj = imb.slice_browser(lp(mdata.Topography.img),h_control_image,param);

