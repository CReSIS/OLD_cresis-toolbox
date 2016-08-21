function mdata = data_loader_prep(param,mdata)
% mdata = data_loader_prep(param,mdata)
%
% Description. Usually this function is called from tomo.collate. Loads
%   layer data.
%
% Inputs:
%   param = struct with processing parameters
%   mdata = contains frame data
%
% Outputs:
%   mdata = contains frame data
%
% See also: tomo.collate
%
% Author: 

fprintf('Load layer data...\n');
%% Load Layers
param_load_layers = param;
param_load_layers.cmd.frms = round([-1,0,1] + mdata{1}.frm);

layer_params = [];
idx = 0;
idx = idx + 1;
layer_params(idx).name = 'surface';
layer_params(idx).source = 'ops';
idx = idx + 1;
layer_params(idx).name = 'bottom';
layer_params(idx).source = 'ops';
layers = opsLoadLayers(param_load_layers,layer_params);

%% Interpolate the data (preserving gaps)
master = [];
master.GPS_time = mdata{1}.GPS_time;
master.Latitude = mdata{1}.Latitude;
master.Longitude = mdata{1}.Longitude;
master.Elevation = mdata{1}.Elevation;
for lay_idx = 1:2
  ops_layer = [];
  ops_layer{1}.gps_time = layers(lay_idx).gps_time;

  ops_layer{1}.type = layers(lay_idx).type;
  ops_layer{1}.quality = layers(lay_idx).quality;
  ops_layer{1}.twtt = layers(lay_idx).twtt;
  ops_layer{1}.type(isnan(ops_layer{1}.type)) = 2;
  ops_layer{1}.quality(isnan(ops_layer{1}.quality)) = 1;
  lay = opsInterpLayersToMasterGPSTime(master,ops_layer,[300 60]);
  layers(lay_idx).twtt_ref = lay.layerData{1}.value{2}.data;
end

Surface = layers(1).twtt_ref;
Bottom = layers(2).twtt_ref;

if 0
  %% Debug
  figure(1); clf;
  imagesc([],mdata{2}.Time,lp(mdata{2}.Data));
  hold on;
  h_plot = plot(Surface,'k-');
  plot(Bottom,'k-');
  hold off;
end

%% Add surface and bottom layers to the echogram files
for img=1:3
  fn_dir = ct_filename_out(param,param.surf_extract.out_dir);
  fn = fullfile(fn_dir,sprintf('Data_img_%02.0f_%s_%03.0f.mat',img, ...
    param.day_seg,mdata{img}.frm));

  save(fn,'-append','Surface','Bottom');
  mdata{img}.Surface = Surface;
  mdata{img}.Bottom = Bottom;
end

return