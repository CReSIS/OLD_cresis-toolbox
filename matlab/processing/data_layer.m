function [hdr] = data_layer(param,hdr,field_name)
% [hdr] = data_layer(param,hdr,field_name)
%
% https://ops.cresis.ku.edu/wiki/index.php/Data_load#data_layer.m
%
% Author: John Paden

if ~isempty(param.load.layer)
  %% Load layer
  layer_param = param.load.layer;
  layer_param.existence_check = false;
  % opsLoadLayers will eliminate invalid frames due to edge conditions for
  % frame 1 and the end frame
  layer_param.frms = param.load.frm-1 : param.load.frm+1;
  layer = opsLoadLayers(param,layer_param);
  
  %% Interpolate layer onto data gps time axes
  master = [];
  master.GPS_time = hdr.gps_time;
  master.Latitude = hdr.records{1}.lat;
  master.Longitude = hdr.records{1}.lon;
  master.Elevation = hdr.records{1}.elev;
  ops_layer = [];
  ops_layer{1}.gps_time = layer.gps_time;
  ops_layer{1}.type = layer.type;
  ops_layer{1}.quality = layer.quality;
  ops_layer{1}.twtt = layer.twtt;
  ops_layer{1}.type(isnan(ops_layer{1}.type)) = 2;
  ops_layer{1}.quality(isnan(ops_layer{1}.quality)) = 1;
  lay = opsInterpLayersToMasterGPSTime(master,ops_layer,[300 60]);
  hdr.(field_name) = lay.layerData{1}.value{2}.data;
end
