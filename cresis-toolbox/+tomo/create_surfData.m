function create_surfData(param,mdata)
% tomo.create_surfData(param,mdata)
%
% Description: Usually this function is called from tomo_collate_task.
%   Using a surface DEM and an ice mask, this function adds an aligned
%   surface dem and ice mask to a file.
%
% Inputs:
%   param: struct from parameter spreadsheet
%    .tomo_collate
%     .geotiff_fn
%     .ice_mask_fn
%   mdata: 3D data file struct
%     .Latitude
%     .Longitude
%     .Topography
%
% Outputs:
%   NONE
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfData,
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

%% Query OPS for surface and bottom information
param_load_layers = param;
param_load_layers.cmd.frms = round([-1,0,1] + param.proc.frm);

if ~isfield(param.tomo_collate,'layer_source') || isempty(param.tomo_collate.layer_source)
  layer_params = [];
  idx = 0;
  idx = idx + 1;
  layer_params(idx).name = 'surface';
  layer_params(idx).source = 'ops';
  idx = idx + 1;
  layer_params(idx).name = 'bottom';
  layer_params(idx).source = 'ops';
  layers = opsLoadLayers(param_load_layers,layer_params);

  %% Interpolate surface and bottom information to mdata
  master = [];
  master.GPS_time = mdata.GPS_time;
  master.Latitude = mdata.Latitude;
  master.Longitude = mdata.Longitude;
  master.Elevation = mdata.Elevation;
  for lay_idx = 1:length(layer_params)
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
elseif strcmp(param.tomo_collate.layer_source,'layerData')
  layer_source_dir = ct_filename_out(param,'','CSARP_layerData');
  layer_source_fn = sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm);
	layerData = load(fullfile(layer_source_dir,layer_source_fn));
  Surface = layerData.layerData{1}.value{2}.data;
  Bottom = layerData.layerData{2}.value{2}.data;
else
  error('Layer Source string is not recognized');
end

if length(Bottom)~=size(mdata.Topography.img,3)
  Bottom = mdata.Bottom;
end
if length(Surface)~=size(mdata.Topography.img,3)
  Surface = mdata.Surface;
end

%% Interpolate Bottom, mdata.twtt from twtt to bins
Bottom_bin = interp1(mdata.Time, 1:length(mdata.Time), Bottom);
if isfield(mdata,'ice_mask')
  ice_mask = mdata.ice_mask;
else
  ice_mask = ones(size(mdata.twtt));
end
Bottom_bin(isnan(Bottom_bin)) = -1;

%% Surface tracking prep: Convert img to double and log-scale and threshold
data = 10*log10(double(mdata.Topography.img));
% data_threshold: log scale data will be clipped to this threshold
data_threshold = param.tomo_collate.data_threshold;
data(data>data_threshold) = data_threshold;

%% Surface tracking prep: Find surface bins, truncate to end of record with guard for image template mu
mu_length = 15;
max_time = mdata.Time(end-mu_length);

twtt_bin = round(interp1(mdata.Time, 1:length(mdata.Time), mdata.twtt));
twtt_bin(isnan(twtt_bin) | twtt_bin > max_time) = length(mdata.Time)-mu_length;

%% Training parameters for image template's mu and sigma
mu = [];
sigma = [];

for rline = 1:size(mdata.Topography.img,3)
  if ~mod(rline-1,100)
    fprintf('  Training %d of %d (%s)\n', rline, size(mdata.Topography.img,3), datestr(now));
  end
  [m, s] = tomo.train_model(data(:,:,rline), ...
    double(twtt_bin(:,rline)));
  mu = [mu; m];
  sigma = [sigma; s];
end

%   mu = repmat([11.6734184113247 11.8357634315107 11.7477015213467 11.5642270542054 11.3655718245298 11.2178010788707 11.11172116154 11.0442549382899 10.9800832668574 10.9047999009164 10.8000063888223],size(mdata.Topography.img,3),1);
%   sigma = repmat([17.9297614680615 18.5178215941504 17.1485050463076 15.8106528912151 14.7936777080171 14.146975962117 13.9673485950864 13.9574519525412 13.5837122364561 13.0310380580007 12.2855990897649],size(mdata.Topography.img,3),1);

mdata.Topography.mu = mu;
mdata.Topography.sigma = sigma;

% in_dir: Directory where 3D image files are at
in_dir = ct_filename_out(param,param.tomo_collate.in_dir);

% combined_fn: Filename with 3D data
combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm));

Topography = mdata.Topography;
save(combined_fn,'-append','Topography');

%% Run detect
detect_surface = zeros(size(mdata.Topography.img,2),size(mdata.Topography.img,3));
for rline = 1:size(mdata.Topography.img,3)
  if ~mod(rline-1,500)
    fprintf('  Detect %d of %d (%s)\n', rline, size(mdata.Topography.img,3), datestr(now));
  end
  detect_surface(:,rline) = tomo.detect(data(:,:,rline), ...
    double(twtt_bin(:,rline)), ...
    double(Bottom_bin(rline)), [], double(ice_mask(:,rline)), ...
    double(mean(mdata.Topography.mu)), double(mean(mdata.Topography.sigma)), ...
    param.tomo_collate.mid,10,param.tomo_collate.smooth_weight,param.tomo_collate.smooth_var,param.tomo_collate.smooth_slope);
end

%% Run extract
if 1
  fprintf('  Extract (%s)\n', datestr(now));
  mu = [-3 -1 3 3 3 4 4 4 4 4 3 3 3 -1 -3];
  sigma = 10*ones(1,15);
  %mu = mean(mdata.Topography.mu);
  %sigma = mean(mdata.Topography.sigma);
  smooth_weight = -1;
  smooth_var = -1;
  smooth_slope = round(20*diff((linspace(-1,1,size(data,2))).^4));
  extract_surface = tomo.refine(data, double(twtt_bin), double(Bottom_bin), ...
    double([]), double(ice_mask), double(mu), double(sigma), ...
    smooth_weight, smooth_var, double(smooth_slope));
  
  extract_surface = reshape(extract_surface,size(mdata.Topography.img,2),size(mdata.Topography.img,3));
else
  extract_surface = detect_surface;
end
  
%% Create surfData
surf = [];
Ndoa = size(mdata.Topography.img,2);

surf(end+1).x = repmat((1:Ndoa).',[1 size(mdata.twtt,2)]);
surf(end).y = twtt_bin;
surf(end).plot_name_values = {'color','black','marker','x'};
surf(end).name = 'ice surface';
surf(end).surf_layer = [];
surf(end).active_layer = 1;
surf(end).mask_layer = [];
surf(end).control_layer = 7;
surf(end).visible = true;

surf(end+1).x =  repmat((1:Ndoa).',[1 size(mdata.twtt,2)]);
surf(end).y = extract_surface;
surf(end).plot_name_values = {'color','blue','marker','^'};
surf(end).name = 'bottom';
surf(end).surf_layer = 1;

surf(end).active_layer = 2;
surf(end).mask_layer = 3;
surf(end).control_layer = 4;
surf(end).visible = true;

surf(end+1).x = repmat((1:Ndoa).',[1 size(mdata.twtt,2)]);
surf(end).y = mdata.ice_mask;
surf(end).plot_name_values = {'color','white','marker','x'};
surf(end).name = 'ice mask';
surf(end).surf_layer = 1;
surf(end).active_layer = 2;
surf(end).mask_layer = 3;
surf(end).control_layer = 4;
surf(end).visible = true;

surf(end+1).x = repmat((1:Ndoa).',[1 size(mdata.twtt,2)]);
surf(end).y = NaN * zeros(size(surf(1).y));
surf(end).y(ceil(Ndoa/2)+1,:) = interp1(mdata.Time,1:length(mdata.Time),Bottom);
surf(end).plot_name_values = {'color','magenta','marker','+'};
surf(end).name = 'bottom gt';
surf(end).surf_layer = 1;
surf(end).active_layer = 2;
surf(end).mask_layer = 3;
surf(end).control_layer = 4;
surf(end).visible = true;

surf(end+1).x =  repmat((1:Ndoa).',[1 size(mdata.twtt,2)]);
surf(end).y = extract_surface;
surf(end).plot_name_values = {'color','yellow','marker','^'};
surf(end).name = 'bottom extract';
surf(end).surf_layer = 1;
surf(end).active_layer = 2;
surf(end).mask_layer = 3;
surf(end).control_layer = 4;
surf(end).visible = false;

surf(end+1).x =  repmat((1:Ndoa).',[1 size(mdata.twtt,2)]);
surf(end).y = detect_surface;
surf(end).plot_name_values = {'color','green','marker','^'};
surf(end).name = 'bottom detect';
surf(end).surf_layer = 1;
surf(end).active_layer = 2;
surf(end).mask_layer = 3;
surf(end).control_layer = 4;
surf(end).visible = false;

surf(end+1).x = repmat((1:Ndoa).',[1 size(mdata.twtt,2)]);
surf(end).y = NaN * zeros(size(surf(1).y));
surf(end).y(ceil(Ndoa/2)+1,:) = interp1(mdata.Time,1:length(mdata.Time),Surface);
surf(end).plot_name_values = {'color','magenta','marker','^'};
surf(end).name = 'surface gt';
surf(end).surf_layer = [];
surf(end).active_layer = 1;
surf(end).mask_layer = [];
surf(end).control_layer = 7;
surf(end).visible = true;

out_dir = ct_filename_out(param,param.tomo_collate.out_dir,'CSARP_surfData');
if ~isdir(out_dir)
  mkdir(out_dir);
end
out_fn_name = sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm);
out_fn = fullfile(out_dir,out_fn_name);
save(out_fn,'surf','-v7.3');

end
