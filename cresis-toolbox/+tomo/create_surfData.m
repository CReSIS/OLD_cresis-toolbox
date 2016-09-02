function create_surfData(param,mdata)
% tomo.create_surfData(param,mdata)
%
% Description: Usually this function is called from tomo.collate_task.
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

%% Convert Surface and Bottom variables from propagation time into image pixels
Bottom_bin = interp1(mdata.Time, 1:length(mdata.Time), mdata.Bottom);
twtt_bin = interp1(mdata.Time, 1:length(mdata.Time), mdata.twtt);
ice_mask = mdata.ice_mask;
Bottom_bin(isnan(Bottom_bin)) = -1;
twtt_bin(isnan(twtt_bin)) = -1;

%% Training parameters mu and sigma
mu = [];
sigma = [];
for rline = 1:size(mdata.Topography.img,3)
%   if ~mod(rline-1,100)
%     fprintf('  Training %d of %d (%s)\n', rline, size(mdata.Topography.img,3), datestr(now));
%   end
  [m, s] = tomo.train_params(10*log10(double(mdata.Topography.img(:,:,rline))), ...
    double(twtt_bin(:,rline)));
  mu = [mu; m];
  sigma = [sigma; s];
end

mdata.Topography.mu = mu;
mdata.Topography.sigma = sigma;

% in_dir: Directory where 3D image files are at
in_dir = ct_filename_out(param,param.tomo_collate.in_dir);

% combined_fn: Filename with 3D data
combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm));

Topography = mdata.Topography;
save(combined_fn,'-append','Topography');

%% Surface tracking prep: Convert img to double and log-scale and threshold
data = 10*log10(double(mdata.Topography.img));
% data_threshold: log scale data will be clipped to this threshold
data_threshold = param.tomo_collate.data_threshold;
data(data>data_threshold) = data_threshold;

%% Run detect
detect_surface = zeros(size(mdata.Topography.img,2),size(mdata.Topography.img,3));
for rline = 1:size(mdata.Topography.img,3)
  if ~mod(rline-1,500)
    fprintf('  Detect %d of %d (%s)\n', rline, size(mdata.Topography.img,3), datestr(now));
  end
  detect_surface(:,rline) = tomo.detect(data(:,:,rline), ...
    double(twtt_bin(:,rline)), ...
    double(Bottom_bin(rline)), [], double(ice_mask(:,rline)), double(mean(mdata.Topography.mu)), double(mean(mdata.Topography.sigma)));
end

%% Run extract
fprintf('  Extract (%s)\n', datestr(now));
extract_surface = tomo.extract(data, double(twtt_bin), double(Bottom_bin), ...
  double([]), double(ice_mask), double(mean(mdata.Topography.mu)), double(mean(mdata.Topography.sigma)));
extract_surface = reshape(extract_surface,size(mdata.Topography.img,2),size(mdata.Topography.img,3));

%% Create surfData
surf = [];
Ndoa = size(mdata.Topography.img,2);

surf(end+1).x = repmat((1:Ndoa).',[1 size(mdata.twtt,2)]);
surf(end).y = interp1(mdata.Time,1:length(mdata.Time),mdata.twtt);
surf(end).plot_name_values = {'color','black','marker','x'};
surf(end).name = 'ice surface';
surf(end).surf_layer = [];
surf(end).active_layer = 1;
surf(end).mask_layer = [];
surf(end).control_layer = [];
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
surf(end).y(ceil(Ndoa/2)+1,:) = interp1(mdata.Time,1:length(mdata.Time),mdata.Bottom);
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

out_dir = ct_filename_out(param,param.tomo_collate.out_dir,'CSARP_surfData');
if ~isdir(out_dir)
  mkdir(out_dir);
end
out_fn_name = sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm);
out_fn = fullfile(out_dir,out_fn_name);
save(out_fn,'surf','-v7.3');

end
