function mdata = add_icemask_surfacedem(param, mdata)
% mdata = tomo.add_icemask_surfacedem(param, mdata)
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
%     .Tomo
%
% Outputs:
%   NONE
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfData,
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

%% Load Geotiff and Ice Mask
dem_res = 10;
global gdem;
if isempty(gdem) || ~ishandle(gdem) || ~isvalid(gdem)
  gdem = dem_class(param,dem_res);
end
gdem.set_res(dem_res);

% Load ice mask
if isfield(param.tomo_collate,'ice_mask_fn') && ~isempty(param.tomo_collate.ice_mask_fn)
  ice_mask_fn = ct_filename_gis(param,param.tomo_collate.ice_mask_fn);
  [~,ice_mask_fn_name,ice_mask_fn_ext] = fileparts(ice_mask_fn);
  if strcmpi(ice_mask_fn_ext,'.tif')
    ice_mask_all.proj = geotiffinfo(ice_mask_fn);
    [ice_mask_all.mask, R, ~] = geotiffread(ice_mask_fn);
    ice_mask_all.X = R(3,1) + R(2,1)*(1:size(ice_mask_all.mask,2));
    ice_mask_all.Y = R(3,2) + R(1,2)*(1:size(ice_mask_all.mask,1));
  else
    ice_mask_all = load(ice_mask_fn);
  end
else
  ice_mask_all = [];
end

% sv_cal_fn: steering vector calibration filename
sv_cal_fn = [ct_filename_ct_tmp(param,'','tomo_collate','theta_cal') '.mat'];
if ~exist(sv_cal_fn,'file')
  sv_cal_fn = [];
end

%% Create DEM
Nx = length(mdata.GPS_time);
dec_idxs = round(linspace(1,Nx,min(Nx,200)));
physical_constants;
[latb,lonb] = bufferm(mdata.Latitude(dec_idxs),mdata.Longitude(dec_idxs),param.tomo_collate.dem_guard/WGS84.semimajor*180/pi);
gdem_str = sprintf('%s:%s:%s_%03d',param.radar_name,param.season_name,param.day_seg,param.load.frm);
if ~strcmpi(gdem_str,gdem.name)
  gdem.set_vector(latb,lonb,gdem_str);
end

gdem.set_vector(latb,lonb,gdem_str);
[DEM,msl,ocean_mask,proj,DEM_x,DEM_y] = gdem.get_vector_mosaic(100);
DEM(ocean_mask) = msl(ocean_mask);
[mdata.x,mdata.y] = projfwd(proj,mdata.Latitude,mdata.Longitude);

%% Interpolate at all the bad value locations using the good data
bad_idxs = find(isnan(DEM));
good_idxs = find(~isnan(DEM));
x_idxs = repmat(1:size(DEM,2),[size(DEM,1) 1]);
y_idxs = repmat((1:size(DEM,1))',[1 size(DEM,2)]);
x_vals = x_idxs(good_idxs);
y_vals = y_idxs(good_idxs);
z_vals = DEM(good_idxs);
x_out = x_idxs(bad_idxs);
y_out = y_idxs(bad_idxs);
z_out = single(griddata(x_vals,y_vals,double(z_vals),x_out,y_out));
if ~isempty(z_out)
  DEM(bad_idxs) = z_out;
end

if 0
  figure(1);clf;
  imagesc(DEM_x,DEM_y,DEM);
  hold on;
  plot(mdata.x,mdata.y,'r');
  keyboard;
end

DEM_x_mesh = repmat(DEM_x,[size(DEM,1) 1]);
DEM_y_mesh= repmat(DEM_y,[1 size(DEM,2)]);

%% First slice
Nx = length(mdata.GPS_time);

theta = mdata.Tomo.theta(:,1); % Theta is constant in each column
if ~isempty(sv_cal_fn)
  theta_cal = load(sv_cal_fn);
  theta = interp1(theta_cal.theta_original, theta_cal.theta, theta, 'linear', 'extrap');
  theta_cal = theta;
end

Nsv = length(theta);
twtt = zeros(Nsv,Nx);
ice_mask = ones(Nsv,Nx);

if all(all(isnan(DEM)))
  warning('Input DEM contains all NaN data for Frame %d.',param.proc.frm);
  twtt(:,:) = NaN;
  Nx = 0;
end
  
DEM_coverage_warning = false;
for rline = 1:Nx
  if ~mod(rline-1,10^floor(log10(Nx)-1))
    fprintf('  %s %d of %d (%s)\n', mfilename, rline, Nx, datestr(now));
  end
  
  dem_guard = param.tomo_collate.dem_guard;
  DEM_mask = DEM_x_mesh > mdata.x(rline)-dem_guard & DEM_x_mesh < mdata.x(rline)+dem_guard ...
    & DEM_y_mesh > mdata.y(rline)-dem_guard & DEM_y_mesh < mdata.y(rline)+dem_guard ...
    & ~isnan(DEM);
  DEM_idxs = find(DEM_mask);
  
  if numel(DEM_idxs)==0
    warning('Range Line %d of Frame %d is not spanned by DEM.',rline,param.proc.frm);
  end
  
  if 0
    set(h_img,'AlphaData',DEM_mask);
  end
  
  % Convert from projection to geodetic (lat,lon,elev)
  [DEM_lat,DEM_lon] = projinv(proj,DEM_x_mesh(DEM_idxs),DEM_y_mesh(DEM_idxs));
  DEM_elev = DEM(DEM_idxs);
  
  % Convert from geodetic (lat,lon,elev) to ECEF (x,y,z)
  physical_constants;
  [DEM_ecef_x,DEM_ecef_y,DEM_ecef_z] = geodetic2ecef(single(DEM_lat)/180*pi,single(DEM_lon)/180*pi,single(DEM_elev),WGS84.ellipsoid);
  
  origin = mdata.param_array.array_proc.fcs{1}{1}.origin(:,rline);
  
  % Convert from ECEF to FCS/SAR
  Tfcs_ecef = [mdata.param_array.array_proc.fcs{1}{1}.x(:,rline), ...
    mdata.param_array.array_proc.fcs{1}{1}.y(:,rline), ...
    mdata.param_array.array_proc.fcs{1}{1}.z(:,rline)];
  Tecef_fcs = inv(Tfcs_ecef);
  
  tmp = Tecef_fcs * [DEM_ecef_x.'-origin(1); DEM_ecef_y.'-origin(2); DEM_ecef_z.'-origin(3)];
  DEM_fcs_x = tmp(1,:);
  DEM_fcs_y = tmp(2,:);
  DEM_fcs_z = tmp(3,:);
  
  if 0
    imagesc(reshape(DEM_fcs_x,[200 200]))
    colorbar;
    
    imagesc(reshape(DEM_fcs_y,[200 200]))
    colorbar;
    
    imagesc(reshape(DEM_fcs_z,[200 200]))
    colorbar;
  end
  
  slice_mask = DEM_fcs_x > -param.tomo_collate.dem_per_slice_guard & DEM_fcs_x < param.tomo_collate.dem_per_slice_guard;
  
  x = DEM_fcs_x(slice_mask);
  y = DEM_fcs_y(slice_mask);
  z = DEM_fcs_z(slice_mask);
  
  if (numel(x)>=3)
    faces = delaunay(double(x),double(y));
    vertices = [double(x).' double(y).' double(z).'];  % vertices stored as Nx3 matrix
    vert1 = vertices(faces(:,1),:);
    vert2 = vertices(faces(:,2),:);
    vert3 = vertices(faces(:,3),:);
  
    orig = [0 0 0];

    intersection = zeros(3,Nsv);

    for theta_idx = 1:length(theta)
      dir = [0 sin(theta(theta_idx)) -cos(theta(theta_idx))];
      [intersect, t] = TriangleRayIntersection(orig, dir, vert1, vert2, vert3);

      intersect_idx = find(intersect);

      if isempty(intersect_idx)
        twtt(theta_idx,rline) = NaN;
        intersection(:,theta_idx) = NaN;
      else
        twtt(theta_idx,rline) = t(intersect_idx(1))/(3e8/2);
        % finds coordinates in approximate center of triangles
        intersection(:,theta_idx) = mean([vert1(intersect_idx(1),:);vert2(intersect_idx(1),:);vert3(intersect_idx(1),:)],1);
      end
    end
  else
    if ~DEM_coverage_warning
      DEM_coverage_warning = true;
      warning('DEM dem_per_slice_guard too small.');
    end
    clear intersection;
    twtt(:,rline) = NaN;
  end

  if exist('ice_mask_all','var')
    if exist('intersection','var')
      % Convert from FCS/SAR to ECEF
      intersection_ecef = Tfcs_ecef * intersection;
      intersection_ecef_x = intersection_ecef(1,:).' + origin(1);
      intersection_ecef_y = intersection_ecef(2,:).' + origin(2);
      intersection_ecef_z = intersection_ecef(3,:).' + origin(3);
      % Convert from ECEF to geodetic
      [intersection_lat,intersection_lon,tri_h] = ecef2geodetic(intersection_ecef_x,intersection_ecef_y,intersection_ecef_z,WGS84.ellipsoid);
      intersection_lat = intersection_lat*180/pi;
      intersection_lon = intersection_lon*180/pi;
      % Convert from geodetic to projection
      [intersection_x,intersection_y] = projfwd(ice_mask_all.proj,intersection_lat,intersection_lon);
      % Get mask coordinates nearest triangle center coordinates
      intersection_x_idx = interp1(ice_mask_all.X,1:length(ice_mask_all.X),intersection_x,'nearest');
      intersection_y_idx = interp1(ice_mask_all.Y,1:length(ice_mask_all.Y),intersection_y,'nearest');
      % Find nan values and set to integer value
      nidx = find(isnan(intersection_x_idx));
      intersection_x_idx(nidx) = 1;
      intersection_y_idx(nidx) = 1;
      % Convert triangle mask coordinates to matrix indices
      mask_idx = (intersection_x_idx-1)*length(ice_mask_all.Y) + intersection_y_idx;
      % Find ice mask for triangle coordinates
      ice_mask(:,rline) = ice_mask_all.mask(mask_idx);
      % Set previously nan valued coordinates to 0 mask
      ice_mask(nidx,rline) = 0;
    end
  end
  
  if 0
    figure(1); clf;
    imagesc([],mdata.Time,lp(mdata.Tomo.img(:,:,rline)))
    hold on
    plot(twtt(:,rline),'k')
    hold off
    keyboard
  end
  
end

if 0
  %% DEBUG
  clf
  imagesc([],theta,ice_mask);
end

%% Add DEM ground truth and ice mask to the echogram files

% in_dir: Directory where 3D image files are at
in_dir = ct_filename_out(param,param.tomo_collate.in_path);
    
% combined_fn: Filename with 3D data
combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.load.frm));

ice_mask = logical(ice_mask);
if param.ct_file_lock
  file_version = '1L';
else
  file_version = '1';
end
save(combined_fn,'-append','twtt','ice_mask','theta','file_version');

if exist('theta_cal','var')
  save(combined_fn,'-append','theta_cal');
else
  theta_cal = theta;
  save(combined_fn,'-append','theta_cal');
end

mdata.theta = theta;
mdata.theta_cal = theta_cal;
mdata.twtt = twtt;
mdata.ice_mask = ice_mask;

end
