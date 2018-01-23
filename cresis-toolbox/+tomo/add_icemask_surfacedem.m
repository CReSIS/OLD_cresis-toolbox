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
%     .Topography
%
% Outputs:
%   NONE
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfData,
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

%% Load Geotiff and Ice Mask

% Load Geotiff projection information
proj = geotiffinfo(param.tomo_collate.geotiff_fn);

% Load Geotiff raster
[DEM, R, ~] = geotiffread(param.tomo_collate.geotiff_fn);

% Load ice mask
if isfield(param.tomo_collate,'ice_mask_fn') && ~isempty(param.tomo_collate.ice_mask_fn)
  ice_mask_all = load(param.tomo_collate.ice_mask_fn);
end

% sv_cal_fn: steering vector calibration filename
sv_cal_fn = [];
if isfield(param.tomo_collate,'sv_cal_fn')
  sv_cal_fn = param.tomo_collate.sv_cal_fn;
end

%% Remove unused DEM data
DEM_x = R(3,1) + R(2,1)*(0:size(DEM,2)-1);
DEM_y = R(3,2) + R(1,2)*(0:size(DEM,1)-1);

[mdata.x,mdata.y] = projfwd(proj,mdata.Latitude,mdata.Longitude);

DEM_threshold = param.tomo_collate.dem_guard;
DEM = DEM(DEM_y > min(mdata.y)-DEM_threshold & DEM_y < max(mdata.y)+DEM_threshold , ...
  (DEM_x > min(mdata.x)-DEM_threshold & DEM_x < max(mdata.x)+DEM_threshold));
DEM_x = DEM_x(DEM_x > min(mdata.x)-DEM_threshold & DEM_x < max(mdata.x)+DEM_threshold);
DEM_y = DEM_y(DEM_y > min(mdata.y)-DEM_threshold & DEM_y < max(mdata.y)+DEM_threshold);

%% Identify all bad values in the DEM
DEM(DEM == param.tomo_collate.geotiff_bad_value) = NaN;

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
DEM_y_mesh= repmat(DEM_y',[1 size(DEM,2)]);

%% Update ocean values with mean sea surface if ocean mask shapefile provided
if any(any(isnan(DEM))) && ~isempty(param.tomo_collate.ocean_mask_fn)

  % Convert border points of DEM to geodetic coordinates
  DEM_lat = cell(4,1);
  DEM_lon = cell(4,1);
  [DEM_lat{1},DEM_lon{1}] = projinv(proj,DEM_x,DEM_y(1)*ones(1,length(DEM_x)));
  [DEM_lat{2},DEM_lon{2}] = projinv(proj,DEM_x,DEM_y(end)*ones(1,length(DEM_x)));
  [DEM_lat{3},DEM_lon{3}] = projinv(proj,DEM_x(1)*ones(1,length(DEM_y)),DEM_y);
  [DEM_lat{4},DEM_lon{4}] = projinv(proj,DEM_x(end)*ones(1,length(DEM_y)),DEM_y);
  % Find geodetic bounding box of DEM from border points
  DEM_lims = [min([DEM_lat{:}]),min([DEM_lon{:}]);max([DEM_lat{:}]),max([DEM_lon{:}])];
  
  % Load ocean mask shape file
  ocean_shp = shaperead(param.tomo_collate.ocean_mask_fn);
  
  % Mask of polygons within DEM geodetic bounding box
  OC_bb_geo = [ocean_shp(:).BoundingBox];
  bb_mask = ~(OC_bb_geo(1,2:2:end)>DEM_lims(2,1) | OC_bb_geo(2,2:2:end)<DEM_lims(1,1) ...
    | OC_bb_geo(1,1:2:end)>DEM_lims(2,2) | OC_bb_geo(2,1:2:end)<DEM_lims(1,2));
  
  % If any polygons are in the DEM limits
  if any(bb_mask)    
    % Further restrict the polygons by checking for bounding box overlap in
    % projected coordinates
    ocean_shp = ocean_shp(bb_mask);
    poly_x = cell(0);
    poly_y = cell(0);
    for shp_idx = 1:length(ocean_shp)
      % convert polygon to projected coordinates
      [x,y] = projfwd(proj,ocean_shp(shp_idx).Y,ocean_shp(shp_idx).X);
      % if polygon is within projected bounding box
      if min(x)<max(DEM_x) && max(x)>min(DEM_x) ...
          && min(y)<max(DEM_y) && max(y)>min(DEM_y)
        % add polygon
        poly_x{end+1} = [x,nan];
        poly_y{end+1} = [y,nan];
      end
    end
    
    % Create ocean mask to determine which points like in the ocean
    if 1
      ocean_mask = true(size(DEM));
    else
      % Only update NaN points
      ocean_mask = isnan(DEM);
    end
    for poly_idx = 1:length(poly_x)
      % find bounding box of polygon
      poly_bb = [min([poly_x{poly_idx}]),min([poly_y{poly_idx}]);max([poly_x{poly_idx}]),max([poly_y{poly_idx}])];
      
      % Mask showing which DEM points are within the polygon bounding box
      if 1
        DEM_bb_mask = find((~(DEM_x_mesh<poly_bb(1,1) | DEM_x_mesh>poly_bb(2,1) ...
          | DEM_y_mesh<poly_bb(1,2) | DEM_y_mesh>poly_bb(2,2))));
      else
        % Only update NaN points
        DEM_bb_mask = find((~(DEM_x_mesh<poly_bb(1,1) | DEM_x_mesh>poly_bb(2,1) ...
          | DEM_y_mesh<poly_bb(1,2) | DEM_y_mesh>poly_bb(2,2))) & isnan(DEM));
      end
      
      % Mask showing which DEM points are in polygon (on land)
      land_mask_tmp = inpolygon(DEM_x_mesh(DEM_bb_mask),DEM_y_mesh(DEM_bb_mask),[poly_x{poly_idx}],[poly_y{poly_idx}]);
      ocean_mask(DEM_bb_mask(land_mask_tmp)) = false;
    end
    
    % If any points lie in the ocean, update these points with mean sea
    % level.
    if any(ocean_mask(:))
      % Load global sea surface
      sea_surface.fn = ct_filename_gis([],fullfile('world','dtu_meansealevel','DTU10MSS_1min.nc'));
      sea_surface.lat = ncread(sea_surface.fn,'lat');
      sea_surface.lon = ncread(sea_surface.fn,'lon');
      sea_surface.elev = ncread(sea_surface.fn,'mss').';
      
      % Convert DEM coordinates to geodetic coordinates
      [DEM_lat_mesh,DEM_lon_mesh] = projinv(proj,DEM_x_mesh,DEM_y_mesh);
      % Interpolate sea elevation to ocean NaN points
      sea_elev = interp2(sea_surface.lat,sea_surface.lon-180,sea_surface.elev.',DEM_lat_mesh(ocean_mask),DEM_lon_mesh(ocean_mask));
      % Replace NaN values in DEM with mean sea surface values
      DEM(ocean_mask) = sea_elev;
    end
    
    % Debug plots
    if 0
      figure; clf;
      imagesc(DEM_x,DEM_y,DEM);
      hold on;
      for idx = 1:length(poly_x)
        plot(poly_x{idx},poly_y{idx},'r');
      end
      plot(mdata.x,mdata.y,'k');
      set(gca,'YDir','normal');
      keyboard
    end
  end
  
end

%% Create a point cloud from the DEM
if 0
  figure(1); clf;
  h_img = imagesc(DEM_x,DEM_y, DEM);
  set(gca,'YDir','normal');
  hold on;
  plot(mdata.x,mdata.y,'k.','LineWidth',2);
  plot(mdata.x(1),mdata.y(1),'ro','LineWidth',2);
  imagesc(OM_x,OM_y,max(max(DEM))*ones(length(OM_y),length(OM_x)));
  hold off;
end

%% First slice
Nx = length(mdata.GPS_time);

theta = mdata.param_combine.array_param.theta;
if ~isempty(sv_cal_fn)
  theta_cal = load(sv_cal_fn);
  theta = theta_cal.theta;
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
  if ~mod(rline-1,500)
    fprintf('  Ice-DEM-Mask %d of %d (%s)\n', rline, Nx, datestr(now));
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
  
  origin = mdata.param_combine.array_param.fcs{1}{1}.origin(:,rline);
  
  % Convert from ECEF to FCS/SAR
  Tfcs_ecef = [mdata.param_combine.array_param.fcs{1}{1}.x(:,rline), ...
    mdata.param_combine.array_param.fcs{1}{1}.y(:,rline), ...
    mdata.param_combine.array_param.fcs{1}{1}.z(:,rline)];
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
  
  if(numel(x)>=3)
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
    imagesc([],mdata.Time,lp(mdata.Topography.img(:,:,rline)))
    hold on
    plot(twtt,'k')
    hold off
    pause
  end
  
end

if 0
  %% DEBUG
  clf
  imagesc([],theta,ice_mask);
end

%% Add DEM ground truth and ice mask to the echogram files

% in_dir: Directory where 3D image files are at
in_dir = ct_filename_out(param,param.tomo_collate.in_dir);
    
% combined_fn: Filename with 3D data
combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm));

ice_mask = logical(ice_mask);
save(combined_fn,'-append','twtt','ice_mask','theta');

if exist('theta_cal','var')
  save(combined_fn,'-v7.3','-append','theta_cal');
else
  theta_cal = theta;
  save(combined_fn,'-v7.3','-append','theta_cal');
end

mdata.theta = theta;
mdata.theta_cal = theta_cal;
mdata.twtt = twtt;
mdata.ice_mask = ice_mask;

end
