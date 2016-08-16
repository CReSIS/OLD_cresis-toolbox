function mdata = DEM_alignment(param,mdata)
% mdata = DEM_alignment(param,mdata)
%
% Description. Usually this function is called from tomo.collate. Estimates
%   signal TWTT time from DEM and creates an ice mask.
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
% Author: John Paden, Jordan Sprick, and Mingze Xu

  fprintf('Finding Ice Mask and Ground Truth...\n');

  %% Load Geotiff and Ice Mask

  % Get the projection information
  proj = geotiffinfo(param.surf_extract.geotiff_fn);
  % Get ice mask
  ice_mask_all = load(param.surf_extract.ice_mask_fn);
  % sv_cal_fn: steering vector calibration filename
  sv_cal_fn = ct_filename_ct_tmp(rmfield(param,'day_seg'),'','sv_calibration','theta_cal.mat');
  
  img = 1;

  %% Remove unused DEM data

  % Read the image
  [DEM, R, tmp] = geotiffread(param.surf_extract.geotiff_fn);

  DEM_x = R(3,1) + R(2,1)*(0:size(DEM,2)-1);
  DEM_y = R(3,2) + R(1,2)*(0:size(DEM,1)-1);

  [mdata{img}.x,mdata{img}.y] = projfwd(proj,mdata{img}.Latitude,mdata{img}.Longitude);

  DEM_threshold = 8e3;
  DEM = DEM(DEM_y > min(mdata{img}.y)-DEM_threshold & DEM_y < max(mdata{img}.y)+DEM_threshold , ...
    (DEM_x > min(mdata{img}.x)-DEM_threshold & DEM_x < max(mdata{img}.x)+DEM_threshold));
  DEM_x = DEM_x(DEM_x > min(mdata{img}.x)-DEM_threshold & DEM_x < max(mdata{img}.x)+DEM_threshold);
  DEM_y = DEM_y(DEM_y > min(mdata{img}.y)-DEM_threshold & DEM_y < max(mdata{img}.y)+DEM_threshold);

  %% Identify all bad values in the DEM
  DEM(DEM == -32767) = NaN;

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
  DEM(bad_idxs) = z_out;


  %% Create a point cloud from the DEM
  if 0
    figure(1); clf;
    h_img = imagesc(DEM_x,DEM_y, DEM);
    set(gca,'YDir','normal');
    hold on;
    plot(mdata{img}.x,mdata{img}.y,'k.','LineWidth',2);
    plot(mdata{img}.x(1),mdata{img}.y(1),'ro','LineWidth',2);
    hold off;
  end

  %% First slice
  Nx = length(mdata{1}.GPS_time);

  theta = mdata{1}.param_combine.array_param.theta;
  if param.surf_extract.theta_calibrated
    theta_cal = load(sv_cal_fn);
    theta = theta_cal.theta;
  end

  Nsv = length(theta);
  twtt = zeros(Nsv,Nx);
  ice_mask = zeros(Nsv,Nx);

  for rline = 1:Nx

    DEM_x_mesh = repmat(DEM_x,[size(DEM,1) 1]);
    DEM_y_mesh= repmat(DEM_y',[1 size(DEM,2)]);

    DEM_mask = DEM_x_mesh > mdata{img}.x(rline)-4e3 & DEM_x_mesh < mdata{img}.x(rline)+4e3 ...
      & DEM_y_mesh > mdata{img}.y(rline)-4e3 & DEM_y_mesh < mdata{img}.y(rline)+4e3;
    DEM_idxs = find(DEM_mask);

    if 0
      set(h_img,'AlphaData',DEM_mask);
    end

    % Convert from projection to geodetic (lat,lon,elev)
    [DEM_lat,DEM_lon] = projinv(proj,DEM_x_mesh(DEM_idxs),DEM_y_mesh(DEM_idxs));
    DEM_elev = DEM(DEM_idxs);

    % Convert from geodetic (lat,lon,elev) to ECEF (x,y,z)
    physical_constants;
    [DEM_ecef_x,DEM_ecef_y,DEM_ecef_z] = geodetic2ecef(DEM_lat/180*pi,DEM_lon/180*pi,DEM_elev,WGS84.ellipsoid);

    origin = mdata{1}.param_combine.array_param.fcs{1}{1}.origin(:,rline);

    % Convert from ECEF to FCS/SAR
    Tfcs_ecef = [mdata{1}.param_combine.array_param.fcs{1}{1}.x(:,rline), ...
      mdata{1}.param_combine.array_param.fcs{1}{1}.y(:,rline), ...
      mdata{1}.param_combine.array_param.fcs{1}{1}.z(:,rline)];
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

    slice_mask = DEM_fcs_x > -120 & DEM_fcs_x < 120;

    x = DEM_fcs_x(slice_mask);
    y = DEM_fcs_y(slice_mask);
    z = DEM_fcs_z(slice_mask);

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

    if 0
      imagesc([],mdata{1}.Time,lp(mdata{1}.Topography.img(:,:,rline)))
      hold on
      plot(twtt,'k')
      hold off
      pause
    end

  end

  if 0
    clf
    imagesc([],theta,ice_mask);
  end

  %% Add DEM ground truth and ice mask to the echogram files
  for img=1:3
    %   fn = fullfile(fn_dir,sprintf('Data_img_%02.0f_20140401_03_%03.0f.mat',img,frame));
    fn_dir = ct_filename_out(param,param.surf_extract.out_dir);
    fn = fullfile(fn_dir,sprintf('Data_img_%02.0f_%s_%03.0f.mat',img, ...
      param.day_seg,mdata{img}.frm));

    save(fn,'-append','twtt','ice_mask');
    mdata{img}.twtt = twtt;
    mdata{img}.ice_mask = ice_mask;
  end
  
  return;