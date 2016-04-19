function [elev,latm,lonm,elevm] = latlon2elevation_profile(lat,lon,use_global_vars,map_interp,plot_en)
% [elev,latm,lonm,elevm] = latlon2elevation_profile(lat,lon,use_global_vars,map_interp,plot_en)
%
% Function for creating an elevation profile from a set of waypoints.
%
% lat: latitude (deg,N)
% lon: longitude (deg,E), same dimensions as lat
% use_global_vars: logical (true will use global variables sea_surface and
%   land_surface to avoid having to reload the large DEMs)
% map_interp: angular sampling rate for map-interpolated outputs (deg)
%
% elev: elevation profile corresponding to lat,lon
% latm,lonm,elevm: interpm interpolated version of elev [requires mapping
%   toolbox]
%
% Author: John Paden

if use_global_vars
  global sea_surface;
  global land_surface;
end

%% Load mean sea level
if ~exist('sea_surface','var') || isempty(sea_surface)
  sea_surface.fn = ct_filename_gis([],fullfile('world','dtu_meansealevel','DTU10MSS_1min.nc'));
  if exist(sea_surface.fn,'file')
    sea_surface.lat = ncread(sea_surface.fn,'lat');
    sea_surface.lon = ncread(sea_surface.fn,'lon');
    sea_surface.elev = ncread(sea_surface.fn,'mss').';
  else
    warning('DTU sea level does not exist: %s\n', sea_surface.fn);
    sea_surface.fn = ct_filename_gis([],fullfile('world','egm96_geoid','WW15MGH.DAC'));
    if ~exist(sea_surface.fn,'file')
      warning('EGM96 geoid does not exist: %s\n', sea_surface.fn);
      sea_surface = [];
    else
      [sea_surface.lat,sea_surface.lon,sea_surface.elev] = egm96_loader(sea_surface.fn);
    end
  end
end

%% Load Greenland or Antarctic DEM
if ~exist('land_surface','var') || isempty(land_surface)
  if lat(1) > 0
    land_surface.fn = ct_filename_gis([],'greenland/DEM/GIMP/gimpdem_90m.tif');
  else
    land_surface.fn = ct_filename_gis([],'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_surface.tif');
  end
  if ~exist(land_surface.fn)
    warning('Land DEM does not exist: %s\n', land_surface.fn);
    land_surface = [];
  else
    [land_surface.dem, land_surface.R, tmp] = geotiffread(land_surface.fn);
    land_surface.dem = double(land_surface.dem);
    land_surface.dem(land_surface.dem == 32767) = NaN;
    land_surface.proj = geotiffinfo(land_surface.fn);
    if 0
      % Debug Plot
      figure; clf;
      land_surface.x = land_surface.R(3,1) + land_surface.R(2,1)*(0:size(land_surface.dem,2)-1);
      land_surface.y = land_surface.R(3,2) + land_surface.R(1,2)*(0:size(land_surface.dem,1)-1);
      imagesc(land_surface.x,land_surface.y,land_surface.dem)
      set(gca,'YDir','normal');
    end
  end
end

%% Interpolate Land DEM and Geoid/mean-sea-surface
if ~isempty(sea_surface)
  sea_dem = interp2(sea_surface.lon,sea_surface.lat,sea_surface.elev,mod(lon,360),lat);
else
  sea_dem = zeros(size(lat));
end
if ~isempty(land_surface)
  [x,y] = projfwd(land_surface.proj,lat,lon);
  land_dem = interp2(land_surface.dem,(x-land_surface.R(3,1))/land_surface.R(2,1)+1, ...
    (y-land_surface.R(3,2))/land_surface.R(1,2)+1);
else
  land_dem = zeros(size(lat));
end

% Merge GIMP and Geoid
elev = land_dem; % Use land DEM first
elev(isnan(elev)) = sea_dem(isnan(elev)); % Fill in NaN with mean sea level

%% Create map interpolated profile
[latm,lonm] = interpm(lat,lon,map_interp);

if ~isempty(sea_surface)
  sea_dem = interp2(sea_surface.lon,sea_surface.lat,sea_surface.elev,mod(lonm,360),latm);
else
  sea_dem = zeros(size(latm));
end
if ~isempty(land_surface)
  [x,y] = projfwd(land_surface.proj,latm,lonm);
  land_dem = interp2(land_surface.dem,(x-land_surface.R(3,1))/land_surface.R(2,1)+1, ...
    (y-land_surface.R(3,2))/land_surface.R(1,2)+1);
else
  land_dem = zeros(size(latm));
end

% Merge GIMP and Geoid
elevm = land_dem; % Use land DEM first
elevm(isnan(elevm)) = sea_dem(isnan(elevm)); % Fill in NaN with mean sea level

%% Debug Plots
if plot_en
  figure; clf;
  scatter(lon, lat, repmat(10,size(lat)), elev, 'Marker', 'x');
  hold on;
  grid on;
  scatter(lonm, latm, repmat(3,size(latm)), elevm(:), 'Marker', '.');
  colormap(jet(256));
  h = colorbar;
  set(get(h,'YLabel'), 'String', 'Elevation (m)');
  
  figure; clf;
  along_trackm = geodetic_to_along_track(latm,lonm);
  plot(along_trackm/1e3, elevm*100/2.54/12);
  grid on;
  xlabel('Along track (km)');
  ylabel('Elevation (ft)');
  
  if ~isempty(land_dem)
    figure; clf;
    land_surface.x = land_surface.R(3,1) + land_surface.R(2,1)*(0:size(land_surface.dem,2)-1);
    land_surface.y = land_surface.R(3,2) + land_surface.R(1,2)*(0:size(land_surface.dem,1)-1);
    buffer_m = 2500; % Buffer around track
    xidxs = land_surface.x > min(x)-buffer_m & land_surface.x < max(x)+buffer_m;
    yidxs = land_surface.y > min(y)-buffer_m & land_surface.y < max(y)+buffer_m;
    imagesc(land_surface.x(xidxs)/1e3,land_surface.y(yidxs)/1e3,land_surface.dem(yidxs,xidxs) * 100/2.54/12);
    xlabel('X (km)');
    ylabel('Y (km)');
    set(gca,'YDir','normal');
    hold on;
    plot(x/1e3,y/1e3,'k');
    h=colorbar;
    set(get(h,'YLabel'),'String','Elevation (ft)');
  end
end

return
