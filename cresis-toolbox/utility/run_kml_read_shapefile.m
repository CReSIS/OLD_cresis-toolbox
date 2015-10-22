% function run_kml_read_shapefile
%
% Example of how to run kml_read_shapefile and plot the results over a geotiff.
%
% Author: Aric Beaver

% This script is intended to load KMZ files and plot them over a geotiff.

kml_fn = '/cresis/scratch2/abeaver/KMZ/icethick/icethick.kml';
% out_fn = fullfile('/cresis/scratch2/abeaver/KMZ/ickthick/','icethick.txt');
% geotiff_fn = '/cresis/scratch2/kpurdon/antarctica_tiff/NEW/nsidc_ant500m_wgs84_elev_m.tif';
geotiff_fn = '/cresis/scratch1/mdce/csarp_support/Landsat-7/antarctica/Antarctica_LIMA_peninsula.tif';


%%%%%%%%%%%%% Load geometry, lat, lon, and idx from KML file %%%%%%%%%%%%%%

tic
kml_data = kml_read_shapefile(kml_fn);
toc

lat = zeros(size(kml_data));
lon = zeros(size(kml_data));
for idx = 1:length(kml_data)
  lat(idx) = kml_data(idx).Y;
  lon(idx) = kml_data(idx).X;
end

%%%%%%%%%%%%%%%%%%%% Plot the KML data over a geotiff %%%%%%%%%%%%%%%%%%%%%

tic
proj = geotiffinfo(geotiff_fn);
[x,y] = projfwd(proj,lat,lon);
x = x/1e3;
y = y/1e3;
[RGB,R,bbox] = geotiffread(geotiff_fn);
R = R/1e3;
toc

figure(10); clf;
mapshow(RGB,R);
axis image off;
hold on;
mapshow(x,y,'DisplayType','point');
hold off;

return;


