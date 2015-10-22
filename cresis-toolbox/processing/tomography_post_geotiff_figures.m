% script tomography_post_geotiff_figures
%
% Takes the .mat file output from tomography_DEM.m and creates
% a geotiff and a 3d image
%
% Author: John Paden

% ================================================================
%% User Settings
% ================================================================

% in_fn: string containing the output from tomography_DEM.m that you want
% to form a geotiff of or 3d image of
in_fn = 'C:\Users\radar\Desktop\NEEM_NGRIP_3d.mat';

% geotiff_en: Flag to create the geotiff file of the DEM
geotiff_en = true;

% image3d_en: Flag to create a fancy picture of the DEM
image3d_en = true;

% ProjectedCSTypeGeoKey: The projection type key to use with geotiff_en.
% The most reliable way to
% set this is to create a geotiff file with the right projection in
% another program and read in the key from that program using geotiffinfo.
% proj = geotiffinfo('/cresis/projects/GIS_data/greenland/Landsat-7/Greenland_natural_150m.tif');
% proj.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey
%  ans =
%       32767
% Some keys are not fully supported and you'll just have to try different
% projections until you find one that works.  Test the output file by
% loading it with Matlab and another program like Arc or Globalmapper
ProjectedCSTypeGeoKey = 32622;
ProjectedCSTypeGeoKey = 3031;

% geotiff_fn: string containing the output geotiff filename
geotiff_fn = 'C:\Users\radar\Desktop\NEEM_NGRIP_3d.tif';

% ================================================================
%% Automated Section
% ================================================================

load(in_fn);

DEM = DEM(2:2:end-1,2:2:end-1);
northAxis = northAxis(2:2:end-1);
eastAxis = eastAxis(2:2:end-1);

if geotiff_en
  %% Create the geotiff file of the DEM
  R = spatialref.MapRasterReference;
  R.XLimWorld = [eastAxis(1) eastAxis(end)];
  R.YLimWorld = [northAxis(1) northAxis(end)];
  R.RasterSize = size(DEM);
  R.RasterInterpretation = 'cells';
  R.ColumnsStartFrom = 'south';
  R.RowsStartFrom = 'west';
  R.RasterInterpretation = 'postings';
  
  key.GTModelTypeGeoKey  = 1;  % Projected Coordinate System (PCS)
  key.GTRasterTypeGeoKey = 2;  % PixelIsPoint
  key.ProjectedCSTypeGeoKey = ProjectedCSTypeGeoKey;
  DEM(isnan(DEM)) = -9999;
  
  fprintf('Creating output %s\n', geotiff_fn);
  geotiffwrite(geotiff_fn, int16(DEM), R, 'GeoKeyDirectoryTag', key);

end

if image3d_en
  %% Create a 3-D surface image (fancy looking figure)
  fig_h = 20;
  figure(fig_h); clf;
  hA2 = axes;
  hC = surf((eastAxis-eastAxis(1))/1e3,(northAxis-northAxis(1))/1e3,double(medfilt2(DEM,[3 11])));
  set(hC,'EdgeAlpha',0.2); grid off;
  set(hA2,'View',[190 80]);
  hold on;
  % h = plot3((neem.east-eastAxis(1))/1e3,(neem.north-northAxis(1))/1e3,neem.DEM,'ko');
  % set(h, 'LineWidth', 3);
  % set(h, 'MarkerSize', 7);
  hold off;
%   hx = xlabel('Easting (km)');
%   set(hx,'Position',[30 15 0]);
%   set(hx,'Rotation',-3);
%   hy = ylabel('Northing (km)');
%   set(hy,'Position',[-3 9 -1500]);
%   set(hy,'Rotation',57);
  zlabel(sprintf('WGS-84 bed\nelevation (m)'));
  grid on;
  hC = colorbar;
  set(get(hC,'YLabel'),'String','WGS-84 bed elevation (m)');
  set(hA2,'Position',[0.12 0.11 0.72 0.815])
  set(hC,'Position',[0.9 0.11 0.022 0.815])
  
  set(fig_h,'PaperOrientation','landscape');
  set(fig_h,'PaperPosition',[0.5 0.5 10 3]);
  
%   print(sprintf('-f%d',fig_h),'-djpeg','-r200','C:\Users\radar\Desktop\NEEM_NGRIP_3D_all.jpg');
  
end

return;
