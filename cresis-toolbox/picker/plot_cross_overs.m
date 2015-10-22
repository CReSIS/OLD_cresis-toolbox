function plot_cross_overs(cross_over_fn,geotiff_fn,param,fig_h)
% plot_cross_overs(cross_over_fn,geotiff_fn,fig_h)
%
%  Takes cross over info from cross_overs_fn and plots on the geotiff.
%  Right now the code plots all cross overs. If you want to limit which
%  cross overs get plotted (e.g. non-zero thickness or angles of
%  intersection closer to 90 deg) then you have to modify the code.  An
%  example of this is included in the code.
%
%  cross_overs_fn = file created by cross_over_analysis.m
%  geotiff_fn = location of geotiff file to plot cross overs on top of
%  param (optional) = params for flightline read in. See
%  run_cross_over_analysis and cross_over_analysis for details.
%  fig_h = figure handle to use when plotting cross overs
%    default is to create a new figure
%
% Example:
%  cross_over_fn = '2009_cross_over.mat';
%  geotiff_fn = '/cresis/scratch1/mdce/csarp_support/Landsat-7/antarctica/Antarctica_LIMA_peninsula.tif';
%  fig_h = 11;
%  plot_cross_overs(cross_over_fn,geotiff_fn,fig_h);
%  figure(fig_h); caxis([0 150]);
%
%  cross_over_fn = '2011_cross_over.mat';
%  geotiff_fn = '/cresis/data3/GIS_data/greenland/Landsat-7/Greenland_natural_90.tif';
%  fig_h = 11;
%  plot_cross_overs(cross_over_fn,geotiff_fn,fig_h);
%  figure(fig_h); caxis([0 150]);
%
% Author: Brady Maasen, John Paden
%
% See also plot_cross_overs, run_cross_over_analysis

tic;

fprintf('Loading the data and geotiff (%.1f sec)\n', toc);
X = load(cross_over_fn);
proj = geotiffinfo(geotiff_fn);
if exist('param','var')
  if param.data_type.type == 1 || param.data_type.type == 3 || param.data_type.type == 4
    fid = fopen(param.fin,'r');
    XFL = textscan(fid,param.data_type.format,'delimiter',param.data_type.delim, ...
      'headerlines',param.data_type.headerlines);
    lat_fl = XFL{param.data_type.col_lat};
    lon_fl = XFL{param.data_type.col_lon};
  end
end

if 0
  % Restrict which cross overs are considered.
  good_idxs = find(abs(X.cross_angle-90) < 60 & X.Thickness(:,1) > 0);
  X.Latitude = X.Latitude(good_idxs,:);
  X.Longitude = X.Longitude(good_idxs,:);
  X.Thickness = X.Thickness(good_idxs,:);
end

% Read the image
[RGB, R, tmp] = geotiffread(geotiff_fn);
R = R/1e3;

if ~exist('fig_h','var') || isempty(fig_h)
  figure; clf;
else
  figure(fig_h); clf;
end
mapshow(RGB, R);
hold on;
if exist('lat_fl','var');
  fprintf('Plotting flightlines (%.1f sec)\n', toc);
  lat_fl = downsample(XFL{param.data_type.col_lat},param.plot_samp);
  lon_fl = downsample(XFL{param.data_type.col_lon},param.plot_samp);
  [x_fl,y_fl] = projfwd(proj,lat_fl,lon_fl);
  x_fl = x_fl/1e3;
  y_fl = y_fl/1e3;
  scatter(x_fl,y_fl,10,'k.')
  hold on;
end
fprintf('Plotting cross overs (%.1f sec)\n', toc);
[x,y] = projfwd(proj,X.Latitude(:,1),X.Longitude(:,1));
x = x/1e3;
y = y/1e3;

thickness = abs(X.Thickness(:,1) - X.Thickness(:,2));
[thickness,sorted_idx] = sort(thickness,1,'descend');
x = x(sorted_idx);
y = y(sorted_idx);
scatter(x,y,500,thickness,'.');
caxis([0 max(thickness)]);
hold off;
xlabel('X (km)');
ylabel('Y (km)');

axis([min(x)-50 max(x)+50 min(y)-50 max(y)+50]);

h = colorbar;
set(get(h,'YLabel'),'String','absolute error (m)');

rms = sqrt(mean(thickness.^2));
rMs = sqrt(median(thickness.^2));
fprintf('\nThe RMS value for thickness error is %f\n',rms);
fprintf('The Root Median Square value for thickness error is %f\n\n', rMs);


end