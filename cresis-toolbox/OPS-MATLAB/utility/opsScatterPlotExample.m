% script opsScatterPlotExample
%
% Example for getting a set of point data and plotting in with Matlab
% scatter function.  This example is for ice bottom and quality levels.
%
% Author: John Paden
error('Copy and paste so that you do not accidently run the code again.')

sys = 'rds';
param = [];
param.properties.location = 'arctic';

% Polygon's at ops3.cresis.ku.edu --> Arctic --> Draw polygon (double click
% to close polygon)

% Helheim polygon
% param.properties.bound = 'POLYGON((-39.5157619560387 66.60645677740496,-39.38100914546357 66.75651423130941,-39.276083937124795 66.91554160493561,-38.75086344123707 66.92169343817854,-38.07075943960155 66.89418468857887,-37.57168206589797 66.79263788284699,-37.478105950021 66.62523970060901,-37.45434464204945 66.41924919176293,-37.36530763195985 66.24544098233096,-37.72331160659728 66.18545358686922,-38.48213912070531 66.19147777690348,-39.15750707699793 66.22101527577271,-39.68089836023789 66.33020636648097,-39.5157619560387 66.60645677740496))';

% Kangerdlugssuaq polygon
% param.properties.bound = 'POLYGON((-32.89308751162966 68.87842019516255,-33.330159367349694 68.87979361675448,-33.57748560590095 68.75852267191158,-33.70663957228397 68.67161240038641,-33.63572567131842 68.55200565799642,-33.31876919364127 68.50453128906366,-32.909363685413744 68.44792509908714,-32.69201822568974 68.46018627042987,-32.53167558704438 68.5416290077124,-32.45332111243092 68.60922340931619,-32.42609559511996 68.68440834748068,-32.8961014755777 68.87865333956601,-32.89308751162966 68.87842019516255))';

% Jakobshavn polygon
% param.properties.bound = 'POLYGON((-50.32303131774781 70.02913831585758,-49.93715335904083 68.45463235949325,-41.73222625446042 68.48970611501662,-41.507788647681096 70.13545713679203,-50.32303131774781 70.02913831585758))';

% Jakobshavn polygon (channel only)
param.properties.bound = 'POLYGON((-49.87839454285566 69.19467094239015,-49.775002993686556 69.06887962878123,-49.457730724927146 69.0166739683372,-49.067262963224636 69.0248718726718,-48.721804586124016 69.11021474395748,-48.486486811488994 69.15897709901661,-48.35032188399231 69.16856094682609,-48.42737255389091 69.27285397608863,-48.97718851464709 69.20223457154142,-49.287828343319354 69.18644546855747,-49.50560454288413 69.20616427532576,-49.7115935470338 69.21366353502647,-49.87839454285566 69.19467094239015))';

fprintf('Getting points (%s)\n', datestr(now));
[status,data] = opsGetPointsWithinPolygon(sys,param);
fprintf('  Done (%s)\n', datestr(now));

tmp = data;% Example for plotting bottom elevationkeep_pnts = ~(strcmpi('2009_Greenland_TO_wise',data.properties.Season.') | isnan(data.properties.Bottom));data.properties.Lat = data.properties.Lat(keep_pnts);data.properties.Lon = data.properties.Lon(keep_pnts);data.properties.Bottom = data.properties.Bottom(keep_pnts);data.properties.Bottom_Quality = data.properties.Bottom_Quality(keep_pnts);data.properties.Elevation = data.properties.Elevation(keep_pnts);return

[data.properties.X,data.properties.Y] = projfwd(proj,data.properties.Lat,data.properties.Lon);figure(3);scatter(data.properties.X(1:10:end)/1e3, data.properties.Y(1:10:end)/1e3, [], data.properties.Elevation(1:10:end) - data.properties.Bottom(1:10:end));axes_elev_only = gca;colorbar;axis normal;figure(1); clf;
geotiff_fn = 'X:/GIS_data/greenland/Landsat-7/Greenland_natural_150m.tif';
[proj,fig_h] = plot_geotiff(geotiff_fn,[],[], 1);
hold on;
scatter(data.properties.X(1:10:end)/1e3, data.properties.Y(1:10:end)/1e3, [], data.properties.Elevation(1:10:end) - data.properties.Bottom(1:10:end));
hold off;
colorbar;
axes_elev = gca;
axis normal;
% Example for plotting quality levels
figure(2); clf;
geotiff_fn = 'X:/GIS_data/greenland/Landsat-7/Greenland_natural_150m.tif';
[proj,fig_h] = plot_geotiff(geotiff_fn,[],[], 2);
hold on;
quality_color = {'g.' 'y.' 'r.'};
for quality_level = [1 2 3]
  X = data.properties.X(1:10:end)/1e3;
  Y = data.properties.Y(1:10:end)/1e3;
  Q = data.properties.Bottom_Quality(1:10:end);
  quality_mask = Q == quality_level;
  if any(quality_mask)
    h(quality_level) = plot(X(quality_mask), Y(quality_mask), quality_color{quality_level});
  end
end
hold off;
axes_quality = gca;
axis normal;
linkaxes([axes_elev axes_quality axes_elev_only],'xy');

return;

% Example for getting a data point
figure(1);
[x,y] = ginput(1);
sub_sample_idxs = 1:10:length(data.properties.X);
dist = sqrt(abs(x - data.properties.X(sub_sample_idxs)/1e3).^2 + abs(y - data.properties.Y(sub_sample_idxs)/1e3).^2);
[min_dist,min_idx] = min(dist);
min_idx = sub_sample_idxs(min_idx);
fprintf('%s %s\n', data.properties.Season{min_idx}, data.properties.Frame{min_idx})



saveas(3,'H:\Kanger_bottom_colorbar.jpg')saveas(1,'H:\Kanger_bottom.jpg')saveas(2,'H:\Kanger_quality.jpg')
