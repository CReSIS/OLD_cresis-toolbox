% script opsScatterPlotExample
%
% Example for getting a set of point data and plotting in with Matlab
% scatter function.  This example is for ice bottom and quality levels.
%
% Author: John Paden
error('Copy and paste so that you do not accidently run the code again.')

%% Load data
% =========================================================================
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

% South East Greenland
% param.properties.bound = 'POLYGON((-37.63785548404194 65.791037592401,-39.47772995382079 66.40743234028808,-40.085584695963625 66.56251349919536,-40.92780583444692 65.95684010485714,-41.55260704096926 65.68255399572651,-42.00296733010212 65.67097831889556,-42.311344235465974 65.30760761366201,-42.73853797959285 64.45932901010737,-42.91333958085607 64.0724679771912,-43.55679770266661 63.65754904128014,-44.08337938609061 63.22611370587406,-44.42122603886692 62.810043497747145,-44.96382787195242 62.36364392520982,-44.973569343178276 61.96876658571534,-44.55240676821975 61.811906447649,-42.1097150701923 61.860645271528256,-41.760048142948605 61.86482228822533,-41.25349780678023 62.402399141442864,-37.50793635358621 65.75787662325365,-37.63785548404194 65.791037592401))';

% Northwest Greenland
% param.properties.bound = 'POLYGON((-71.0053362001043 76.16057398437214,-73.73259157503747 78.19959005672,-71.49950892720081 79.0874461988172,-68.42508855606353 79.59367567218798,-63.63151564199848 79.63285155615856,-58.22908776190653 78.57762161232762,-53.405119599821695 76.64563469686088,-51.5891094053563 75.55134105793663,-51.83985829375671 74.7060616370308,-50.994722479146866 73.96867672079739,-49.45684962554926 72.47232127429558,-48.39243992308465 71.46958353724628,-47.11705266779315 70.20890091915787,-47.06555852802419 69.97294629319876,-51.38645084070283 69.84971418122329,-52.51527759419749 70.21179896954024,-71.10657609316321 76.18943830966413,-71.0053362001043 76.16057398437214))';

fprintf('Getting points (%s)\n', datestr(now));
[status,data] = opsGetPointsWithinPolygon(sys,param);
fprintf('  Done (%s)\n', datestr(now));

%% Example to remove bad data
if 0
  keep_pnts = ~(strcmpi('2009_Greenland_TO_wise',data.properties.Season.') | isnan(data.properties.Bottom));
  data.properties.Lat = data.properties.Lat(keep_pnts);
  data.properties.Lon = data.properties.Lon(keep_pnts);
  data.properties.Bottom = data.properties.Bottom(keep_pnts);
  data.properties.Bottom_Quality = data.properties.Bottom_Quality(keep_pnts);
  data.properties.Elevation = data.properties.Elevation(keep_pnts);
end

%% Plot data
% =========================================================================

geotiff_fn = 'X:/GIS_data/greenland/Landsat-7/Greenland_natural_150m.tif';
proj = geotiffinfo(geotiff_fn);
[data.properties.X,data.properties.Y] = projfwd(proj,data.properties.Lat,data.properties.Lon);

% Scatter plot of bed elevations
figure(1); clf;
scatter(data.properties.X(1:10:end)/1e3, data.properties.Y(1:10:end)/1e3, [], data.properties.Elevation(1:10:end) - data.properties.Bottom(1:10:end));
axes_elev_only = gca;
h_colorbar = colorbar;
set(get(h_colorbar,'YLabel'), 'String', 'WGS-84 elevation (m)');
axis normal;
xlabel('X (km)');
ylabel('Y (km)');

% Scatter plot of bed elevations
h_fig = figure(2); clf;
[proj,fig_h] = plot_geotiff(geotiff_fn,[],[], h_fig);
hold on;
bed_elevation = data.properties.Elevation(1:10:end) - data.properties.Bottom(1:10:end);
scatter(data.properties.X(1:10:end)/1e3, data.properties.Y(1:10:end)/1e3, [], bed_elevation);
axes_elev = gca;
h_colorbar = colorbar;
caxis([min(bed_elevation(:)) max(bed_elevation(:))]);
set(get(h_colorbar,'YLabel'), 'String', 'WGS-84 elevation (m)');
axis normal;
xlabel('X (km)');
ylabel('Y (km)');

% Example for plotting quality levels
h_fig = figure(3); clf;
[proj,fig_h] = plot_geotiff(geotiff_fn,[],[], h_fig);
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
axes_quality = gca;
axis normal;
xlabel('X (km)');
ylabel('Y (km)');
legend('Good','Medium','Bad');

linkaxes([axes_elev axes_quality axes_elev_only],'xy');

return;

%% Example for getting a data point
figure(1);
[x,y] = ginput(1);
sub_sample_idxs = 1:10:length(data.properties.X);
dist = sqrt(abs(x - data.properties.X(sub_sample_idxs)/1e3).^2 + abs(y - data.properties.Y(sub_sample_idxs)/1e3).^2);
[min_dist,min_idx] = min(dist);
min_idx = sub_sample_idxs(min_idx);
fprintf('%s %s\n', data.properties.Season{min_idx}, data.properties.Frame{min_idx})

%% Example for saving the figures to a file
saveas(3,'H:\Kanger_bottom_colorbar.jpg');
saveas(1,'H:\Kanger_bottom.jpg');
saveas(2,'H:\Kanger_quality.jpg')
