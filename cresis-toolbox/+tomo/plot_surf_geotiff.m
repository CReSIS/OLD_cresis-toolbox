function [fig_h, ax, cax] = plot_surf_geotiff(param, param_override, gt, prop)
% [fig_h, ax, cax] = tomo.plot_surf_geotiff(param, param_override, gt, prop)
%
% Generate a DEM image containing the reference frame bottom, 
%  the frames that intersect (cross-over) the reference,
%  and a map of the surrounding area from the GeoTIFF.
%
% The param struct is expected to contain the fields
%  radar_name, season_name, DEM_source, location,
%  day_seg, frm, and geotiff_fn.
%
% A GeoTIFF object "gt" can be passed in. Format:
%   gt.proj = geotiffinfo(param.geotiff_fn);
%   [gt.RGB_bg, gt.R_bg, ~] = geotiffread(param.geotiff_fn);
%
% A 'notitle' argument may be passed in.
%  In this case, no title will be shown above the image.
%
% Currently, only 2014 Greeland P3 DEM files will be accepted.
%
% Handles for the generated figure, axis, and color axis are returned,
%  respectively.
%
% For usage examples, see tomo.run_plot_surf_geotiff.m 
%  and tomo.create_movie.m
%
% Authors: Victor Berger, John Paden
%
% See also: tomo.create_movie.m, tomo.run_plot_surf_geotiff.m

physical_constants;
param = merge_structs(param, param_override);

if ~exist('gt','var')
  fprintf('\nLoading GeoTIFF, this may take a minute.. (%s)\n', datestr(now));
  gt.proj = geotiffinfo(param.geotiff_fn);
  [gt.RGB_bg, gt.R_bg, ~] = geotiffread(param.geotiff_fn);
  fprintf('  Done (%s)', datestr(now));
end

fig_h = figure;
ax = gca;

imagesc((gt.R_bg(3,1) + gt.R_bg(2,1)*(0:size(gt.RGB_bg,2)-1))/1e3, ...
  (gt.R_bg(3,2) + gt.R_bg(1,2)*(0:size(gt.RGB_bg,1)-1))/1e3, gt.RGB_bg,'Parent',ax);

set(ax,'YDir','normal');
axis(ax, 'equal')
hold on;
xlabel('X (km)');
ylabel('Y (km)');

if ~isfield(param, 'topDEM') || ~isfield(param, 'bottomDEM')
  try
    param.topDEM    = load(fullfile(ct_filename_out(param, param.DEM_source, ''), sprintf('%s_%03d_top',param.day_seg, param.frm)));
    param.bottomDEM = load(fullfile(ct_filename_out(param, param.DEM_source, ''), sprintf('%s_%03d_bottom',param.day_seg, param.frm)));
  catch ME
    warning('Attention: DEM for reference frame %s not found. Skipping this frame.', sprintf('%s_%03d',param.day_seg, param.frm));
    keyboard
  end
end

lat = param.bottomDEM.points.lat(32, :);
lon = param.bottomDEM.points.lon(32, :);
lat = lat(1:10:end);
lon = lon(1:10:end);

[x,y] = projfwd(gt.proj,lat,lon);

if max(x)/1e3-min(x)/1e3 > max(y)/1e3-min(y)/1e3
  dif = max(x)/1e3-min(x)/1e3;
  axis([min(x)/1e3-2, max(x)/1e3+2, (min(y)/1e3-2)-dif/2, (max(y)/1e3+2)+dif/2]);
else
  dif = max(y)/1e3-min(y)/1e3;
  axis([(min(x)/1e3-2)-dif/2, (max(x)/1e3+2)+dif/2, min(y)/1e3-2, max(y)/1e3+2]);
end

%% Create a buffer around the flowline
max_angle = 500 / WGS84.semimajor;
max_angle_deg = max_angle * 180/pi;

[latbuffer, lonbuffer] = bufferm(lat, lon, max_angle_deg);

if any(isnan(latbuffer)) || any(isnan(lonbuffer))
  fprintf('\nAttention: frame %s has a self-intersection (this is usually not a problem)\n', sprintf('%s_%03d',param.day_seg, frm_idx));
  latbuffer = latbuffer(find(isnan(latbuffer), 1, 'last')+1:end);
  lonbuffer = lonbuffer(find(isnan(lonbuffer), 1, 'last')+1:end);
end

%% Convert bufferm to WKT polygon
WKT_polygon = [sprintf('POLYGON(('), ...
  sprintf('%0.14g %0.14g,',[lonbuffer(1:end-1).'; latbuffer(1:end-1).']), ...
  sprintf('%0.14g %0.14g',[lonbuffer(end).'; latbuffer(end).']), sprintf('))')];

%% Get radar points from OPS
OPSparams = [];
OPSparams.properties.location = param.location;
OPSparams.properties.bound = WKT_polygon;
fprintf('\nLoading radar points from OPS, this may take a minute.. (%s)\n', datestr(now));

try
  [~,data] = opsGetFramesWithinPolygon(param.radar_name,OPSparams);
  fprintf('  Done (%s)', datestr(now));
catch ME
  fprintf('\nFailed to load cross-overs for frame %s, so none will be shown.\n', sprintf('%s_%03d',param.day_seg, param.frm));
  drawCO = false;
end

if exist('data', 'var') && isfield(data, 'frame')
  %% HACK -> restrict cross-over frames to 2014 Gr P3
  crossover_frms = {};
  for data_idx = 1:length(data.frame)
    if ~isempty(strfind(data.frame{data_idx}, '2014'))
      crossover_frms{end+1} = data.frame(data_idx);
    end
  end
  %%
  
  crossover_DEM     = cell(1, length(crossover_frms));
  crossover_badmask = ones(1, length(crossover_frms));
  for crossover_idx = 1:length(crossover_frms)
    ds = param.day_seg;
    param.day_seg = crossover_frms{crossover_idx}{1};
    param.day_seg = param.day_seg(1:end-4);
    
    try
      crossover_DEM{crossover_idx}.topDEM = load(fullfile(ct_filename_out(param, param.DEM_source, ''), sprintf('%s_top',crossover_frms{crossover_idx}{1})));
      crossover_DEM{crossover_idx}.bottomDEM = load(fullfile(ct_filename_out(param, param.DEM_source, ''), sprintf('%s_bottom',crossover_frms{crossover_idx}{1})));
      crossover_DEM{crossover_idx}.NaNDEM = crossover_DEM{crossover_idx}.bottomDEM.DEM;
      crossover_DEM{crossover_idx}.NaNDEM(crossover_DEM{crossover_idx}.NaNDEM == -32767) = NaN;
    catch ME
      fprintf('\nFailed to load cross-over frame %s', crossover_frms{crossover_idx}{1});
      crossover_badmask(crossover_idx) = 0;
    end
    param.day_seg = ds;
  end
end

if exist('data', 'var') && isfield(data, 'frame')
  for surf_idx = 1:length(crossover_frms)
    if ~crossover_badmask(surf_idx) || strcmp(crossover_frms{surf_idx}{1},sprintf('%s_%03d',param.day_seg, param.frm))
      continue;
    end
    
    imagesc(crossover_DEM{surf_idx}.bottomDEM.xaxis/1e3,crossover_DEM{surf_idx}.bottomDEM.yaxis/1e3, ...
      crossover_DEM{surf_idx}.NaNDEM, 'Parent', ax, 'AlphaData', ~isnan(crossover_DEM{surf_idx}.NaNDEM));
    
    camlight
    shading(ax, 'interp');
  end
end

param.NaNDEM = param.bottomDEM.DEM;
param.NaNDEM(param.NaNDEM == -32767) = NaN;
imagesc(param.bottomDEM.xaxis/1e3,param.bottomDEM.yaxis/1e3,param.NaNDEM, ...
  'Parent',ax, 'AlphaData', ~isnan(param.NaNDEM));

colormap(demcmap(param.bottomDEM.points.elev, 24));
cb2 = colorbar;
cb2.Label.String = 'Bed elevation (WGS-84, m)';
plot(param.bottomDEM.points.x(32, :)/1e3, param.bottomDEM.points.y(32, :)/1e3, 'k', 'LineWidth', 1 , 'Parent', ax);
legend(ax,'Flight line');
shading(ax, 'interp');

cax = caxis;

if ~exist('prop', 'var') || ~strcmp(prop, 'notitle')
  title(sprintf('%s %s:  %s', param.radar_name, param.season_name, ...
    sprintf('%s_%03d',param.day_seg, param.frm)), ...
    'Interpreter','none','FontWeight','normal','FontSize', 15, 'Parent', ax);
end

min_x = min( inf, min(param.bottomDEM.points.x(:)));
max_x = max(-inf, max(param.bottomDEM.points.x(:)));
min_y = min( inf, min(param.bottomDEM.points.y(:)));
max_y = max(-inf, max(param.bottomDEM.points.y(:)));

figure_dots_per_km = 20;
map_buffer = 2e3;

map_min_x = min_x-map_buffer;
map_max_x = max_x+map_buffer;
map_min_y = min_y-map_buffer;
map_max_y = max_y+map_buffer;

axis(ax, [map_min_x map_max_x map_min_y map_max_y]/1e3);
set(ax,'Units','pixels');
map_axes = get(ax,'Position');
set(fig_h,'Units','pixels');
map_pos = get(fig_h,'Position');

map_new_axes = map_axes;
map_new_axes(3) = round(figure_dots_per_km*(map_max_x-map_min_x)/1e3);
map_new_axes(4) = round(figure_dots_per_km*(map_max_y-map_min_y)/1e3);
map_pos(3) = map_new_axes(3) + map_pos(3)-map_axes(3);
map_pos(4) = map_new_axes(4) + map_pos(4)-map_axes(4);

set(fig_h,'Position',map_pos);
set(ax,'Position',map_new_axes);
set(fig_h,'PaperPositionMode','auto');

end