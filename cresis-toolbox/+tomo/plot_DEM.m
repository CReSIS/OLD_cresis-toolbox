function [fig_h, ax, cax, data, crossover_frms, crossover_DEM, crossover_badmask] = plot_DEM(param, varargin)
% [fig_h, ax, cax, data, crossover_frms, ...
%  crossover_DEM, crossover_badmask] = plot_DEM(param, varargin)
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
% A 'title' argument may be passed in.
%  In the case of ('title', 'on') a title will be shown above the image.
%  Default = 'off';
%
% A 'crossover' argument may be passed in.
%  In the case of ('crossover', 'on') crossovers will be shown.
%  Default = 'off';
%
% Arguments for 'width' and 'height' (in pixels) may be passed in.
%
% Currently, only 2014 Greenland P3 DEM files will be accepted.
%
% Handles for the generated figure, axis, and color axis are returned,
%  respectively.
%
% For usage examples, see tomo.run_plot_DEM.m
%  and tomo.create_movie.m
%
% Authors: Victor Berger, John Paden
%
% See also: tomo.create_movie.m, tomo.run_plot_DEM.m

p = inputParser;
addRequired(p, 'param');
addOptional(p, 'param_override', '',  @(x)(isstruct(x)));
addOptional(p, 'gt', '', @(x)(isstruct(x)));
addOptional(p, 'title', 'off')
addOptional(p, 'width', 560, @isnumeric);
addOptional(p, 'height', 420, @isnumeric);
addOptional(p, 'crossover', 'off', @(x)(strcmp(x, 'on') || (strcmp(x, 'off'))));
addOptional(p, 'fill', 'on')
addOptional(p, 'cmaplim', '', @isnumeric);
parse(p, param, varargin{:});

physical_constants;

if exist('param_override', 'var')
  param = merge_structs(param, param_override);
end

if isempty(p.Results.gt)
  fprintf('\nLoading GeoTIFF, this may take a minute.. (%s)\n', datestr(now));
  gt.proj = geotiffinfo(param.geotiff_fn);
  [gt.RGB_bg, gt.R_bg, ~] = geotiffread(param.geotiff_fn);
  fprintf('  Done (%s)', datestr(now));
else
  gt = p.Results.gt;
end

if strcmp(p.Results.crossover, 'on')
  drawCO = 1;
else
  drawCO = 0;
end

fig_h          = figure;
ax             = gca;
fig_h.Units    = 'pixels';
fig_h.Position = [0 0 p.Results.width p.Results.height];

imagesc((gt.R_bg(3,1) + gt.R_bg(2,1)*(0:size(gt.RGB_bg,2)-1))/1e3, ...
  (gt.R_bg(3,2) + gt.R_bg(1,2)*(0:size(gt.RGB_bg,1)-1))/1e3, gt.RGB_bg,'Parent',ax);

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

set(ax,'YDir','normal');
hold on;
xlabel('X (km)');
ylabel('Y (km)');

if drawCO
  %% Create a buffer around the flowline
  max_angle = 500 / WGS84.semimajor;
  max_angle_deg = max_angle * 180/pi;
  
  [latbuffer, lonbuffer] = bufferm(lat, lon, 10*max_angle_deg);
  
  if any(isnan(latbuffer)) || any(isnan(lonbuffer))
    fprintf('\nAttention: frame %s has a self-intersection (this is usually not a problem)\n', sprintf('%s_%03d',param.day_seg, param.frm));
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
    [~,data] = opsGetPointsWithinPolygon(param.radar_name,OPSparams);
    fprintf('  Done (%s)', datestr(now));
  catch ME
    fprintf('\nFailed to load cross-overs for frame %s, so none will be shown.\n', sprintf('%s_%03d',param.day_seg, param.frm));
    drawCO = false;
  end
  
  if exist('data', 'var') && isfield(data.properties, 'Frame')
    crossover_frms = {};
    
    %% HACK -> restrict cross-over frames to 2014 Gr P3
    if 1
      keep_pnts = strcmpi('2014_Greenland_P3',data.properties.Season.');
      data.properties.Frame = data.properties.Frame(keep_pnts);
      crossover_frms = unique(data.properties.Frame);
    end
    
    crossover_DEM     = cell(1, length(crossover_frms));
    crossover_badmask = ones(1, length(crossover_frms));
    for crossover_idx = 1:length(crossover_frms)
      ds = param.day_seg;
      param.day_seg = crossover_frms{crossover_idx};
      param.day_seg = param.day_seg(1:end-4);
      
      try
        crossover_DEM{crossover_idx}.topDEM = load(fullfile(ct_filename_out(param, param.DEM_source, ''), sprintf('%s_top',crossover_frms{crossover_idx})));
        crossover_DEM{crossover_idx}.bottomDEM = load(fullfile(ct_filename_out(param, param.DEM_source, ''), sprintf('%s_bottom',crossover_frms{crossover_idx})));
        crossover_DEM{crossover_idx}.NaNDEM = crossover_DEM{crossover_idx}.bottomDEM.DEM;
        crossover_DEM{crossover_idx}.NaNDEM(crossover_DEM{crossover_idx}.NaNDEM == -32767) = NaN;
      catch ME
        fprintf('\nFailed to load cross-over frame %s', crossover_frms{crossover_idx});
        crossover_badmask(crossover_idx) = 0;
      end
      param.day_seg = ds;
    end
  end
  
  if drawCO && exist('data', 'var') && isfield(data.properties, 'Frame')
    for surf_idx = 1:length(crossover_frms)
      if ~crossover_badmask(surf_idx) || strcmp(crossover_frms{surf_idx},sprintf('%s_%03d',param.day_seg, param.frm))
        continue;
      end
      
      imagesc(crossover_DEM{surf_idx}.bottomDEM.xaxis/1e3,crossover_DEM{surf_idx}.bottomDEM.yaxis/1e3, ...
        crossover_DEM{surf_idx}.NaNDEM, 'Parent', ax, 'AlphaData', ~isnan(crossover_DEM{surf_idx}.NaNDEM));
      
      camlight
      shading(ax, 'interp');
    end
  end
end

param.NaNDEM = param.bottomDEM.DEM;
param.NaNDEM(param.NaNDEM == -32767) = NaN;
imagesc(param.bottomDEM.xaxis/1e3,param.bottomDEM.yaxis/1e3,param.NaNDEM, ...
  'Parent',ax, 'AlphaData', ~isnan(param.NaNDEM));

if isempty(p.Results.cmaplim)
  colormap(demcmap(param.bottomDEM.points.elev, 240));
else
  colormap(demcmap(p.Results.cmaplim, 240));
end
cb2 = colorbar;
cb2.Label.String = 'Bed elevation (WGS-84, m)';

plot(param.bottomDEM.points.x(round(size(param.bottomDEM.points.x, 1)/2), :)/1e3, ...
  param.bottomDEM.points.y(round(size(param.bottomDEM.points.y, 1)/2), :)/1e3, ...
  'r', 'LineWidth', 1 , 'Parent', ax);
legend(ax,'Flight line');
shading(ax, 'interp');

cax = caxis;

if strcmp(p.Results.title, 'on')
  title(sprintf('%s %s:  %s', param.radar_name, param.season_name, ...
    sprintf('%s_%03d',param.day_seg, param.frm)), ...
    'Interpreter','none','FontWeight','normal','FontSize', 15, 'Parent', ax);
end

axis(ax, axis_equal(ax, x, y, 'buffer', 2.5));

end