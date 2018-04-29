function map_info = publish_map(cmd,param,map_info)
% map_info = publish_map(cmd,param,map_info)
%
% Function for plotting maps to go with echograms.
%
% cmd = string
%  'setup': loads geotiffs and sets of figure(s)
%  'plot': plots the particular frame onto the figure
%  'delete': deletes all temporary plots from plot command
%
% param = structure controlling how setup and plot is done
%  .location: string ('Greenland', 'Antarctica', 'Canada', 'Norway', or 'Arctic')
%  .type: string ('combined','singles','contour','contour-singles')
%    combined = combines all three map plots into one figure
%    singles = plot each map into a separate figure
%    contour = plots contour maps only (useful for north pole)
%    contour-singles = same as contour, but in separate figures
%  .fig_hand: base figure handle to use (handles will be fig_hand,
%    fig_hand+1, etc
%  .resample = resample the images (important when creating ps, eps or
%    pdf outputs of the figures because the full resolution will
%    be stored to the file if you don't resample)
%  .decimate_seg = decimate the trajectory to reduce the file size
%  .year/.day/.month = necessary to retreive correct ASCAT geoTIFF for sea
%   ice maps. 
%
% Uses some global variables to reduce memory consumption
%
% Example:
%
%  param.location = 'Greenland';
%  param.type = 'singles';
%  param.fig_hand = 11;
%  map_info = publish_map('setup',param);
%  fn = 'Z:\mdce\icards\1993_Greenland_P3\CSARP_standard\19930624_01\Data_19930624_01_014.mat';
%  load(fn);
%  [X,Y] = projfwd(map_info.proj,Latitude,Longitude);
%  param.day_seg_x = X;
%  param.day_seg_y = Y;
%  param.frame_X = X;
%  param.frame_Y = Y;
%  param.map_title = '19930624_01_014';
%  param.decimate_seg = true;
%  map_info = publish_map('plot',param,map_info);
%  saveas(11,'~/tmp/test_region.fig');
%  saveas(12,'~/tmp/test_mag.fig');
%  saveas(13,'~/tmp/test_contour.fig');
%
% Author: John Paden

% Store large map variables in global variables to help reduce
% memory usage
global h_region_CData;
global h_region_XData;
global h_region_YData;
global h_zoom_CData;
global h_zoom_XData;
global h_zoom_YData;

if strcmpi(cmd,'setup')
  % ======================================================================
  % ======================================================================
  % SETUP COMMAND
  % ======================================================================
  % ======================================================================
  
  if strcmpi(param.location,'Greenland')
    hemisphere = 'north';
    geotiff_fn = fullfile(ct_filename_gis(param,'greenland'),'Landsat-7','mzl7geo_90m_lzw.tif');
    cmap = gray(256);
    map_axis = [-1000 1000 -3500 -500];
    contour_position = [0.12 0.72 0.15 0.25];
    map_sub_title = 'Polar Stereograph 70N/-45E';
    
  elseif strcmpi(param.location,'Antarctica')
    if strcmpi(param.type,'ASCAT')
      geotiff_fn = fullfile(ct_filename_gis(param,[]),'antarctica','ASCAT',sprintf('ASCAT_%04d%02d%02d.tif', param.year, param.month, param.day));
      map_axis = [-3000 3500 -3000 4000];
      hemisphere = 'south';
      map_sub_title = 'Polar Stereograph -71/0E';
      map_sub_title2 = sprintf('ASCAT \\sigma_0 at 40\\circ incidence for %02d/%02d/%04d', param.month, param.day, param.year);
    else
      hemisphere = 'south';
      coastline = shaperead('landareas', 'UseGeoCoords', true,...
        'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});
      geotiff_fn = fullfile(ct_filename_gis(param,'antarctica'),'Landsat-7','Antarctica_LIMA.tif');
      cmap = gray(256);
      map_axis = [-3000 3000 -2500 2500];
      contour_position = [0.12 0.72 0.3 0.25];
      map_sub_title = 'Polar Stereograph -71N/0E';
    end
    
  elseif strcmpi(param.location,'Canada')
    hemisphere = 'north';
    geotiff_fn = fullfile(ct_filename_gis(param,'canada'),'Landsat-7','Canada_150m.tif');
    cmap = gray(256);
    map_axis = [-3000 1000 -3499 -500];
    contour_position = [0.12 0.72 0.3 0.25];
    map_sub_title = 'Polar Stereograph 70N/-45E';
    
  elseif strcmpi(param.location,'Norway')
    hemisphere = 'north';
    geotiff_fn = fullfile(ct_filename_gis(param,'norway'),'norway.tif');
    cmap = gray(256);
    map_axis = [544 762 8794 8976];
    contour_position = [0.12 0.72 0.3 0.25];
    map_sub_title = '';
     
  elseif strcmpi(param.location,'Arctic')
     if strcmpi(param.type,'ASCAT')
      geotiff_fn = fullfile(ct_filename_gis(param,[]),'arctic','ASCAT',sprintf('ASCAT_%04d%02d%02d.tif', param.year, param.month, param.day));
      map_axis = [-3500 3500 -5000 5000];
      hemisphere = 'north';
      map_sub_title = 'Polar Stereograph 70N/-45E';
      map_sub_title2 = sprintf('ASCAT \\sigma_0 at 40\\circ incidence for %02d/%02d/%04d', param.month, param.day, param.year);
     else
      hemisphere = 'north';
      geotiff_fn = fullfile(ct_filename_gis(param,'arctic'),'NaturalEarth_Data','Arctic_NaturalEarth.tif');
      cmap = jet(256);
      cmap(255,:) = [0.5 0.5 0.5];
      cmap(256,:) = [0.5 0.5 0.5];
      map_axis = [-3000 1400 -2000 2000];
      contour_position = [0.1200    0.7100    0.2700    0.2700];
      map_sub_title = 'Polar Stereograph 70N/-45E';
     end
    
  else
    % If this is not one of the standard locations, this must be a geotiff
    % filename and the type must be geotiff... if not, it is an error.
    if ~strcmpi(param.type,'geotiff')
      error('Unsupported location %s\n', param.location);
    end
  end

  if strcmpi(param.type,'geotiff')
    geotiff_fn = ct_filename_gis(param,param.geotiff{1});
    hemisphere = [];
    map_sub_title = param.geotiff{2};
  end

  if strcmp(hemisphere,'north')
    coastlines = shaperead('landareas','UseGeoCoords', true);
    coastline.Lat = [];
    coastline.Lon = [];
    for shape_idx = 1:length(coastlines)
      if any(coastlines(shape_idx).Lat > 60)
        coastline.Lat = cat(2,coastline.Lat,coastlines(shape_idx).Lat);
        coastline.Lon = cat(2,coastline.Lon,coastlines(shape_idx).Lon);
      end
    end
  end

  geotiff_fn_idx = find(param.type == ',',1) + 1;
  if ~isempty(geotiff_fn_idx)
    geotiff_fn = fullfile(ct_filename_gis(param,''),param.type(geotiff_fn_idx:end));
    param.type = param.type(1:geotiff_fn_idx-2);
  end

  if strcmpi(param.type, 'ASCAT')
    if length(param.fig_hand)>=1 && ishandle(param.fig_hand(1))
      clf(param.fig_hand(1));
    else
      if length(param.fig_hand)>=1 && isnumeric(param.fig_hand(1))
        param.fig_hand(1) = figure(param.fig_hand(1));
      else
        param.fig_hand(1) = figure();
      end
    end
    map_info.fig_hand(1) = param.fig_hand(1);
    set(param.fig_hand(1),'Renderer','painters');
    if ~exist(geotiff_fn,'file')
      error('%s not found, run ascat_to_geotiff for this segment.', geotiff_fn);
    end
    % Get the projection information
    proj = geotiffinfo(geotiff_fn);
    if strcmpi(param.location,'Arctic')
      warning('Applying fix to bad geotiff created by ascat_to_geotiff');
      proj.ProjParm(1) = 70;
    end
    % Read the image
    [RGB, R, tmp] = geotiffread(geotiff_fn);
    % =======================================================================
    % Large raster image map
    ah_region = axes('parent',param.fig_hand(1));
    %set(ah_region,'Position',[0.50000 0.1000 0.8*2/3 0.8000]);
    set(ah_region, 'Position',[0.54 0.15 0.8*2/3 0.75]);
    R = R/1e3;
    h_region = mapshow(RGB, R,'Parent',ah_region);
    h_region_CData = get(h_region,'CData');
    h_region_XData = get(h_region,'XData');
    h_region_YData = get(h_region,'YData');
    h_region_YData = linspace(h_region_YData(1),h_region_YData(end),size(h_region_CData,1));
    ah_region_orig_axis = axis(ah_region);
    xlabel(ah_region, 'X (km)')
    ylabel(ah_region, 'Y (km)')
    % =======================================================================
    % Zoomed out map...
    ah_map = axes('parent',param.fig_hand(1));
    h_region = mapshow(RGB, R,'Parent',ah_map);    
    axis(ah_map, map_axis)
    %set(ah_map,'Position',[0.01 0.1 0.5 .8])
    if strcmpi(param.location, 'Antarctica')
      set(ah_map, 'Position',[0.05 0.18 0.49 .6]);
      set(param.fig_hand(1),'PaperPosition',[0.5 0.5 8 6]);
    else
      set(ah_map, 'Position',[0.05 0.11 0.49 .78]);
      set(param.fig_hand(1),'PaperPosition',[0.5 0.5 8 6]);
    end
    xlabel(ah_map,'X (km)')
    ylabel(ah_map,'Y (km)')
    set(param.fig_hand(1),'PaperOrientation','Portrait');
    title(map_sub_title2);
    h_region_CData = get(h_region,'CData');
    h_region_XData = get(h_region,'XData');
    h_region_YData = get(h_region,'YData');
    h_region_YData = linspace(h_region_YData(1),h_region_YData(end),size(h_region_CData,1));
    ah_zoom = ah_region;
    h_zoom = [];
    
  elseif strcmpi(param.type, 'geotiff')
    if length(param.fig_hand)>=1 && ishandle(param.fig_hand(1))
      clf(param.fig_hand(1));
    else
      if length(param.fig_hand)>=1 && isnumeric(param.fig_hand(1))
        param.fig_hand(1) = figure(param.fig_hand(1));
      else
        param.fig_hand(1) = figure();
      end
    end
    map_info.fig_hand(1) = param.fig_hand(1);
    set(param.fig_hand(1),'Renderer','painters');
    % Get the projection information
    proj = geotiffinfo(geotiff_fn);
    % Read the image
    [RGB, R, tmp] = geotiffread(geotiff_fn);
    if size(RGB,3) == 3 && strcmp(class(RGB),'uint16') && max(RGB(:)) <= 255
        RGB = uint8(RGB);
    end
    if strcmpi(class(RGB),'int16') || strcmpi(class(RGB),'single')
        RGB = double(RGB);
    end
    % =======================================================================
    % Raster image map
    ah_region = axes('parent',param.fig_hand(1));
    %set(ah_region,'Position',[0.50000 0.1000 0.8*2/3 0.8000]);
    set(ah_region, 'Position',[0.1 0.1 0.8 0.8]);
    R = R/1e3;
    h_region = mapshow(RGB, R,'Parent',ah_region);
    h_region_CData = get(h_region,'CData');
    h_region_XData = get(h_region,'XData');
    h_region_YData = get(h_region,'YData');
    h_region_YData = linspace(h_region_YData(1),h_region_YData(end),size(h_region_CData,1));
    ah_region_orig_axis = axis(ah_region);
    xlabel(ah_region, 'X (km)')
    ylabel(ah_region, 'Y (km)')

    set(param.fig_hand(1),'PaperOrientation','Portrait');
    ah_map = [];
    ah_zoom = [];
    h_zoom = [];

 elseif strcmpi(param.type,'combined')
    if length(param.fig_hand)>=1 && ishandle(param.fig_hand(1))
      clf(param.fig_hand(1));
    else
      if length(param.fig_hand)>=1 && isnumeric(param.fig_hand(1))
        param.fig_hand(1) = figure(param.fig_hand(1));
      else
        param.fig_hand(1) = figure();
      end
    end
    map_info.fig_hand(1) = param.fig_hand(1);
    set(param.fig_hand(1),'Renderer','painters');
    % =======================================================================
    % Large raster image map
    ah_region = axes('parent',param.fig_hand(1));
    set(ah_region,'Position',[0.45000 0.1000 0.8*2/3 0.8000]);
    % Get the projection information
    proj = geotiffinfo(geotiff_fn);
    % Read the image
    [RGB, R, tmp] = geotiffread(geotiff_fn);
    R = R/1e3;
    % Plot the image
    if size(RGB,3) == 1
      % Special code for Arctic sea ice maps
      RGB(RGB == 251) = 254;
      RGB(RGB == 255) = 255;
      colormap(cmap);
      h_region = mapshow(RGB, cmap, R,'Parent',ah_region);
    else
      h_region = mapshow(RGB, R,'Parent',ah_region);
    end
    h_region_CData = get(h_region,'CData');
    h_region_XData = get(h_region,'XData');
    h_region_YData = get(h_region,'YData');
    h_region_YData = linspace(h_region_YData(1),h_region_YData(end),size(h_region_CData,1));
    ah_region_orig_axis = axis(ah_region);
    xlabel(ah_region, 'X (km)')
    ylabel(ah_region, 'Y (km)')
    
    % =======================================================================
    % Magnified raster image map
    ah_zoom = axes('parent',param.fig_hand(1));
    set(ah_zoom,'Position',[0.1 0.12 0.5*2/3 0.5])
    % Plot the image
    % Plot the image
    if size(RGB,3) == 1
      % Special code for Arctic sea ice maps
      RGB(RGB == 251) = 254;
      RGB(RGB == 255) = 255;
      colormap(cmap);
      h_zoom = mapshow(RGB, cmap, R,'Parent',ah_zoom);
    else
      h_zoom = mapshow(RGB, R,'Parent',ah_zoom);
    end
    h_zoom_CData = get(h_zoom,'CData');
    h_zoom_XData = get(h_zoom,'XData');
    h_zoom_YData = get(h_zoom,'YData');
    h_zoom_YData = linspace(h_zoom_YData(1),h_zoom_YData(end),size(h_zoom_CData,1));
    ah_zoom_orig_axis = axis(ah_zoom);
    xlabel(ah_zoom,'X (km)')
    ylabel(ah_zoom,'Y (km)')
    
    % =======================================================================
    % Contour map
    ah_map = axes('parent',param.fig_hand(1));
    [ant_X,ant_Y] = projfwd(proj,coastline.Lat,coastline.Lon);
    %h_coast = plot(ah_map, ant_X/1000, ant_Y/1000,'k');
    cur_idx = 1;
    for patch_idx = find(isnan(ant_X))
      h_coast = patch(ant_X(cur_idx:patch_idx-1)/1000, ant_Y(cur_idx:patch_idx-1)/1000,[0.7 0.7 0.7]);
      cur_idx = patch_idx + 1;
    end
    set(ah_map,'Color',[0.5 0.5 0.8]);
    axis(ah_map, map_axis)
    set(ah_map,'Position',contour_position)
    xlabel(ah_map,'X (km)')
    ylabel(ah_map,'Y (km)')
    
    set(param.fig_hand(1),'PaperOrientation','Portrait');
    set(param.fig_hand(1),'PaperPosition',[0.5 0.5 8 6]);
    
  elseif strcmpi(param.type,'singles')
    if length(param.fig_hand)>=1 && ishandle(param.fig_hand(1))
      clf(param.fig_hand(1));
    else
      if length(param.fig_hand)>=1 && isnumeric(param.fig_hand(1))
        param.fig_hand(1) = figure(param.fig_hand(1));
      else
        param.fig_hand(1) = figure();
      end
    end
    map_info.fig_hand(1) = param.fig_hand(1);
    set(param.fig_hand(1),'Renderer','painters');
    % =======================================================================
    % Large raster image map
    ah_region = axes('parent',param.fig_hand(1));
    % Get the projection information
    proj = geotiffinfo(geotiff_fn);
    % Read the image
    [RGB, R, tmp] = geotiffread(geotiff_fn);
    R = R/1e3;
    % Plot the image
    h_region = mapshow(ah_region, RGB, R);
    h_region_CData = get(h_region,'CData');
    h_region_XData = get(h_region,'XData');
    h_region_YData = get(h_region,'YData');
    h_region_YData = linspace(h_region_YData(1),h_region_YData(end),size(h_region_CData,1));
    ah_region_orig_axis = axis(ah_region);
    xlabel(ah_region, 'X (km)')
    ylabel(ah_region, 'Y (km)')
    set(param.fig_hand(1),'PaperOrientation','Landscape');
    set(param.fig_hand(1),'PaperPosition',[0.5 0.5 8 6]);
    
    if length(param.fig_hand)>=2 && ishandle(param.fig_hand(2))
      clf(param.fig_hand(2));
    else
      if length(param.fig_hand)>=2 && isnumeric(param.fig_hand(2))
        param.fig_hand(2) = figure(param.fig_hand(2));
      else
        param.fig_hand(2) = figure();
      end
    end
    map_info.fig_hand(2) = param.fig_hand(2);
    set(param.fig_hand(2),'Renderer','painters');
    % =======================================================================
    % Magnified raster image map
    ah_zoom = axes('parent',param.fig_hand(2));
    % Plot the image
    h_zoom = mapshow(ah_zoom, RGB, R);
    h_zoom_CData = get(h_zoom,'CData');
    h_zoom_XData = get(h_zoom,'XData');
    h_zoom_YData = get(h_zoom,'YData');
    h_zoom_YData = linspace(h_zoom_YData(1),h_zoom_YData(end),size(h_zoom_CData,1));
    ah_zoom_orig_axis = axis(ah_zoom);
    xlabel(ah_zoom,'X (km)')
    ylabel(ah_zoom,'Y (km)')
    set(param.fig_hand(2),'PaperOrientation','Landscape');
    set(param.fig_hand(2),'PaperPosition',[0.5 0.5 8 6]);
    
    if length(param.fig_hand)>=3 && ishandle(param.fig_hand(3))
      clf(param.fig_hand(3));
    else
      if length(param.fig_hand)>=3 && isnumeric(param.fig_hand(3))
        param.fig_hand(3) = figure(param.fig_hand(3));
      else
        param.fig_hand(3) = figure();
      end
    end
    map_info.fig_hand(3) = param.fig_hand(3);
    set(param.fig_hand(3),'Renderer','painters');
    % =======================================================================
    % Contour map
    ah_map = axes('parent',param.fig_hand(3));
    [ant_X,ant_Y] = projfwd(proj,coastline.Lat,coastline.Lon);
    h_coast = plot(ah_map, ant_X/1000, ant_Y/1000,'k');
    axis(ah_map, map_axis)
    xlabel(ah_map,'X (km)')
    ylabel(ah_map,'Y (km)')
    set(param.fig_hand(3),'PaperOrientation','Landscape');
    set(param.fig_hand(3),'PaperPosition',[0.5 0.5 8 6]);
    
  elseif strcmpi(param.type,'contour')
    if length(param.fig_hand)>=1 && ishandle(param.fig_hand(1))
      clf(param.fig_hand(1));
    else
      if length(param.fig_hand)>=1 && isnumeric(param.fig_hand(1))
        param.fig_hand(1) = figure(param.fig_hand(1));
      else
        param.fig_hand(1) = figure();
      end
    end
    map_info.fig_hand(1) = param.fig_hand(1);
    set(param.fig_hand(1),'Renderer','painters');
    % =======================================================================
    % Regional map
    
    % Get the projection information
    proj = geotiffinfo(geotiff_fn);
    ah_region = axes('parent',param.fig_hand(1));
    set(ah_region,'Position',[0.5000 0.1000 0.8*0.58 0.8000]);
    [ant_X,ant_Y] = projfwd(proj,coastline.Lat,coastline.Lon);
    cur_idx = 1;
    for patch_idx = find(isnan(ant_X))
      h_region = patch(ant_X(cur_idx:patch_idx-1)/1000, ant_Y(cur_idx:patch_idx-1)/1000,[0.7 0.7 0.7]);
      cur_idx = patch_idx + 1;
    end
    set(ah_region,'Color',[0.5 0.5 0.8]);
    axis(ah_region, map_axis)
    xlabel(ah_region,'X (km)')
    ylabel(ah_region,'Y (km)')
    set(param.fig_hand(1),'PaperOrientation','Portrait');
    set(param.fig_hand(1),'PaperPosition',[0.5 0.5 8 6]);
    
    % =======================================================================
    % Contour map
    
    % Get the projection information
    ah_map = axes('parent',param.fig_hand(1));
    [ant_X,ant_Y] = projfwd(proj,coastline.Lat,coastline.Lon);
    cur_idx = 1;
    for patch_idx = find(isnan(ant_X))
      patch(ant_X(cur_idx:patch_idx-1)/1000, ant_Y(cur_idx:patch_idx-1)/1000,[0.7 0.7 0.7]);
      cur_idx = patch_idx + 1;
    end
    set(ah_map,'Color',[0.5 0.5 0.8]);
    axis(ah_map, map_axis)
    set(ah_map,'Position',contour_position)
    xlabel(ah_map,'X (km)')
    ylabel(ah_map,'Y (km)')
    
    ah_zoom = 0;
    h_zoom = 0;
    
  elseif strcmpi(param.type,'contour-singles')
    if length(param.fig_hand)>=1 && ishandle(param.fig_hand(1))
      clf(param.fig_hand(1));
    else
      if length(param.fig_hand)>=1 && isnumeric(param.fig_hand(1))
        param.fig_hand(1) = figure(param.fig_hand(1));
      else
        param.fig_hand(1) = figure();
      end
    end
    map_info.fig_hand(1) = param.fig_hand(1);
    set(param.fig_hand(1),'Renderer','painters');
    % =======================================================================
    % Regional map
    
    % Get the projection information
    proj = geotiffinfo(geotiff_fn);
    ah_region = axes('parent',param.fig_hand(1));
    [ant_X,ant_Y] = projfwd(proj,coastline.Lat,coastline.Lon);
    cur_idx = 1;
    for patch_idx = find(isnan(ant_X))
      h_region = patch(ant_X(cur_idx:patch_idx-1)/1000, ant_Y(cur_idx:patch_idx-1)/1000,[0.7 0.7 0.7]);
      cur_idx = patch_idx + 1;
    end
    set(ah_region,'Color',[0.5 0.5 0.8]);
    axis(ah_region, map_axis)
    xlabel(ah_region,'X (km)')
    ylabel(ah_region,'Y (km)')
    set(param.fig_hand(1),'PaperOrientation','Landscape');
    set(param.fig_hand(1),'PaperPosition',[0.5 0.5 8 6]);
    
    if length(param.fig_hand)>=2 && ishandle(param.fig_hand(2))
      clf(param.fig_hand(2));
    else
      if length(param.fig_hand)>=2 && isnumeric(param.fig_hand(2))
        param.fig_hand(2) = figure(param.fig_hand(2));
      else
        param.fig_hand(2) = figure();
      end
    end
    map_info.fig_hand(2) = param.fig_hand(2);
    set(param.fig_hand(2),'Renderer','painters');
    % =======================================================================
    % Contour map
    
    % Get the projection information
    ah_map = axes('parent',param.fig_hand(2));
    [ant_X,ant_Y] = projfwd(proj,coastline.Lat,coastline.Lon);
    cur_idx = 1;
    for patch_idx = find(isnan(ant_X))
      patch(ant_X(cur_idx:patch_idx-1)/1000, ant_Y(cur_idx:patch_idx-1)/1000,[0.7 0.7 0.7]);
      cur_idx = patch_idx + 1;
    end
    set(ah_map,'Color',[0.5 0.5 0.8]);
    axis(ah_map, map_axis)
    xlabel(ah_map,'X (km)')
    ylabel(ah_map,'Y (km)')
    set(param.fig_hand(2),'PaperOrientation','Landscape');
    set(param.fig_hand(2),'PaperPosition',[0.5 0.5 8 6]);
    
    ah_zoom = 0;
    h_zoom = 0;
    
  else
    error('Unsupported type %s', param.type);
  end
  map_info.ah_map = ah_map;
  map_info.ah_zoom = ah_zoom;
  map_info.ah_region = ah_region;
  map_info.h_region = h_region;
  map_info.h_zoom = h_zoom;
  map_info.proj = proj;
  map_info.map_sub_title = map_sub_title;
  
elseif strcmpi(cmd,'plot')
  % ======================================================================
  % ======================================================================
  % PLOT COMMAND
  % ======================================================================
  % ======================================================================
  if strcmpi(param.type,'contour') || strcmpi(param.type, 'contour-singles')
    % Rename variables for convenience
    ah_map = map_info.ah_map;
    ah_region = map_info.ah_region;
    h_region = map_info.h_region;
    day_seg_x = param.day_seg_x/1000;
    day_seg_y = param.day_seg_y/1000;
    for frm_idx = 1:length(param.frame_X)
      frame_X{frm_idx} = param.frame_X{frm_idx}/1000;
      frame_Y{frm_idx} = param.frame_Y{frm_idx}/1000;
    end
    min_x = min(day_seg_x);
    max_x = max(day_seg_x);
    min_y = min(day_seg_y);
    max_y = max(day_seg_y);
    
    % =====================================================================
    % Plot the whole segment on to the regional raster map
    hold(ah_region,'on')
    if param.decimate_seg
      plot_idxs = unique(round(linspace(1,length(day_seg_x),1000)));
    else
      plot_idxs = 1:length(day_seg_x);
    end
    set(1,'Renderer','painters');
    h_seg = plot(ah_region,day_seg_x(plot_idxs), ...
      day_seg_y(plot_idxs),'b.');
    % Make window 10% larger than required, plus extra for legend
    mid_x = (max_x+min_x)/2;
    range_x = 1.1*0.5*(max_x-min_x);
    mid_y = (max_y+min_y)/2;
    range_y= 1.2*1.1*0.5*(max_y-min_y);
    if range_x < 100
      range_x = 100;
    end
    if range_y < 100*3/2
      range_y = 100*3/2;
    end
    if range_x < range_y*2/3
      range_x = range_y*2/3;
    elseif range_y < range_x*3/2
      range_y = range_x*3/2;
    end
    corners = [mid_x-range_x mid_x+range_x mid_y-range_y*1.4/1.2 mid_y+range_y/1.2];
    axis(ah_region, corners);
    if ~isempty(param.map_title)
      map_title = sprintf('%s\n%s', param.map_title, map_info.map_sub_title);
    else
      map_title = map_info.map_sub_title;
    end
    map_info.h_title = title(ah_region, map_title,'Interpreter','none');
    
    % =====================================================================
    % Plot the bounding box on the contour map
    hold(ah_map,'on');
    h_contour = plot(ah_map, [corners(1) corners(2) corners(2) corners(1) corners(1)], ...
      [corners(3) corners(3) corners(4) corners(4) corners(3)],'LineWidth',4);
    hold(ah_map,'off');
    
    % =====================================================================
    % Plot the current frame(s) on to the regional raster map
    hold(ah_region,'on')
    for frm_idx = 1:length(frame_X)
      if param.decimate_seg
        plot_idxs = unique(round(linspace(1,length(frame_X{frm_idx}),200)));
      else
        plot_idxs = 1:length(frame_X{frm_idx});
      end
      h_frm(frm_idx) = plot(ah_region,frame_X{frm_idx}(plot_idxs),frame_Y{frm_idx}(plot_idxs),'.r');
      h_start(frm_idx) = plot(ah_region,frame_X{frm_idx}(1),frame_Y{frm_idx}(1),'xg','LineWidth',4);
    end
    hold(ah_region,'off')
    
    ah_legend = legend([h_seg,h_frm(1) h_start(1)],'Segment','Frame','Start','Location','Southeast');
    
    map_info.h_seg = h_seg;
    map_info.h_contour = 0;
    map_info.h_frm_mag = 0;
    map_info.h_start_mag = 0;
    map_info.h_frm = h_frm;
    map_info.h_start = h_start;
  
  elseif strcmpi(param.type,'ASCAT')
    % Rename variables for convenience
    ah_map = map_info.ah_map;
    ah_region = map_info.ah_region;
    h_region = map_info.h_region;
    day_seg_x = param.day_seg_x/1000;
    day_seg_y = param.day_seg_y/1000;
    for frm_idx = 1:length(param.frame_X)
      frame_X{frm_idx} = param.frame_X{frm_idx}/1000;
      frame_Y{frm_idx} = param.frame_Y{frm_idx}/1000;
    end
    min_x = min(day_seg_x);
    max_x = max(day_seg_x);
    min_y = min(day_seg_y);
    max_y = max(day_seg_y);
    
    % =====================================================================
    % Plot the whole segment on to the regional raster map
    hold(ah_region,'on')
    if param.decimate_seg
      plot_idxs = unique(round(linspace(1,length(day_seg_x),1000)));
    else
      plot_idxs = 1:length(day_seg_x);
    end
    h_seg = plot(ah_region,day_seg_x(plot_idxs), ...
      day_seg_y(plot_idxs),'k.');
    % Make window 10% larger than required, plus extra for legend
    mid_x = (max_x+min_x)/2;
    range_x = 1.1*0.5*(max_x-min_x);
    mid_y = (max_y+min_y)/2;
    range_y= 1.2*1.1*0.5*(max_y-min_y);
    if range_x < 100
      range_x = 100;
    end
    if range_y < 100*3/2
      range_y = 100*3/2;
    end
    if range_x < range_y*2/3
      range_x = range_y*2/3;
    elseif range_y < range_x*3/2
      range_y = range_x*3/2;
    end
    corners = [mid_x-range_x mid_x+range_x mid_y-range_y*1.4/1.2 mid_y+range_y/1.2];
    axis(ah_region, corners);
    set(ah_region, 'Position', [0.53 0.1 0.8*2/3 0.8]);
    if ~isempty(param.map_title)
      map_title = sprintf('%s\n%s', param.map_title, map_info.map_sub_title);
    else
      map_title = map_info.map_sub_title;
    end
    map_info.h_title = title(ah_region, map_title,'Interpreter','none');
    
    % =====================================================================
    % Plot the bounding box on the contour map
    hold(ah_map,'on');
    h_contour = plot(ah_map, [corners(1) corners(2) corners(2) corners(1) corners(1)], ...
      [corners(3) corners(3) corners(4) corners(4) corners(3)],'LineWidth',4, 'Color', 'k');
    hold(ah_map,'off');
    
    % =====================================================================
    % Plot the current frame(s) on to the regional raster map
    hold(ah_region,'on')
    for frm_idx = 1:length(frame_X)
      if param.decimate_seg
        plot_idxs = unique(round(linspace(1,length(frame_X{frm_idx}),200)));
      else
        plot_idxs = 1:length(frame_X{frm_idx});
      end
      h_frm(frm_idx) = plot(ah_region,frame_X{frm_idx}(plot_idxs),frame_Y{frm_idx}(plot_idxs),'.r');
      h_start(frm_idx) = plot(ah_region,frame_X{frm_idx}(1),frame_Y{frm_idx}(1),'xg','LineWidth',4);
    end
    hold(ah_region,'off')
    
    ah_legend = legend([h_seg,h_frm(1) h_start(1)],'Segment','Frame','Start','Location','Southeast');
    
    map_info.h_seg = h_seg;
    map_info.h_contour = 0;
    map_info.h_frm_mag = 0;
    map_info.h_start_mag = 0;
    map_info.h_frm = h_frm;
    map_info.h_start = h_start;
         
  elseif strcmpi(param.type,'geotiff')
    % Rename variables for convenience
    ah_map = map_info.ah_map;
    ah_region = map_info.ah_region;
    h_region = map_info.h_region;
    day_seg_x = param.day_seg_x/1000;
    day_seg_y = param.day_seg_y/1000;
    for frm_idx = 1:length(param.frame_X)
      frame_X{frm_idx} = param.frame_X{frm_idx}/1000;
      frame_Y{frm_idx} = param.frame_Y{frm_idx}/1000;
    end
    min_x = min(frame_X{frm_idx});
    max_x = max(frame_X{frm_idx});
    min_y = min(frame_Y{frm_idx});
    max_y = max(frame_Y{frm_idx});
    
    % =====================================================================
    % Plot the whole segment on to the regional raster map
    hold(ah_region,'on')
    if param.decimate_seg
      plot_idxs = unique(round(linspace(1,length(day_seg_x),1000)));
    else
      plot_idxs = 1:length(day_seg_x);
    end
    h_seg = plot(ah_region,day_seg_x(plot_idxs), ...
      day_seg_y(plot_idxs),'k.');
    % Make window an extra 3 km on each side
    mid_x = (max_x+min_x)/2;
    range_x = (max_x-min_x)+6;
    mid_y = (max_y+min_y)/2;
    range_y= (max_y-min_y)+6;
    % Make sure aspect ratio is good
    if range_x * 6/8 > range_y
      range_y = range_x*6/8;
    else
      range_x = range_y*8/6;
    end
    corners = [mid_x-range_x/2 mid_x+range_x/2 mid_y-range_y/2 mid_y+range_y/2];
    axis(ah_region, corners);
    % A hack to get Matlab to display axis properly and use up the whole
    % figure space
    axis(ah_region,'normal'); axis(ah_region,'equal');
    axis(ah_region, corners);
    
    if ~isempty(param.map_title)
      map_title = sprintf('%s\n%s', param.map_title, map_info.map_sub_title);
    else
      map_title = map_info.map_sub_title;
    end
    map_info.h_title = title(ah_region, map_title,'Interpreter','none');
        
    % =====================================================================
    % Plot the current frame(s) on to the regional raster map
    hold(ah_region,'on')
    for frm_idx = 1:length(frame_X)
      if param.decimate_seg
        plot_idxs = unique(round(linspace(1,length(frame_X{frm_idx}),200)));
      else
        plot_idxs = 1:length(frame_X{frm_idx});
      end
      h_frm(frm_idx)= plot(ah_region,frame_X{frm_idx}(plot_idxs),frame_Y{frm_idx}(plot_idxs),'.r');
      h_start(frm_idx) = plot(ah_region,frame_X{frm_idx}(1),frame_Y{frm_idx}(1),'xg','LineWidth',4);
    end
    hold(ah_region,'off')
    
    ah_legend = [];
    
    map_info.h_seg = h_seg;
    map_info.h_contour = 0;
    map_info.h_frm_mag = 0;
    map_info.h_start_mag = 0;
    map_info.h_frm = h_frm;
    map_info.h_start = h_start;
       
  elseif strcmpi(param.type,'combined') || strcmpi(param.type, 'singles')
    % Rename variables for convenience
    ah_map = map_info.ah_map;
    ah_zoom = map_info.ah_zoom;
    ah_region = map_info.ah_region;
    h_region = map_info.h_region;
    h_zoom = map_info.h_zoom;
    day_seg_x = param.day_seg_x/1000;
    day_seg_y = param.day_seg_y/1000;
    for frm_idx = 1:length(param.frame_X)
      frame_X{frm_idx} = param.frame_X{frm_idx}/1000;
      frame_Y{frm_idx} = param.frame_Y{frm_idx}/1000;
    end
    min_x = min(day_seg_x);
    max_x = max(day_seg_x);
    min_y = min(day_seg_y);
    max_y = max(day_seg_y);
    
    % =====================================================================
    % Plot the whole segment on to the regional raster map
    hold(ah_region,'on')
    if param.decimate_seg
      plot_idxs = unique(round(linspace(1,length(day_seg_x),1000)));
    else
      plot_idxs = 1:length(day_seg_x);
    end
    h_seg = plot(ah_region,day_seg_x(plot_idxs), ...
      day_seg_y(plot_idxs),'b.');
    % Make window 10% larger than required, plus extra for legend
    mid_x = (max_x+min_x)/2;
    range_x = 1.1*0.5*(max_x-min_x);
    mid_y = (max_y+min_y)/2;
    range_y= 1.2*1.1*0.5*(max_y-min_y);
    if range_x < 100
      range_x = 100;
    end
    if range_y < 100*3/2
      range_y = 100*3/2;
    end
    if range_x < range_y*2/3
      range_x = range_y*2/3;
    elseif range_y < range_x*3/2
      range_y = range_x*3/2;
    end
    corners = [mid_x-range_x mid_x+range_x mid_y-range_y*1.4/1.2 mid_y+range_y/1.2];
    axis(ah_region, corners);
    if ~isempty(param.map_title)
      map_title = sprintf('%s\n%s', param.map_title, map_info.map_sub_title);
    else
      map_title = map_info.map_sub_title;
    end
    map_info.h_title = title(ah_region, map_title,'Interpreter','none');
    
    % =====================================================================
    % Plot the bounding box on the contour map
    hold(ah_map,'on');
    h_contour = plot(ah_map, [corners(1) corners(2) corners(2) corners(1) corners(1)], ...
      [corners(3) corners(3) corners(4) corners(4) corners(3)],'LineWidth',4);
    hold(ah_map,'off');
    
    % =====================================================================
    % Plot the current frame(s) on to the regional raster map
    hold(ah_region,'on')
    for frm_idx = 1:length(frame_X)
      if param.decimate_seg
        plot_idxs = unique(round(linspace(1,length(frame_X{frm_idx}),200)));
      else
        plot_idxs = 1:length(frame_X{frm_idx});
      end
      h_frm(frm_idx) = plot(ah_region,frame_X{frm_idx}(plot_idxs),frame_Y{frm_idx}(plot_idxs),'.r');
      h_start(frm_idx) = plot(ah_region,frame_X{frm_idx}(1),frame_Y{frm_idx}(1),'xg','LineWidth',4);
    end
    hold(ah_region,'off')
    
    ah_legend = legend(ah_region,'Segment','Frame','Start','Location','Southeast');
    
    % =====================================================================
    % Plot the current/first frame on to the magnified raster map
    frame_X = frame_X{1};
    frame_Y = frame_Y{1};
    if param.decimate_seg
      plot_idxs = unique(round(linspace(1,length(frame_X),200)));
    else
      plot_idxs = 1:length(frame_X);
    end
   
    hold(ah_zoom,'on')
    h_frm_mag = plot(ah_zoom,frame_X(plot_idxs),frame_Y(plot_idxs),'.r');
    h_start_mag = plot(ah_zoom,frame_X(1),frame_Y(1),'xg','LineWidth',4);
    hold(ah_zoom,'off')

    mid_x = (max(frame_X)+min(frame_X))/2;
    range_x = 1.05*0.5*(max(frame_X)-min(frame_X));
    mid_y = (max(frame_Y)+min(frame_Y))/2;
    range_y= 1.05*0.5*(max(frame_Y)-min(frame_Y));
    if range_x < range_y*2/3
      range_x = range_y*2/3;
    elseif range_y < range_x*3/2
      range_y = range_x*3/2;
    end
    if range_x == 0
      range_x = 1;
    end
    if range_y == 0
      range_y = 1;
    end
    corners = [mid_x-range_x mid_x+range_x mid_y-range_y mid_y+range_y];
    axis(ah_zoom, corners);
    
    if param.resample
      % Resample region map because it will be saved in the pdf/eps/fig at full
      % resolution otherwise
      PaperPosition = get(map_info.fig_hand(1),'PaperPosition');
      Position = get(ah_region,'Position');
      resolution = 200;
      pixels = round(PaperPosition(3:4) .* Position(3:4) * resolution);
      ah_region_axis = axis(ah_region);
      new_XData = linspace(ah_region_axis(1),ah_region_axis(2),pixels(1));
      new_YData = linspace(ah_region_axis(3),ah_region_axis(4),pixels(2));
      new_CData = zeros(pixels(2),pixels(1),3,'uint8');
      warning off
      for col = 1:size(new_CData,2)
        % Determine the column in the source image to use for interpolation
        col_source = round(interp1(h_region_XData([1 end]), ...
          [1 size(h_region_CData,2)],new_XData(col)));
        if isnan(col_source)
          new_CData(:,col,:) = 255;
        else
          % Interpolate each column of the output resampled image
          for band = 1:3
            newCol = interp1(h_region_YData, ...
              double(h_region_CData(:,col_source,band)),new_YData);
            newCol(isnan(newCol)) = 255;
            new_CData(:,col,band) = newCol;
          end
        end
      end
      warning on
      
      set(h_region,'CData',new_CData);
      set(h_region,'XData',new_XData);
      set(h_region,'YData',new_YData);
      
      % Resample zoom map because it will be saved in the pdf/eps at full
      % resolution otherwise
      PaperPosition = get(map_info.fig_hand(1),'PaperPosition');
      Position = get(ah_zoom,'Position');
      resolution = 200;
      pixels = round(PaperPosition(3:4) .* Position(3:4) * resolution);
      ah_zoom_axis = axis(ah_zoom);
      new_XData = linspace(ah_zoom_axis(1),ah_zoom_axis(2),pixels(1));
      new_YData = linspace(ah_zoom_axis(3),ah_zoom_axis(4),pixels(2));
      new_CData = zeros(pixels(2),pixels(1),3,'uint8');
      warning off
      for col = 1:size(new_CData,2)
        % Determine the column in the source image to use for interpolation
        col_source = round(interp1(h_zoom_XData([1 end]), ...
          [1 size(h_zoom_CData,2)],new_XData(col)));
        if isnan(col_source)
          new_CData(:,col,:) = 255;
        else
          % Interpolate each column of the output resampled image
          for band = 1:3
            newCol = interp1(h_zoom_YData, ...
              double(h_zoom_CData(:,col_source,band)),new_YData);
            newCol(isnan(newCol)) = 255;
            new_CData(:,col,band) = newCol;
          end
        end
      end
      warning on
      
      set(h_zoom,'CData',new_CData);
      set(h_zoom,'XData',new_XData);
      set(h_zoom,'YData',new_YData);
    end
    
    map_info.h_seg = h_seg;
    map_info.h_contour = h_contour;
    map_info.h_frm_mag = h_frm_mag;
    map_info.h_start_mag = h_start_mag;
    map_info.h_frm = h_frm;
    map_info.h_start = h_start;
  end
  
elseif strcmpi(cmd,'delete')
  if isfield(map_info,'h_seg') && ishandle(map_info.h_seg)
    % If any handle exists, then they all should exist, so delete them
    % all
    delete(map_info.h_seg);
    if map_info.h_contour ~= 0
      delete(map_info.h_contour);
    end
    if map_info.h_frm_mag ~= 0
      delete(map_info.h_frm_mag);
    end
    if map_info.h_start_mag ~= 0
      delete(map_info.h_start_mag);
    end
    delete(map_info.h_frm);
    delete(map_info.h_start);
  end
  
else
  error('Unsupported command %s', cmd);
end

return;






