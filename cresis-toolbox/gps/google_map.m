classdef google_map
  % google_map class
  %
  % World coordinates description:
  %  x ranges from 0 to 256 and corresponds linearly to -180 to +180
  %   longitude
  %  y ranges from 0 to 256 and corresponds to -85.051128779806632 to +85.051128779806632
  %  https://developers.google.com/maps/documentation/javascript/coordinates#world-coordinates
  %
  % Authors: Rohan Choudhari, John Paden
  
  %% public properties
  properties
    key = 'AIzaSyCNexiP6WcIda8ZEa2MnwznWrGotDoLu0w';
    tile_size = 256;
    w = 1280;
    h = 1280;
    scale = 2;
    zoom_range = [0 18];
    wms_obj;
    zoom;
    c_lat;
    c_lon;
  end
  
  methods
    %% google_map constructor
    function obj = google_map(data)
      obj.wms_obj = WebMapServer('http://maps.googleapis.com/maps/api');
    end
    
    %% google_map destructor
    function delete(obj)
    end
    
    %% request_google_map
    % wc_x_min: world coordinate scalar numeric x-axis minimum
    % wc_x_max: world coordinate scalar numeric x-axis maximum
    % wc_y_min: world coordinate scalar numeric y-axis minimum
    % wc_y_max: world coordinate scalar numeric y-axis maximum
    %
    % Returns the Google map with the most precision that captures the
    % extent specified by the input. The output will be clipped to the
    % region specified by the input.
    function [A,x_axis,y_axis] = request_google_map(obj, wc_x_min, wc_x_max, wc_y_min, wc_y_max)
      if wc_x_max < wc_x_min
        tmp = wc_x_min;
        wc_x_min = wc_x_max;
        wc_x_max = tmp;
      end
      if wc_y_max < wc_y_min
        tmp = wc_y_min;
        wc_y_min = wc_y_max;
        wc_y_max = tmp;
      end
      if wc_y_min < 0
        wc_y_min = 0;
      end
      if wc_y_max > 256
        wc_y_max = 256;
      end
      if wc_x_min < -256
        wc_x_min = -256;
      end
      if wc_x_max > 512
        wc_x_max = 512;
      end
      % Determine the zoom level
      x_extent = wc_x_max-wc_x_min;
      y_extent = wc_y_max-wc_y_min;
      max_extent = max(x_extent,y_extent);
      obj.zoom = floor(log2(256/max_extent));
      if obj.zoom < obj.zoom_range(1)
        obj.zoom = obj.zoom_range(1);
      end
      if obj.zoom > obj.zoom_range(2)
        obj.zoom = obj.zoom_range(2);
      end
      
      % Determine the center coordinate
      c_wc_x = (wc_x_min+wc_x_max)/2;
      c_wc_y = (wc_y_min+wc_y_max)/2;
      
      % Correct the center coordinate for what Google will really return.
      % (All requested center coordinates are rounded to the nearest pixel
      % at the specified zoom level.)
      dwc = 1/2^obj.zoom;
      c_wc_x = round(c_wc_x/dwc)*dwc;
      c_wc_y = round(c_wc_y/dwc)*dwc;
      
      % Find the lat,lon for the center coordinate
      [obj.c_lat,obj.c_lon] = google_map.world_to_latlon(c_wc_x,c_wc_y);
      
      % Get Google map image (A is obj.h x obj.w x 3 RGB image)
      A = request_google_mapc(obj, obj.c_lat, obj.c_lon, obj.zoom);
      
      % Create the corresponding x and y axes for the image A
      dwc = 1/2^obj.zoom/obj.scale;
      x_axis = c_wc_x + (-obj.w/2:obj.w/2-1)*dwc;
      y_axis = c_wc_y + (-obj.h/2:obj.h/2-1).'*dwc;
      
      % Truncate image to the input region
      x_mask = x_axis >= wc_x_min & x_axis <= wc_x_max;
      y_mask = y_axis >= wc_y_min & y_axis <= wc_y_max;
      A = A(y_mask,x_mask,:);
      x_axis = x_axis(x_mask);
      y_axis = y_axis(y_mask);
    end
    
    %% request_google_mapc
    % c_lat: scalar numeric latitude, center pixel coordinate
    % c_lon: scalar numeric longitude, center pixel coordinate
    % zoom: nonnegative integer representing zoom level (0 = no zoom)
    % A: image returned from Google
    % Google will only use 6 decimals of accuracy for c_lat and c_lon
    function A = request_google_mapc(obj, c_lat, c_lon, zoom)
      % Request PNG map (PNG is the default format)
      % - Only the first six digits of precision for lon and lat are used
      %   by Google
      url = sprintf('http://maps.googleapis.com/maps/api/staticmap?center=%0.6f,%0.6f&zoom=%d&scale=%d&size=%dx%d&maptype=satellite&key=%s', ...
        c_lat, c_lon, zoom, obj.scale, obj.w, obj.h, obj.key);
      A = obj.wms_obj.getMap(url);
    end
     
  end
  
  methods(Static)
    
    %% greenland
    function [wc_x_min,wc_x_max,wc_y_min,wc_y_max] = greenland()
      wc_x_min = 48.75;
      wc_x_max = 150;
      wc_y_min = 256-247.5;
      wc_y_max = 256-167.75;
    end
    
    %% antarctica
    function [wc_x_min,wc_x_max,wc_y_min,wc_y_max] = antarctica()
      wc_x_min = 29.25;
      wc_x_max = 130.25;
      wc_y_min = 256-89.5;
      wc_y_max = 256-9.75;
    end
   
    %% latlon_to_world
    function [wc_x, wc_y] = latlon_to_world(lat, lon)
      wc_x = [];
      wc_y = [];
      for idx = 1:length(lat)
        siny = sin(lat(idx)*pi/180);
        siny = min(max(siny, -0.9999), 0.9999);
        wc_x(idx) = 256*(0.5 + lon(idx)/360);
        wc_y(idx) = 256*(0.5 - log((1 + siny)/(1 - siny))/(4*pi));
      end
      wc_x = mod(wc_x,256);
    end
    
    %% world_to_latlon
    function [lat, lon] = world_to_latlon(wc_x, wc_y)
      lat = [];
      lon = [];
      for idx = 1:length(wc_x)
        lat(idx) = -(((wc_y(idx)/256) - 0.5)*4*pi)/2;
        lat(idx) = asind(tanh(lat(idx)));
        
        lon(idx) = ((wc_x(idx)/256)-0.5)*360;
      end
    end
    
  end

end
