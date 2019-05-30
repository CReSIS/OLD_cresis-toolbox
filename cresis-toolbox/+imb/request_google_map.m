% request_google_map returns a Google Map centered at (c_lat, c_lon) amd at
% zoom level 'zoom'
function obj = request_google_map(obj, c_lat, c_lon, zoom)
  
  obj.c_lat = c_lat;
  obj.c_lon = c_lon;
  obj.zoom = zoom;

  %% Getting the raster map from Google Maps API
  wms_obj = WebMapServer('http://maps.googleapis.com/maps/api');
  url = strcat('http://maps.googleapis.com/maps/api/staticmap?center=',num2str(c_lat),',',num2str(c_lon),'&zoom=',num2str(zoom),'&scale=2&size=', num2str(obj.w),'x',num2str(obj.h),'&maptype=satellite&key=',obj.key);
  obj.A = wms_obj.getMap(url);

  %% Storing important coordinates in the object
  
  % Lat, lon in world and pixel coordinates
  [obj.c_wc_x, obj.c_wc_y] = imb.latlon_to_world(obj.c_lat, obj.c_lon);
  [obj.c_px_x, obj.c_px_y] = imb.world_to_pixel(obj.c_wc_x, obj.c_wc_y, obj.zoom);
  
  % Calculating the pixel coordinates of the map vertices
  obj.top_left_px_x = obj.c_px_x-320;
  obj.top_left_px_y = obj.c_px_y-320;

  obj.bottom_left_px_x = obj.c_px_x-320;
  obj.bottom_left_px_y = obj.c_px_y+320;

  obj.top_right_px_x = obj.c_px_x+320;
  obj.top_right_px_y = obj.c_px_y-320;

  obj.bottom_right_px_x = obj.c_px_x+320;
  obj.bottom_right_px_y = obj.c_px_y+320;
  
  % Converting the above calculated pixel coordinates of the vertices to
  % world coordinates
  [obj.top_left_wc_x, obj.top_left_wc_y] = imb.pixel_to_world(obj.top_left_px_x, obj.top_left_px_y, obj.zoom);
  [obj.bottom_left_wc_x, obj.bottom_left_wc_y] = imb.pixel_to_world(obj.bottom_left_px_x, obj.bottom_left_px_y, obj.zoom);
  [obj.top_right_wc_x , obj.top_right_wc_y] = imb.pixel_to_world(obj.top_right_px_x, obj.top_right_px_y, obj.zoom);
  [obj.bottom_right_wc_x, obj.bottom_right_wc_y] = imb.pixel_to_world(obj.bottom_right_px_x, obj.bottom_right_px_y, obj.zoom);

end