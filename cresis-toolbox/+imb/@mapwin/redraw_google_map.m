% Called by imb/mapwin/key_press to pan around the map or zoom in/out in the picker
% Updates the Map in the figure and the axis data
function redraw_google_map(obj, x_min, x_max, y_min, y_max)
  
  %% Updating map depending on whether user wants to pan or zoom
  if ~obj.googleObj.zoom_in_out == 0 
    % If zoom in/out is selected
    
    % Getting the midpoint of the zoom box
    x_m = (x_min+x_max)/2;
    y_m = (y_min+y_max)/2;

    % Converting the midpoint from world coordinates to lat lon
    [new_c_lat, new_c_lon] = imb.world_to_latlon(x_m, y_m);

    % Determining the new zoom level
    area_zoom_box = (x_max - x_min)*(y_max - y_min);
    map_area = (obj.map_panel.h_axes.XLim(2) - obj.map_panel.h_axes.XLim(1)) * (obj.map_panel.h_axes.YLim(2) - obj.map_panel.h_axes.YLim(1));
    r = map_area/area_zoom_box;
    factor = 1;

    if isinf(r)
      % User didn't draw a zoom box. It was a left click to zoom in/out
      % Incrementing or decrementing the zoom level
      if obj.googleObj.zoom_in_out == 1
        zoom = obj.googleObj.zoom+1;
      else
        zoom = obj.googleObj.zoom-1;
      end
    else
      % r is a finite value implies the user drew a zoom box
      % Determining the zoom level relative to r
      if r>2 && r<60
        factor = 2;
      elseif r>60 && r<130
        factor = 3;
      elseif r>130
        factor = 4;
      end
      zoom = obj.googleObj.zoom+factor;  
    end
    
    % Resetting the zoom_in_out attribute
    obj.googleObj.zoom_in_out = 0;
  else
    % Pan is selected
    
    % Zoom level remains the same
    zoom = obj.googleObj.zoom;
    
    % 1 - Up
    % 2 - Right
    % 3 - Down
    % 4 - Left
    % Calculating the new lat, lon depending on the arrow key pressed
    if obj.googleObj.pan == 1
      new_c_lat = obj.googleObj.c_lat+1;
      new_c_lon = obj.googleObj.c_lon;
    elseif obj.googleObj.pan == 2
      new_c_lat = obj.googleObj.c_lat;
      new_c_lon = obj.googleObj.c_lon+1;
    elseif obj.googleObj.pan == 3
      new_c_lat = obj.googleObj.c_lat-1;
      new_c_lon = obj.googleObj.c_lon;
    elseif obj.googleObj.pan == 4
      new_c_lat = obj.googleObj.c_lat;
      new_c_lon = obj.googleObj.c_lon-1;
    end
  end
  
  % Getting the new google map centered at the new lat lon
  obj.googleObj = imb.request_google_map(obj.googleObj, new_c_lat, new_c_lon, zoom);
  A = obj.googleObj.A;
  
  % Updating the map and the axes
  
  A = flipud(A);
  
  % New Axis limits and values
  xaxis = [obj.googleObj.bottom_left_wc_x obj.googleObj.bottom_right_wc_x];
  yaxis = [obj.googleObj.bottom_left_wc_y obj.googleObj.top_left_wc_y];
  
  set(obj.map_panel.h_image,'XData', xaxis, ...
    'YData', yaxis, ...
    'CData', A, ...
    'Visible', 'on');

  set(obj.map_panel.h_axes, 'Xlim', sort(xaxis([1 end])), ...
  'Ylim', sort(yaxis([1 end])), ...
  'YDir', 'reverse', ...
  'Visible', 'on');
end