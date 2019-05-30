% Converts given world coordinates to lat lon
function [lat, lon] = world_to_latlon(wc_x, wc_y)
      lat = [];
      lon = [];
      for idx = 1:length(wc_x)
        lat(idx) = -(((wc_y(idx)/256) - 0.5)*4*pi)/2;
        lat(idx) = asind(tanh(lat(idx)));
        
        lon(idx) = ((wc_x(idx)/256)-0.5)*360;
      end
end