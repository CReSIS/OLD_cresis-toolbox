% Converts given lat lon values to world coordinates
function [wc_x, wc_y] = latlon_to_world(lat, lon)
  wc_x = [];
  wc_y = [];
  for idx = 1:length(lat)
    siny = sin(lat(idx)*pi/180);
    siny = min(max(siny, -0.9999), 0.9999);
    wc_x(idx) = 256*(0.5 + lon(idx)/360);
    wc_y(idx) = 256*(0.5 - log((1 + siny)/(1 - siny))/(4*pi));
  end
end
