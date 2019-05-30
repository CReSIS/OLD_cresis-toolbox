% Converts the given world coordinates to pixel coordinates
function [px_x, px_y] = world_to_pixel(wc_x, wc_y, zoom)
  px_x = [];
  px_y = [];
  for idx = 1:length(wc_x)
    px_x(idx) = floor(wc_x(idx)*(2^zoom));
    px_y(idx) = floor(wc_y(idx)*(2^zoom));
  end
end
