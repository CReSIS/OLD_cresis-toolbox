% Converts given pixel coordinates to world coordinates
function [wc_x, wc_y] = pixel_to_world(px_x, px_y, zoom)
  wc_x = [];
  wc_y = [];
  for idx = 1:length(px_x)
    wc_x(idx) = (px_x(idx)/(2^zoom));
    wc_y(idx) = (px_y(idx)/(2^zoom));
  end
end
