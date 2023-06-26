% =====================================================================
% Support function for finding min and max x/y coords from init/final click
% Also ensures that mins/maxes don't exceed current axis limits
% =====================================================================
function [x_min x_max y_min y_max] = sort_clicks(x_ext,y_ext,xi,yi,xf,yf)
  
  x_min = min(xi,xf);
  if x_min < x_ext(1)
    x_min = x_ext(1);
  end
  x_max = max(xi,xf);
  if x_max > x_ext(end)
    x_max = x_ext(end);
  end
  y_min = min(yi,yf);
  if y_min < y_ext(1)
    y_min = y_ext(1);
  end
  y_max = max(yi,yf);
  if y_max > y_ext(end)
    y_max = y_ext(end);
  end

return;
