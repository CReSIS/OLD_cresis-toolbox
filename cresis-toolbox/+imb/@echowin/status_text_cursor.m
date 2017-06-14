function [status_str] = status_text_cursor(obj)
% [status_str] = echowin.status_text_cursor(obj)
%
% Builds a string containing information about the current cursor position 
% to be printed to the echogram window's status text bar. 

lat = interp1(obj.eg.gps_time,obj.eg.latitude,obj.cursor.gps_time,'linear','extrap');
lon = interp1(obj.eg.gps_time,obj.eg.longitude,obj.cursor.gps_time,'linear','extrap');

date_str = [datestr(epoch_to_datenum(obj.cursor.gps_time),'yyyy-mm-dd HH:MM:SS') sprintf('.%02.0f',mod(obj.cursor.gps_time,1)*100)];

% Current layers
cur_layers = find(obj.left_panel.layer_panel.selected_layers).';

if isempty(cur_layers)
  status_str = sprintf('%.6f N, %.6f E, X:%.3f km, Y:%.3f km, %s',...
    lat,lon,obj.cursor.x,obj.cursor.y,date_str);
else
  physical_constants;
  elev = interp1(obj.eg.gps_time,obj.eg.elevation,obj.cursor.gps_time,'linear','extrap');
  [~,unique_idxs] = unique(obj.eg.layer.x{1});
  warning off;
  surf_y = interp1(obj.eg.layer.x{1}(unique_idxs), ...
    obj.eg.layer.y{1}(unique_idxs),obj.cursor.gps_time);
  [~,unique_idxs] = unique(obj.eg.layer.x{cur_layers(1)});
  try
    layer_y = interp1(obj.eg.layer.x{cur_layers(1)}(unique_idxs), ...
      obj.eg.layer.y{cur_layers(1)}(unique_idxs),obj.cursor.gps_time);
    if ~all(isfinite(layer_y))
      error('');
    end
    depth = (layer_y-surf_y)*c/2/sqrt(er_ice);
    status_str = sprintf('%.6f N, %.6f E, X:%.3f km, Y:%.3f km, %s, %.2f Elevation, %.2f Depth',...
      lat,lon,obj.cursor.x,obj.cursor.y,date_str,elev - surf_y*c/2 - depth, depth);
  catch
    % There may be insufficient layer data to interpolate, so just print
    % out all the non-layer information
    status_str = sprintf('%.6f N, %.6f E, X:%.3f km, Y:%.3f km, %s',...
      lat,lon,obj.cursor.x,obj.cursor.y,date_str);
  end
  warning on;
end

return;