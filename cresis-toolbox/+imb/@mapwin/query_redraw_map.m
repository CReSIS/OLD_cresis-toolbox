function query_redraw_map(obj,x_min,x_max,y_min,y_max)
% mapwin.query_redraw_map(obj,x_min,x_max,y_min,y_max)
%
% mapwin class support function
% Called anytime the map window needs to get a new map from the map server.
% E.g. button_up and button_scroll zoom and search functions, key_press
%   pan function, search_callback.m, update_echowin_flightlines.m and
%   map_update_cursors_callback.m in the event that the map must be
%   adjusted to keep the flightline/cursor in view
%
% x_min, x_max, y_min, y_max: Scalars that specify the new map limits
%   (in km) to get

%% Scale the requested limits to meters
obj.cur_request.XLim = sort([x_min x_max]*obj.map_scale);
obj.cur_request.YLim = sort([y_min y_max]*obj.map_scale);

%% Force requested limits to maintain the aspect ratio of the figure
old_u = get(obj.map_panel.h_axes,'units');
set(obj.map_panel.h_axes,'Units','pixels')
PixelBounds = round(get(obj.map_panel.h_axes,'Position'));
set(obj.map_panel.h_axes,'Position',PixelBounds);
height = round((PixelBounds(4))*1 - 0);
width = round((PixelBounds(3))*1 - 0);
aspect_ratio = height/width;
set(obj.map_panel.h_axes,'units',old_u);

% Grow the limits so that the requested region as as big or bigger
% than the requested limits
if aspect_ratio*diff(obj.cur_request.XLim) > diff(obj.cur_request.YLim)
  growth = aspect_ratio*diff(obj.cur_request.XLim) - diff(obj.cur_request.YLim);
  obj.cur_request.YLim(1) = obj.cur_request.YLim(1) - growth/2;
  obj.cur_request.YLim(2) = obj.cur_request.YLim(2) + growth/2;
else
  growth = diff(obj.cur_request.YLim)/aspect_ratio - diff(obj.cur_request.XLim);
  obj.cur_request.XLim(1) = obj.cur_request.XLim(1) - growth/2;
  obj.cur_request.XLim(2) = obj.cur_request.XLim(2) + growth/2;
end

if obj.map_source == 0
  %% Build the new WMS query, submit it and then retrieve the result
  obj.cur_request.ImageHeight =  height;
  obj.cur_request.ImageWidth  = width;
  modrequest = strcat(obj.cur_request.RequestURL,'&viewparams=',obj.seasons_modrequest,obj.season_group_ids_modrequest);
  A = obj.wms.getMap(modrequest);
  R = obj.cur_request.RasterRef;
  R = R/1e3;
  
  %% Redraw the map
  x_data = R(3,1) + [0 size(A,2)*R(2,1)];
  y_data = R(3,2) + [0 size(A,1)*R(1,2)];
  
else
  [A,x_data,y_data] = obj.google_map.request_google_map(obj.cur_request.XLim(1), obj.cur_request.XLim(2), 256-obj.cur_request.YLim(1), 256-obj.cur_request.YLim(2));
  A = flipud(A);
  y_data = sort(256-y_data);

end

set(obj.map_panel.h_image,'XData',x_data);
set(obj.map_panel.h_image,'YData',y_data);
set(obj.map_panel.h_image,'CData',A);
set(obj.map_panel.h_axes, {'XLim','YLim'}, {sort(x_data([1 end])) sort(y_data([1 end]))});

return;