function picker_map_button(src,event)

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

[x,y,but] = get_mouse_info(hui.mapfig.handle,hui.mapfig.axes.handle);
fprintf('Map: x = %.3f, y = %.3f, but = %d\n', x, y, but);

if but == 1
  % ===================================================================
  % Left mouse button: select a frame
  
  % Find the closest point
  [min_dist min_idx] = find_closest_point(x,y);
  picker_map(2, gCtrl.source.F(min_idx));

elseif but ==4
  picker_pick(1);

else  
  % ===================================================================
  % Right mouse button: move cursors
  
  % ------------------------------------------------------------------
  % Find the closest point for the pick frame and put marker on it
  [min_dist min_idx] = find_closest_point(x,y,gCtrl.source.cur_pick);
  
  % Update the cursor data structures
  gCtrl.pick.cur_idx = min_idx;
  gCtrl.pick.cursor = gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time(gCtrl.pick.cur_idx);
  fprintf('  Pick: GPS %s\n', datestr(epoch_to_datenum(gCtrl.pick.cursor)));

  % Update the map cursor
  set(hui.mapfig.pick_cursor.h,'XData', ...
    gCtrl.source.geo{gCtrl.source.cur_pick}.X(gCtrl.pick.cur_idx));
  set(hui.mapfig.pick_cursor.h,'YData', ...
    gCtrl.source.geo{gCtrl.source.cur_pick}.Y(gCtrl.pick.cur_idx));
  
  % Update the pick figure cursor
  if ishandle(hui.pickfig.handle)
    % If the pick window is open...
    cursor_xpos = interp1(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time, ...
      1:length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time), gCtrl.pick.cursor, 'linear','extrap');
    set(hui.pickfig.cursor.h,'XData',[cursor_xpos cursor_xpos]);
    set(hui.pickfig.cursor.h,'YData',[gCtrl.pick.time([1 end])*1e6]);
  end
  
  % ------------------------------------------------------------------
  % Find the closest point for the view frame and put marker on it
  [min_dist min_idx] = find_closest_point(x,y,gCtrl.source.cur_view);
  
  % Update the cursor data structures
  gCtrl.view.cur_idx = min_idx;
  gCtrl.view.cursor = gCtrl.source.geo{gCtrl.source.cur_view}.GPS_time(gCtrl.view.cur_idx);
  fprintf('  View: GPS %s\n', datestr(epoch_to_datenum(gCtrl.view.cursor)));

  % Update the map cursor
  set(hui.mapfig.view_cursor.h,'XData', ...
    gCtrl.source.geo{gCtrl.source.cur_view}.X(gCtrl.view.cur_idx));
  set(hui.mapfig.view_cursor.h,'YData', ...
    gCtrl.source.geo{gCtrl.source.cur_view}.Y(gCtrl.view.cur_idx));
  
  % Update the view figure cursor
  if ishandle(hui.viewfig.handle)
    % If the view window is open...
    cursor_xpos = interp1(gCtrl.source.geo{gCtrl.source.cur_view}.GPS_time, ...
      1:length(gCtrl.source.geo{gCtrl.source.cur_view}.GPS_time), gCtrl.view.cursor, 'linear','extrap');
    set(hui.viewfig.cursor.h,'XData',[cursor_xpos cursor_xpos]);
    set(hui.viewfig.cursor.h,'YData',[gCtrl.view.time([1 end])*1e6]);
  end

end

return;

% =====================================================================
% Support function for finding closest point overall or in specific frame
% =====================================================================
function [min_dist min_idx] = find_closest_point(x,y,frm)

global hui; % hui: user interface handles
global gCtrl; % Geobase: contains all the geographic info

if exist('frm','var')
  [min_dist min_idx] = min(sqrt((x - gCtrl.source.geo{frm}.X).^2 ...
    + (y - gCtrl.source.geo{frm}.Y).^2));
else
  [min_dist min_idx] = min(sqrt((x - gCtrl.source.X).^2 + (y - gCtrl.source.Y).^2));
end

return;
