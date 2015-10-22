function picker_view_button(src,event)
% pick_view_button(src,event)
%
% Support function for picker_view (the view figure window)
%
% Author: John Paden, Aric Beaver

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

[x,y,but] = get_mouse_info(hui.viewfig.handle,hui.viewfig.axes.handle);

x_time = interp1(1:length(gCtrl.source.geo{gCtrl.source.cur_view}.GPS_time), ...
  gCtrl.source.geo{gCtrl.source.cur_view}.GPS_time, x, 'linear','extrap');
fprintf('View: x = %.2f (GPS %s), y = %.2f us, but = %d\n', x, ...
  datestr(epoch_to_datenum(x_time),'yyyymmdd HH:MM:SS.FFF'), y, but);

if but == 3
  % ===================================================================
  % Right mouse button: move cursors

  % Find the gps time of click
  cursor_time = interp1(1:length(gCtrl.source.geo{gCtrl.source.cur_view}.GPS_time), ...
    gCtrl.source.geo{gCtrl.source.cur_view}.GPS_time, x, 'linear','extrap');
  
  % ------------------------------------------------------------------
  % Find the closest index
  [min_dist min_idx] = min(abs(cursor_time - gCtrl.source.geo{gCtrl.source.cur_view}.GPS_time));
  
  % Update the cursor data structures
  gCtrl.view.cur_idx = min_idx;
  gCtrl.view.cursor = gCtrl.source.geo{gCtrl.source.cur_view}.GPS_time(gCtrl.view.cur_idx);

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
  
  physical_constants;
  rline_idx = min_idx;
  surface_time = gCtrl.source.layers{gCtrl.source.cur_view}.layerData{1}.value{2}.data(rline_idx);
  bottom_time = gCtrl.source.layers{gCtrl.source.cur_view}.layerData{2}.value{2}.data(rline_idx);
  thickness_time = bottom_time-surface_time;
  if ~isfinite(surface_time) && ~isfinite(bottom_time)
    % Print nothing    
  elseif isfinite(surface_time) && isfinite(bottom_time)
    fprintf('  Surface: %.2fus, Bottom: %.2fus, Thickness: %.2fus (%.2fm)\n',...
      surface_time*1e6,bottom_time*1e6,thickness_time*1e6,thickness_time*c/2/sqrt(er_ice));
  elseif isfinite(surface_time) && ~isfinite(bottom_time)
    fprintf('  Surface: %.2fus, Bottom: n/a, Thickness: n/a\n',surface_time*1e6);
  elseif ~isfinite(surface_time) && isfinite(bottom_time)
    fprintf('  Surface: n/a, Bottom: %.2fus, Thickness: n/a\n',bottom_time*1e6);
  else
    % Do nothing?
  end
end

return
