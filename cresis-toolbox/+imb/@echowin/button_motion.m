function button_motion(obj,src,event)
% button_motion(obj,src, event)
%
% Called by echowin when mouse cursor is moved. Handles tracking mouse position
% and mouse cursor changing when outside of axes

set(obj.h_fig,'Units','normalized');
mouse_pos = get(obj.h_fig,'CurrentPoint');
set(obj.right_panel.handle,'Units','normalized');
uipanel_pos = get(obj.right_panel.handle,'Position');
set(obj.right_panel.handle,'Units','Points');
if obj.busy_mode
  set(obj.h_fig,'Pointer','watch');
else
  % check if inside status bar
  set(obj.right_panel.status_panel.handle,'Units','normalized');
  status_pos = get(obj.right_panel.status_panel.handle,'Position');
  status_h = status_pos(4);
  if mouse_pos(1) < uipanel_pos(1) || mouse_pos(2) < status_h
    set(obj.h_fig,'Pointer','Arrow');
    return;
  elseif obj.zoom_mode
    set(obj.h_fig,'Pointer','custom');
  end
  if ~obj.busy_mode
    % print current mouse position to status bar
    axis_pos = get(obj.right_panel.axes.handle,'CurrentPoint');
    axis_pos = axis_pos(1,:);
    xlim = sort(get(obj.right_panel.axes.handle,'XLim'));
    ylim = sort(get(obj.right_panel.axes.handle,'YLim'));
    below_x = axis_pos(1) < xlim(1);
    above_x = axis_pos(1) > xlim(2);
    below_y = axis_pos(2) < ylim(1);
    above_y = axis_pos(2) > ylim(2);
    % correct edges
    if below_x
      x = xlim(1);
    elseif above_x
      x = xlim(2);
    else
      x = axis_pos(1);
    end
    if below_y
      y = ylim(1);
    elseif above_y
      y = ylim(2);
    else
      y = axis_pos(2);
    end
    % get gps time
    gps_t = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,x,'linear','extrap');
    % get x and y positions
    try
      xpos = interp1(obj.eg.map_gps_time,obj.eg.map_x*1e3,gps_t,'linear','extrap');
      ypos = interp1(obj.eg.map_gps_time,obj.eg.map_y*1e3,gps_t,'linear','extrap');
    catch ME
      warning(sprintf('Duplicate values in map''s GPS time, should be addressed (%s)',obj.eg.cur_sel.day_seg));
      [unique_gps unique_idxs] = unique(obj.eg.map_gps_time);
      xpos = interp1(unique_gps,obj.eg.map_x(unique_idxs)*1e3,gps_t,'linear','extrap');
      ypos = interp1(unique_gps,obj.eg.map_y(unique_idxs)*1e3,gps_t,'linear','extrap');
    end
    loc = obj.eg.cur_sel.location;
    [lat,lon] = projinv(obj.eg.projmat,xpos,ypos);
    yaxis_unit = get(obj.left_panel.yaxisPM,'Value');
    switch yaxis_unit
      case 1
        y_unit = 'us';
      case 2
        y_unit = 'm';
      case 3
        y_unit = 'm';
      case 4
        y_unit = '';
    end
    set(obj.right_panel.status_panel.mouseCoordText,'String',sprintf('%8.3fN %8.3fW Y=%8.4f%s',lat,lon,y,y_unit));
  end
end

return

