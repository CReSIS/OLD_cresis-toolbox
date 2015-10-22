function plot_cursors(obj)
% echowin.plot_cursors(obj,h_obj,event)
%
% (Re)draws the cursor in the echogram window. Called by shift clicking in the
% echogram window as well as when shift clicking on the map positions the
% cursor within the limits of an echogram window currently opened.

% modify existing cursor
if isempty(obj.cursor.gps_time)
  set(obj.cursor.h,{'XData','YData'},{[],[]});
else
  cursor_xpos = interp1(obj.eg.image_gps_time,obj.eg.image_xaxis,obj.cursor.gps_time,'nearest','extrap');
  set(obj.cursor.h,{'XData','YData'},{[cursor_xpos cursor_xpos],...
    obj.eg.image_yaxis([1 end])});
  set(obj.cursor.h,'visible','on');
end

end
