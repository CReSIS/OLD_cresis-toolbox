function plot_cursors(obj)
% echowin.plot_cursors(obj)
%
% (Re)draws the cursor in the echogram window. Called by shift clicking in the
% echogram window as well as when shift clicking on the map positions the
% cursor within the limits of an echogram window currently opened.

% modify existing cursor
if isempty(obj.cursor.gps_time)
  set(obj.cursor.h,{'XData','YData'},{[],[]});
else
  [~,rline] = min(abs(obj.eg.image_gps_time - obj.cursor.gps_time));
  obj.update_cursor(obj.eg.image_xaxis(rline),NaN,true);
end

end
