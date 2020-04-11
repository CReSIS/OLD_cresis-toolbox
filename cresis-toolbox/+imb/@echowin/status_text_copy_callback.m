function status_text_copy(obj,source,event)
% echowin.status_text_copy(obj,source,event)
%
% Copies the status text to the clipboard

% Used to copy status text:
%str = get(obj.right_panel.status_panel.statusText,'String');

if ~isempty(obj.cursor.gps_time)
  % Determine which frame the cursor is in
  cur_frm = find(obj.cursor.gps_time >= obj.eg.start_gps_time,1,'last');
  if isempty(cur_frm)
    cur_frm = 1;
  end
  
  status_str = status_text_cursor(obj);
  
  str = sprintf('%s_%03d: %s', obj.eg.cur_sel.day_seg, cur_frm, status_str);
  
  clipboard('copy',str);
end
