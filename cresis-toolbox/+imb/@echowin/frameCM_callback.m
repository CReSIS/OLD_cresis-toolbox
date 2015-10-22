function frameCM_callback(obj,hObj,event)
% echowin.frameCM_callback(obj,hObj,event)
%
% Copies to clipboard the currently selected frame.

cur_sel = get(obj.left_panel.frameLB,'Value');
frame_list = get(obj.left_panel.frameLB,'String');
try
  str = frame_list{cur_sel};
catch
  str = '';
end
if ~isempty(str)
  clipboard('copy',str);
end

end
