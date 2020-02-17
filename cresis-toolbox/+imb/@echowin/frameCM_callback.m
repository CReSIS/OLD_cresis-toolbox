function frameCM_callback(obj,hObj,event)
% echowin.frameCM_callback(obj,hObj,event)
%
% Copies to clipboard the currently selected frame.

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);

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
