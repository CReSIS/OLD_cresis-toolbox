function paramPB_callback(obj,hObj,event)

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);

% Clicking this button always make the tool parameters figure visible
if ~isempty(obj.tool.list{get(obj.left_panel.toolPM,'Value')}.h_fig)
  if ~obj.tool.accessed
    this_pos = get(obj.h_fig,'Position'); %echowin position
    param_pos = get(obj.tool.list{get(obj.left_panel.toolPM,'Value')}.h_fig,'Position');
    new_pos = [this_pos(1) this_pos(2)+this_pos(4)-param_pos(4) param_pos(3) param_pos(4)];
    set(obj.tool.list{get(obj.left_panel.toolPM,'Value')}.h_fig,'Position',new_pos);
    set(obj.tool.list{get(obj.left_panel.toolPM,'Value')}.h_fig,'Visible','on');
    obj.tool.accessed = true;
    obj.tool.visible = true;
  else
    set(obj.tool.list{get(obj.left_panel.toolPM,'Value')}.h_fig,'Visible','on');
    figure(obj.tool.list{get(obj.left_panel.toolPM,'Value')}.h_fig);
    obj.tool.visible = true;
  end
end
