function paramPB_callback(obj,hObj,event)

% Clicking this button always make the tool parameters figure visible
if ~isempty(obj.tool_list{get(obj.left_panel.toolPM,'Value')}.h_fig)
  if ~obj.tool_accessed
    set(obj.tool_list{get(obj.left_panel.toolPM,'Value')}.h_fig,'Units','pixels');
    set(obj.h_fig,'Units','Pixels');
    this_pos = get(obj.h_fig,'Position'); %echowin position
    param_pos = get(obj.tool_list{get(obj.left_panel.toolPM,'Value')}.h_fig,'Position');
    new_pos = [this_pos(1) this_pos(2)+this_pos(4)-param_pos(4) param_pos(3) param_pos(4)];
    set(obj.tool_list{get(obj.left_panel.toolPM,'Value')}.h_fig,'Position',new_pos);
    set(obj.tool_list{get(obj.left_panel.toolPM,'Value')}.h_fig,'Visible','on');
    obj.tool_accessed = true;
    obj.tool_visible = true;
  else
    set(obj.tool_list{get(obj.left_panel.toolPM,'Value')}.h_fig,'Visible','on');
    obj.tool_visible = true;
  end
  % Give focus to the tool param window
  try
    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    javaFrame = get(obj.tool_list{get(obj.left_panel.toolPM,'Value')}.h_fig,'JavaFrame');
    javaFrame.getAxisComponent.requestFocus;
  catch
    obj.status_text_set(sprintf('Focus error, click inside echogram window before using key shortcuts'),'replace');
  end
end

return

