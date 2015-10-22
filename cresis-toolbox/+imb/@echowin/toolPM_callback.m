function toolPM_callback(obj,hObj,event)

tool_idx = get(obj.left_panel.toolPM,'Value');
tmp = obj.tool_list{tool_idx}; obj.left_click = @tmp.left_click;
tmp = obj.tool_list{tool_idx}; obj.left_click_and_drag = @tmp.left_click_and_drag;
tmp = obj.tool_list{tool_idx}; obj.right_click_and_drag = @tmp.right_click_and_drag;

% don't run the following code if the new tool doesn't have a param
% interface (quality)
if ~isempty(obj.tool_list{get(obj.left_panel.toolPM,'Value')}.h_fig)
  if obj.tool_accessed
    % grab top left position of old tool window
    set(obj.tool_list{obj.old_tool_idx}.h_fig,'Units','Pixels');
    param_pos = get(obj.tool_list{obj.old_tool_idx}.h_fig,'Position');
    % get width and height of param window
    set(obj.tool_list{tool_idx}.h_fig,'Units','Pixels');
    param_win = get(obj.tool_list{tool_idx}.h_fig,'Position');
    param_w = param_win(3);
    param_h = param_win(4);
    if obj.tool_visible
      % close the old tool window if it's open
      % (the other case is when a tool was open, but then a tool with no
      % param figure was selected, such as the browse tool, and now a tool
      % with a param figure has been selcted)
      if strcmp(get(obj.tool_list{obj.old_tool_idx}.h_fig,'Visible'),'on')
        set(obj.tool_list{obj.old_tool_idx}.h_fig,'Visible','off');
      end
      % position the new tool window
      new_pos = [param_pos(1) param_pos(2)+param_pos(4)-param_h...
        param_w param_h];
      set(obj.tool_list{tool_idx}.h_fig,'Units','pixels');
      set(obj.tool_list{tool_idx}.h_fig,'Position',new_pos);
      % show the new tool window
      set(obj.tool_list{tool_idx}.h_fig,'Visible','on');
      figure(obj.h_fig);
    else
      % just position the new tool window
      new_pos = [param_pos(1) param_pos(2)+param_pos(4)-param_h...
        param_w param_h];
      set(obj.tool_list{tool_idx}.h_fig,'Units','pixels');
      set(obj.tool_list{tool_idx}.h_fig,'Position',new_pos);
    end
  end
  
  % save this tool idx so that when the tool is changed next time, the old
  % tool's window position can be grabbed
  obj.old_tool_idx = tool_idx;
else
  if strcmp(get(obj.tool_list{obj.old_tool_idx}.h_fig,'Visible'),'on')
    % just close the old tool window
    set(obj.tool_list{obj.old_tool_idx}.h_fig,'Visible','off');
    
  end
end

return
