function toolPM_callback(obj,hObj,event)
% toolPM_callback(obj,hObj,event)
%
% Switch tools based on user selection in tool popup menu

tool_idx = get(obj.left_panel.toolPM,'Value');
tmp = obj.tool_list{tool_idx}; obj.left_click = @tmp.left_click;
tmp = obj.tool_list{tool_idx}; obj.left_click_and_drag = @tmp.left_click_and_drag;
tmp = obj.tool_list{tool_idx}; obj.right_click_and_drag = @tmp.right_click_and_drag;

% Hide the old tool window
set(obj.tool_list{obj.old_tool_idx}.h_fig,'Visible','off');

% Check to make sure tool has a parameter window
if ~isempty(obj.tool_list{get(obj.left_panel.toolPM,'Value')}.h_fig)
  % If any tool's parameter window has been accessed, then show the current
  % tool's parameter window in the same spot. If no parameter window has
  % been opened there is nothing to be done yet.
  if obj.tool_accessed
    % Grab top left position of old tool window
    set(obj.tool_list{obj.old_tool_idx}.h_fig,'Units','Pixels');
    param_pos = get(obj.tool_list{obj.old_tool_idx}.h_fig,'Position');
    % Get width and height of this tool's param window
    set(obj.tool_list{tool_idx}.h_fig,'Units','Pixels');
    param_win = get(obj.tool_list{tool_idx}.h_fig,'Position');
    param_w = param_win(3);
    param_h = param_win(4);
    % Update the new tool's parameter window location
    new_pos = [param_pos(1) param_pos(2)+param_pos(4)-param_h...
      param_w param_h];
    set(obj.tool_list{tool_idx}.h_fig,'Units','pixels');
    set(obj.tool_list{tool_idx}.h_fig,'Position',new_pos);
    if obj.tool_visible
      % Show the new tool's parameter window
      set(obj.tool_list{tool_idx}.h_fig,'Visible','on');
      %figure(obj.h_fig);
    end
  end
  
  % Save the current tool's index so that the next tool selected with a
  % parameter window can open its parameter window in the same spot as this
  % tool's window was.
  obj.old_tool_idx = tool_idx;
end
