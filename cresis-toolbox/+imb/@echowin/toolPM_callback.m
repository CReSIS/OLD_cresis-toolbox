function toolPM_callback(obj,hObj,event)
% toolPM_callback(obj,hObj,event)
%
% Switch tools based on user selection in tool popup menu

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);

tool_idx = get(obj.left_panel.toolPM,'Value');
tmp = obj.tool.list{tool_idx}; obj.tool.left_click_fh = @tmp.left_click;
tmp = obj.tool.list{tool_idx}; obj.tool.left_click_and_drag_fh = @tmp.left_click_and_drag;
tmp = obj.tool.list{tool_idx}; obj.tool.right_click_fh = @tmp.right_click;
tmp = obj.tool.list{tool_idx}; obj.tool.right_click_and_drag_fh = @tmp.right_click_and_drag;

% Hide the old tool window
set(obj.tool.list{obj.tool.old_idx}.h_fig,'Visible','off');

% Check to make sure tool has a parameter window
if ~isempty(obj.tool.list{get(obj.left_panel.toolPM,'Value')}.h_fig)
  % If any tool's parameter window has been accessed, then show the current
  % tool's parameter window in the same spot. If no parameter window has
  % been opened there is nothing to be done yet.
  if obj.tool.accessed
    % Grab top left position of old tool window
    param_pos = get(obj.tool.list{obj.tool.old_idx}.h_fig,'Position');
    % Get width and height of this tool's param window
    param_win = get(obj.tool.list{tool_idx}.h_fig,'Position');
    param_w = param_win(3);
    param_h = param_win(4);
    % Update the new tool's parameter window location
    new_pos = [param_pos(1) param_pos(2)+param_pos(4)-param_h...
      param_w param_h];
    set(obj.tool.list{tool_idx}.h_fig,'Position',new_pos);
    if obj.tool.visible
      % Show the new tool's parameter window
      set(obj.tool.list{tool_idx}.h_fig,'Visible','on');
      figure(obj.h_fig);
    end
  end
  
  % Save the current tool's index so that the next tool selected with a
  % parameter window can open its parameter window in the same spot as this
  % tool's window was.
  obj.tool.old_idx = tool_idx;
end
