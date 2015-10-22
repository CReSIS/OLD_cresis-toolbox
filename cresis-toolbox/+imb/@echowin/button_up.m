function button_up(obj,src,event)
% echowin.button_up(obj,src, event)
%
% Called by echowin when a mouse press is released. Handles zoom in/out,
% application of tools, and marker updating.

% Make sure that click is on the right side panel
mouse_pos = get(obj.h_fig,'CurrentPoint');

% Check to make sure mouse clicked inside of obj.right_panel.handle
%   Since extends the full y-length, just check to the right of minimum x
set(obj.right_panel.handle,'Units','normalized');
uipanel_pos = get(obj.right_panel.handle,'Position');
set(obj.right_panel.handle,'Units','Points');
if mouse_pos(1) <= uipanel_pos(1)
  return
end

% Check to make sure mouse clicked inside of obj.right_panel.axes.panel_h
%   Since extends the full x-length of obj.right_panel.handle, just check 
%   up of minimum y
set(obj.right_panel.axes.panel_h,'Units','normalized');
uipanel_pos = get(obj.right_panel.axes.panel_h,'Position');
set(obj.right_panel.axes.panel_h,'Units','Points');
if mouse_pos(2) <= uipanel_pos(2)
  return
end

[x,y,but] = get_mouse_info(obj.h_fig,obj.right_panel.axes.handle);
%fprintf('Echo Button Up: x = %.3f, y = %.3f, but = %d\n', x, y, but);

if isempty(obj.click_x)
  obj.click_x = x;
  obj.click_y = y;
end

cur_axis = axis(obj.right_panel.axes.handle);
if y < cur_axis(3) || y > cur_axis(4)
  clicked_outside_axis = true;
else
  clicked_outside_axis = false;
end

% ===================================================================
% Tool Mode
% Left click: Apply tool to a point
% Alt + Left click and drag: Apply tool to a region
% Right click and drag: Delete region
% Zoom Mode
% Left click: Zoom at point
% Left click and drag: Zoom to region
% Right click: Zoom out at point
% Double click: Zoom reset
% Any Mode
% Scroll: Zooms in/out (handled by button_scroll.m)
% Ctrl + any click: Select layer
% Shift + any click: Set marker point
if obj.control_pressed
  %% Ctrl + any click: Select closest layer and crossover

  % Select Layer
  distances = zeros(1,length(obj.layer_h)/2);
  for idx = 1:length(distances)
    layer_data_x = cat(2,get(obj.layer_h(2*idx-1),'Xdata'),get(obj.layer_h(2*idx),'Xdata'));
    layer_data_y = cat(2,get(obj.layer_h(2*idx-1),'Ydata'),get(obj.layer_h(2*idx),'Ydata'));
    if ~isempty(layer_data_x)
      [tmp pnt_idx] = min(abs(layer_data_x-x));
      distances(idx) = sqrt((layer_data_y(pnt_idx)-y)^2+(layer_data_x(pnt_idx)-x)^2);
    else
      distances(idx) = inf;
    end
  end
  [tmp min_idx] = min(distances);
  obj.layerLB_sync('sel',min_idx);
  
  % Select crossover
  distances = abs(obj.eg.crossovers.x_curUnit - x).^2 + abs(obj.eg.crossovers.y_curUnit - y).^2;
  [tmp min_idx] = min(distances);
  obj.eg.crossovers.gui.set_selected(min_idx);
  
  % Update the visibility of selected layers/crossovers
  obj.set_visibility();
  return;
elseif obj.shift_pressed
  %% Shift + any click: Set cursor
  % get absolute cursor position
  obj.cursor.gps_time = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,...
    x,'linear','extrap');
  obj.cursor.x = interp1(obj.eg.map_gps_time,obj.eg.map_x,obj.cursor.gps_time,'linear','extrap');
  obj.cursor.y = interp1(obj.eg.map_gps_time,obj.eg.map_y,obj.cursor.gps_time,'linear','extrap');
  
  % (re)draw cursor in this echogram window
  obj.plot_cursors();

  str = obj.status_text_cursor();
  obj.status_text_set(str,'replace');
  
  % tell the map to redraw all cursors
  notify(obj,'update_cursors');
  return;
end

if obj.zoom_mode
  if but == 4
    %% Double click: Zoom reset
    obj.redraw(-inf,inf,-inf,inf,struct('clipped',true));
    
  elseif but == 1 && x~=obj.click_x && y~=obj.click_y
    %% Left click and drag: Zoom to region
    
    xaxis = get(obj.right_panel.axes.handle,'XLim');
    yaxis = get(obj.right_panel.axes.handle,'YLim');
    [x_min x_max y_min y_max] = imb.sort_clicks(xaxis,yaxis,obj.click_x,obj.click_y,x,y);
    
    % check if selection is valid(sort_clicks returns weird values if init/
    % final clicks are both outside the same axis boundary)
    if x_min > x_max || x_max < x_min || y_min > y_max || y_max < y_min
      return;
    end   
    
    % Convert x_min, x_max to GPS time
    xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,[x_min x_max],'linear','extrap');
    
    % Redraw echogram imagesc with new axis
    obj.redraw(xlims(1),xlims(2),y_min,y_max);
    
  elseif but == 1
    %% Left click: Zoom at point
    if clicked_outside_axis
      return;
    end
    zooms = -1.5;
    cur_axis = axis(obj.right_panel.axes.handle);
    y_extent = cur_axis(4) - cur_axis(3);
    x_extent = cur_axis(2) - cur_axis(1);
    xlims = [x - x_extent*2^zooms, x + x_extent*2^zooms];
    ylims = [y - y_extent*2^zooms, y + y_extent*2^zooms];
    
    % Convert x_min, x_max to GPS time
    xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,xlims,'linear','extrap');
    
    % Draw data with new axis, but do not allow new data to be loaded (i.e.
    % clip new axis to limits of loaded data
    obj.redraw(xlims(1),xlims(2),ylims(1),ylims(2),struct('clipped',true));
    
  elseif but == 3
    %% Right click: Zoom out at point
    zooms = -0.5;
    cur_axis = axis(obj.right_panel.axes.handle);
    y_extent = cur_axis(4) - cur_axis(3);
    x_extent = cur_axis(2) - cur_axis(1);
    xlims = [x - x_extent*2^zooms, x + x_extent*2^zooms];
    ylims = [y - y_extent*2^zooms, y + y_extent*2^zooms];
    
    % Convert x_min, x_max to GPS time
    xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,xlims,'linear','extrap');
    
    % Draw data with new axis, but do not allow new data to be loaded (i.e.
    % clip new axis to limits of loaded data
    obj.redraw(xlims(1),xlims(2),ylims(1),ylims(2),struct('clipped',true));
  end
else
  %% Populate param structure
  % ====================================================================
  
  % Current quality
  param.cur_quality = get(obj.left_panel.qualityPM,'Value');
  % Current layers
  param.cur_layers = find(obj.left_panel.layer_panel.selected_layers).';
  param.layer.x = obj.eg.layer.x_curUnit{1};
  param.layer.y = obj.eg.layer.y_curUnit;
  param.layer.type = obj.eg.layer.type;
  param.layer.qual = obj.eg.layer.qual;
  % Echogram and layer information
  param.image_x = get(obj.left_panel.imagewin.img,'XData');
  param.image_y = get(obj.left_panel.imagewin.img,'YData');
  param.image_c = get(obj.left_panel.imagewin.img,'CData');
  
  %% Find proper param window position (for browse tool)
  if ~obj.tool_accessed
    param.keep_tool_pos = false;
    set(obj.h_fig,'Units','Pixels');
    this_pos = get(obj.h_fig,'Position');
    param.tool_x = this_pos(1);
    param.tool_y = this_pos(2)+this_pos(4); % still need to subtract the tool window's height for proper y pos
    % if browse tool is currently selected, set tool_accessed to be true
    tool_idx = get(obj.left_panel.toolPM,'Value');
    if tool_idx == 4
      obj.tool_accessed = true;
    end
  else
    param.keep_tool_pos = true;
  end
  
  %% Apply command based on mouse/keyboard click type
  if but == 1 % left click
    if obj.alt_pressed && x~=obj.click_x && y~=obj.click_y
      %% Alt + Left click and drag: Apply tool to a region
      param.x = sort([obj.click_x x]);
      param.y = sort([obj.click_y y]);
      cmds = obj.left_click_and_drag(param);
    else
      %% Left click: Apply tool to a point
      if clicked_outside_axis
        return;
      end
      param.x = x;
      param.y = y;
      cmds = obj.left_click(param);
    end
    if ~isempty(cmds)
      % Push the new command(s) to the stack
      obj.undo_stack.push(obj.cmds_convert_units(cmds));
    end
  elseif but==3
    if x~=obj.click_x && y~=obj.click_y
      %% Right click and drag: Apply tool to a region
      set(obj.right_panel.echoCM,'visible','off');
      param.x = sort([obj.click_x x]);
      param.y = sort([obj.click_y y]);
      param.point_path_id = obj.eg.map_id;
      cmds = obj.right_click_and_drag(param);
      if ~isempty(cmds)
        % Push the new command(s) to the stack
        obj.undo_stack.push(obj.cmds_convert_units(cmds));
      end
    end
  end
end

return

