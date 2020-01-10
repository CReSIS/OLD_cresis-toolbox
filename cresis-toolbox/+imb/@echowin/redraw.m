function redraw(obj,x_min,x_max,y_min,y_max,param)
% redraw(obj,x_min,x_max,y_min,y_max,param)
%
% Member function of imb.echowin class.
%
% Replots the echogram because of a change in the x-extent from left/right
% arrow keys, zoom in/out, or from a selection from the frame list box.

% need: image_xaxis_gpstime[1,end]

if ~exist('param','var')
  param = struct();
end
if ~isfield(param,'clipped')
  % param.clipped = 0
  %   Will load data only when the image limits were already at the edge
  %   of the loaded data
  % param.clipped = 1 [default]
  %   Will not load new data even if request is outside bounds of loaded
  %   data
  % param.clipped = 2
  %   Will always load new data if needed
  % param.clipped = 3
  %   Force a complete reload of the data
  param.clipped = 1;
end

%% Check x_min and x_max against segment boundaries
if x_min == -inf
  x_min = obj.eg.start_gps_time(1);
end
if x_max == inf
  x_max = obj.eg.stop_gps_time(end);
end
% if x_min < obj.eg.start_gps_time(1)
%   x_range = x_max - x_min;
%   x_min = obj.eg.start_gps_time(1);
%   x_max = x_min + x_range;
% end
% if x_max > obj.eg.stop_gps_time(end)
%   x_range = x_max - x_min;
%   x_max = obj.eg.stop_gps_time(end);
%   x_min = max(obj.eg.start_gps_time(1), x_max - x_range);
% end
if y_min < min(obj.eg.image_yaxis([1 end]))
  y_range = y_max - y_min;
  y_min = min(obj.eg.image_yaxis([1 end]));
  y_max = y_min + y_range;
end
if y_max > max(obj.eg.image_yaxis([1 end]))
  y_range = y_max - y_min;
  y_max = max(obj.eg.image_yaxis([1 end]));
  y_min = max(min(obj.eg.image_yaxis([1 end])), y_max - y_range);
end

%% Check to see if the new x_extent is outside the currently loaded data
%   To avoid unnecessarily loading data, just move the x-limits to
%   the edge of the currently loaded frames UNLESS the x-limits were
%   already on the edge, in which case we'll load the new data.
dx = obj.eg.gps_time(2)-obj.eg.gps_time(1);
min_gps_time = obj.eg.gps_time(1);
max_gps_time = obj.eg.gps_time(end);
if obj.eg.frame_idxs(1) == 1
  % First frame in segment so adjust beginning to make sure all layer points
  % will be displayed.
  min_gps_time = obj.eg.start_gps_time(1)-dx;
elseif obj.eg.frame_idxs(end) == length(obj.eg.stop_gps_time)
  % Last frame of segment so adjust end to make sure all layer points will 
  % be displayed.
  max_gps_time = obj.eg.stop_gps_time(end)+dx;
end
cur_xlims = xlim(obj.right_panel.axes.handle);
cur_xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,cur_xlims,'linear','extrap');
load_new_data = 0;
desire_frame_idxs = obj.eg.frame_idxs;
if param.clipped == 2
  frame_idx = find(obj.eg.start_gps_time*(1-10*eps) <= x_min, 1,'last'); % Catch rounding errors with 10*eps
  load_new_data = frame_idx ~= obj.eg.frame_idxs(1);
  if diff([obj.eg.frame_idxs(1),frame_idx]) < 0
    load_new_data = load_new_data * -1;
  end
  % Try to use as much of the currently existing loaded data as possible
  %     if load_new_data > 0 && load_new_data < 2*obj.default_params.max_frames-1
  %       load_new_data = max(0,load_new_data - obj.default_params.max_frames + 1);
  %     end
  if load_new_data
    desire_frame_idxs = frame_idx + (0:obj.default_params.max_frames-1);
    desire_frame_idxs = desire_frame_idxs(desire_frame_idxs >= 1 ...
      & desire_frame_idxs <= length(obj.eg.frame_names));
  end
elseif param.clipped == 3
  load_new_data = 1;
  desire_frame_idxs = obj.eg.frame_idxs(1) + (0:obj.default_params.max_frames-1);
  desire_frame_idxs = desire_frame_idxs(desire_frame_idxs >= 1 ...
    & desire_frame_idxs <= length(obj.eg.frame_names));
else
  if x_min < min_gps_time
    % Request is for new data before currently loaded frames
    if param.clipped == 1 || param.clipped == 0 && cur_xlims(1) > min_gps_time
      % Just move the x-limits to the edge of the currently loaded frames
      x_range = x_max - x_min;
      x_min = min_gps_time;
      x_max = x_min + x_range;
    else
      % Load previous data
      load_new_data = -1;
    end
  elseif ~load_new_data && x_max > max_gps_time
    % Request is for new data after currently loaded frames
    if param.clipped == 1 || param.clipped == 0 && cur_xlims(2) < max_gps_time
      % Just move the x-limits to the edge of the currently loaded frames
      x_range = x_max - x_min;
      x_max = max_gps_time;
      x_min = max(min_gps_time,x_max - x_range);
    else
      % Load previous data
      load_new_data = 1;
    end
  end
  if load_new_data
    if (load_new_data > 0 && obj.eg.frame_idxs(end) == length(obj.eg.frame_names))...
        || (load_new_data < 0 && obj.eg.frame_idxs(1) == 1)
      load_new_data = false;
    end
    if load_new_data > 0
      desire_frame_idxs = obj.eg.frame_idxs(end) - 1 + load_new_data + (0:obj.default_params.max_frames-1);
      desire_frame_idxs = desire_frame_idxs(desire_frame_idxs >= 1 ...
        & desire_frame_idxs <= length(obj.eg.frame_names));
    elseif load_new_data < 0
      desire_frame_idxs = sort(obj.eg.frame_idxs(1)+ 1 + load_new_data + sign(load_new_data).*(0:obj.default_params.max_frames-1));
      desire_frame_idxs = desire_frame_idxs(desire_frame_idxs >= 1 ...
        & desire_frame_idxs <= length(obj.eg.frame_names));
    end
  end
end

%% Update the current frame selection in frameLB (load_echogram uses this)
% First determine which frame is on the left most side of the box
cur_frm = find(x_min >= obj.eg.start_gps_time,1,'last');
if isempty(cur_frm)
  cur_frm = 1;
end
  
% Update the selection in the frame LB and the source LB
obj.update_frame_and_sourceLB(cur_frm);

if ~load_new_data
  %% Just updating the axis and not loading data
  
  % Make sure axis call does not stretch beyond image limits
  if x_min < min_gps_time
    x_range = x_max - x_min;
    x_min = min_gps_time;
    x_max = x_min + x_range;
  end
  if x_max > max_gps_time
    x_range = x_max - x_min;
    x_max = max_gps_time;
    x_min = max(min_gps_time, x_max - x_range);
  end
  
  % Convert from GPS time to x-axis
  xlims = interp1(obj.eg.image_gps_time,obj.eg.image_xaxis,[x_min x_max],'linear','extrap');
  axis(obj.right_panel.axes.handle,[xlims y_min y_max]);
  
  % Change dynamic range
  obj.left_panel.imagewin.update_caxis();
  
else
  obj.busy_mode = true;
  set(obj.h_fig,'Pointer','watch');
  obj.status_text_set(sprintf('(%s) Redrawing...', datestr(now,'HH:MM:SS')),'replace');
  drawnow;
  fprintf('START ECHOWIN REDRAW (%s)\n',datestr(now,'HH:MM:SS'));
  
  load_new_data = 0;
  clipped = param.clipped;
  
  old_frame_idxs = obj.eg.frame_idxs;
  [x_min,x_max,y_min,y_max] = obj.load_echogram(desire_frame_idxs,clipped,x_min,x_max,y_min,y_max);
  
  if length(old_frame_idxs) == length(desire_frame_idxs) ...
      && all(old_frame_idxs == desire_frame_idxs)
    % No frames changed
    % Update echogram surface if there are enough good points from OPS
    % Find good surface points
    good_mask = ~isnan(obj.eg.layer.y{obj.eg.layer_id==1});
    if sum(good_mask) > 2
      % There are surface layer points in the database, overwrite the surface
      % with these
      obj.eg.surface = interp1(obj.eg.map_gps_time(good_mask),obj.eg.layer.y{obj.eg.layer_id==1}(good_mask),obj.eg.gps_time);
      obj.eg.surface = interp_finite(obj.eg.surface,0);
    end
    obj.plot_echogram(x_min,x_max,y_min,y_max);
    obj.plot_layers();
    obj.plot_cursors();
    obj.plot_crossovers();
  else
    obj.load_flightline();
    obj.load_layers();
    obj.plot_echogram(x_min,x_max,y_min,y_max);
    obj.plot_layers();
    obj.crossovers_en = obj.eg.crossovers.gui.crossovers_en();
    obj.load_crossovers();
    obj.plot_cursors();
  end
  
  % Since we have reloaded the layers, we must resync to our undo stack
  % to re-execute any un-saved changes
  cmds_list = obj.undo_stack.get_save_cmds(false);
  obj.cmds_execute(cmds_list,'redo');
  
  %% Set layer and cross over visibility
  obj.set_visibility();
  
  %% Plot cursor (if valid)
  xlimits = xlim(obj.right_panel.axes.handle);
  start_gps = obj.eg.image_gps_time(find(obj.eg.image_xaxis>=xlimits(1),1));
  stop_gps = obj.eg.image_gps_time(find(obj.eg.image_xaxis<=xlimits(2),1,'last'));
  
  if ~isempty(obj.cursor.gps_time)
    obj.plot_cursors();
  end
  
  %% Cleanup
  fprintf('DONE ECHOWIN REDRAW (%s)\n',datestr(now,'HH:MM:SS'));
  
  obj.status_text_set(sprintf(' done. (%s)', datestr(now,'HH:MM:SS')),'append');
  obj.busy_mode = false;
  if obj.zoom_mode
    set(obj.h_fig,'Pointer','custom');
  else
    set(obj.h_fig,'Pointer','Arrow');
  end
  
  fprintf(' DONE (%s)\n',datestr(now,'HH:MM:SS'));
end
  
% Plot new selection on flight path (calls
% mapwin.update_echowin_flightlines)
notify(obj,'update_echowin_flightline');

return;
