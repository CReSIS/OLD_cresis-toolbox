function convert_data_to_image(obj,x_min,x_max,y_min,y_max)
% convert_data_to_image(obj,x_min,x_max,y_min,y_max)
% 
% Member function of imb.echowin class.
% 
% Plots data in obj.eg native coordinates (twtt and range-lines) onto
% the echogram window which may have a different coordinate system selected:
%  1. Echogram
%  2. Layer data
%  3. Cursor
%  4. Cross overs
% Called from draw, redraw, x-axis, and y-axis popup menu call backs.

physical_constants;
% ======================================================================
%% Convert the data along the x-axis according to the units
xaxis_choice = get(obj.left_panel.xaxisPM,'Value');
if xaxis_choice == 1 % rangeline
  % update image_xaxis and image_gps_time
  obj.eg.image_xaxis = 1:length(obj.eg.gps_time);
  obj.eg.image_gps_time = obj.eg.gps_time;
  % convert crossover gps times
  obj.eg.crossovers.x_curUnit = interp1(obj.eg.image_gps_time,1:length(obj.eg.image_gps_time),obj.eg.crossovers.gps_time);
  % update image_data according to xaxis_gpstime
  obj.eg.image_data = obj.eg.data;
  % update x label
  obj.eg.x_label = 'Range Line';
elseif xaxis_choice == 2 % Along track
  along_track = geodetic_to_along_track(obj.eg.latitude, obj.eg.longitude, obj.eg.elevation);
  % Set along track sampling, estimate best along-track pixel size, dx, from the data
  dx = median(diff(along_track));
  along_track_uniform = along_track(1): dx :along_track(end);
  % update image_xaxis and image_gps_time
  obj.eg.image_xaxis = along_track_uniform/1000;
  obj.eg.image_gps_time = interp1(along_track,obj.eg.gps_time,...
    along_track_uniform,'linear');
  % convert crossover gps times
  obj.eg.crossovers.x_curUnit = interp1(along_track_uniform,obj.eg.image_gps_time,...
    obj.eg.crossovers.gps_time,'linear');
  % update image_data according to image_gps_time
  obj.eg.image_data = interp1(obj.eg.gps_time,...
    obj.eg.data.',obj.eg.image_gps_time,'linear').';
  % update x label
  obj.eg.x_label = 'Along track (km)';
elseif xaxis_choice == 3 % GPS time
  tmp_gpstime = obj.eg.gps_time;
  % Set gps time sampling, estimate best gps time pixel size, dg, from the data
  dg = median(diff(tmp_gpstime));
  gps_time_uniform = tmp_gpstime(1):dg:tmp_gpstime(end);
  % update image_xaxis and image_gps_time
  obj.eg.image_gps_time = gps_time_uniform;
  obj.eg.image_xaxis = gps_time_uniform;
  % convert crossover gps times
  obj.eg.crossovers.x_curUnit = obj.eg.crossovers.gps_time;
  % update display_data according to xaxis_gpstime
  obj.eg.image_data = interp1(obj.eg.gps_time,...
    obj.eg.data.',obj.eg.image_gps_time,'linear').';
  % update x label
  obj.eg.x_label = 'GPS time';
end

% ======================================================================
%% Convert the data along the y-axis according to the units
% perform y-axis conversion (from twtt)
yaxis_choice = get(obj.left_panel.yaxisPM,'Value');
if yaxis_choice == 1 % TWTT
  % update yaxis and yaxis_time
  obj.eg.image_yaxis = obj.eg.time*1e6;
  % update y label
  obj.eg.y_label = 'Two-way Propagation (us)';
  obj.eg.y_order = 'reverse';
  % convert layerPnts_y from TWTT to current unit
  layer_y_curUnit = cell(1,length(obj.eg.layer.y));
  for idx = 1:length(obj.eg.layer.y)
    layer_y_curUnit{idx} = obj.eg.layer.y{idx}*1e6;
  end
  % convert crossover values too
  % calculate elevation offsets for crossovers first
  % raw elevation_1/elevation_2 in meters, elev_offset in twtt
  obj.eg.crossovers.elev_offset = 2*(obj.eg.crossovers.cross_elev-obj.eg.crossovers.source_elev)./c.*1e6;
  % calculate crossover twtts with elevation differences compensated
  obj.eg.crossovers.y_curUnit = obj.eg.crossovers.twtt*1e6 + obj.eg.crossovers.elev_offset;
elseif yaxis_choice == 2 % WGS_84 Elevation
  elevation = interp1(obj.eg.gps_time,...
    obj.eg.elevation,obj.eg.image_gps_time,'linear');
  time = obj.eg.time;
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surface,obj.eg.image_gps_time,'linear');
  physical_constants;
  elev_max = max(elevation - time(1)*c/2);
  elev_min = min(elevation - surface*c/2 - (time(end)-surface)*c/(sqrt(er_ice)*2));
  dt = time(2) - time(1);
  drange = dt * c/(sqrt(er_ice)*2);
  elev_uniform = (elev_max:-drange:elev_min).';
  % update image_data
  Nt = size(obj.eg.image_data,1);
  obj.eg.image_data = [obj.eg.image_data;...
    zeros(length(elev_uniform)-Nt,size(obj.eg.image_data,2))];
  warning('off','MATLAB:interp1:NaNinY')
  for idx = 1:length(surface)
    range = min(time,surface(idx)) * c/2 ...
      + max(0,time-surface(idx)) * c/(sqrt(er_ice)*2);
    elev = elevation(idx) - range;
    obj.eg.image_data(:,idx) = interp1(elev,...
      obj.eg.image_data(1:Nt,idx),elev_uniform,'linear');
  end
  warning('on','MATLAB:interp1:NaNinY')
  % update image_yaxis
  obj.eg.image_yaxis = elev_uniform;
  % update y label
  obj.eg.y_label = 'Elevation (m)';
  obj.eg.y_order = 'normal';
  % Convert layerPnts_y from TWTT to WGS_84 Elevation
  layer_y_curUnit = obj.eg.layer.y;
  for idx = 1: length(layer_y_curUnit)
    elevation = interp1(obj.eg.gps_time,...
      obj.eg.elevation,...
      obj.eg.layer.x{idx},'linear');
    surface = interp1(obj.eg.gps_time,...
      obj.eg.surface,...
      obj.eg.layer.x{idx},'linear');
    for pnt_idx = 1:length(layer_y_curUnit{idx})
      range = min(layer_y_curUnit{idx}(pnt_idx),surface(pnt_idx))*c/2 ...
        +max(0,layer_y_curUnit{idx}(pnt_idx)-surface(pnt_idx)) * c/(sqrt(er_ice)*2);
      layer_y_curUnit{idx}(pnt_idx) = elevation(pnt_idx) - range;
    end
  end
  % convert crossovers too
  elevation = obj.eg.crossovers.source_elev;
  obj.eg.crossovers.elev_offset = 2*(obj.eg.crossovers.cross_elev-obj.eg.crossovers.source_elev)./c;
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surface,...
    obj.eg.crossovers.gps_time,'linear');
  % correct surface using elevation offset info (since we are using surface
  % for frame1 and not frame2 and the twtt of the crossover is from frame2)
  surface = surface - obj.eg.crossovers.elev_offset;
  range = min(obj.eg.crossovers.twtt,surface).*c/2 ...
        + max(0,obj.eg.crossovers.twtt-surface).*c/(sqrt(er_ice)*2);
  obj.eg.crossovers.y_curUnit = elevation - range;
elseif yaxis_choice == 3 % Depth/Range
  time = obj.eg.time;
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surface,obj.eg.image_gps_time,'linear');
  range_min = time(1)*c/2;
  range_max = max(surface*c/2 + (time(end)-surface)*c/(sqrt(er_ice)*2));
  dt = time(2) - time(1);
  drange = dt * c/(sqrt(er_ice)*2);
  range_uniform = (range_min:drange:range_max).';
  % update image_data
  Nt = size(obj.eg.image_data,1);
  obj.eg.image_data = [obj.eg.image_data;...
    zeros(length(range_uniform)-Nt,size(obj.eg.image_data,2))];
  for idx = 1:length(surface)
    range = min(time,surface(idx)) * c/2 ...
      + max(0,time-surface(idx)) * c/(sqrt(er_ice)*2);
    obj.eg.image_data(:,idx) = interp1(range,...
      obj.eg.image_data(1:Nt,idx),range_uniform,'linear');
  end
  % update image_yaxis
  obj.eg.image_yaxis = range_uniform;
  % update y label
  obj.eg.y_label = 'Depth/Range(m)';
  obj.eg.y_order = 'reverse';
  % convert layerPnts_y from TWTT to depth
  layer_y_curUnit = obj.eg.layer.y;
  for idx = 1:length(layer_y_curUnit)
    surface = interp1(obj.eg.gps_time,...
      obj.eg.surface,...
      obj.eg.layer.x{idx},'linear');
    for pnt_idx = 1:length(layer_y_curUnit{idx})
      layer_y_curUnit{idx}(pnt_idx) = min(layer_y_curUnit{idx}(pnt_idx),surface(pnt_idx))*c/2 ...
        +max(0,layer_y_curUnit{idx}(pnt_idx)-surface(pnt_idx)) * c/(sqrt(er_ice)*2);
    end
  end
  % convert crossovers too
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surface,...
    obj.eg.crossovers.gps_time,'linear');
  obj.eg.crossovers.elev_offset = 2*(obj.eg.crossovers.cross_elev-obj.eg.crossovers.source_elev)./c;
  % correct surface using elevation offset info (since we are using surface
  % for frame1 and not frame2 and the twtt of the crossover is from frame2)
  % doesn't work right now for some reason
  surface = surface + obj.eg.crossovers.elev_offset;
  obj.eg.crossovers.y_curUnit = min(obj.eg.crossovers.twtt,surface).*c/2 ...
        + max(0,obj.eg.crossovers.twtt-surface).*c/(sqrt(er_ice)*2);
elseif yaxis_choice == 4 % Range bin
  % update image_yaxis
  obj.eg.image_yaxis = 1:length(obj.eg.time);
  % update y label
  obj.eg.y_label = 'Range bin';
  obj.eg.y_order = 'reverse';
  % convert gCtrl layerPnts_y from TWTT to range bin
  layer_y_curUnit = obj.eg.layer.y;
  for idx = 1:length(layer_y_curUnit)
    layer_y_curUnit{idx} = interp1(obj.eg.time,...
      1:length(obj.eg.time),...
      layer_y_curUnit{idx},'linear');
  end
  % convert crossovers too
  obj.eg.crossovers.elev_offset = 2*(obj.eg.crossovers.cross_elev-obj.eg.crossovers.source_elev)./c;
  % calculate crossover twtts with elevation differences compensated
  corrected_twtt = obj.eg.crossovers.twtt + obj.eg.crossovers.elev_offset;
  obj.eg.crossovers.y_curUnit = interp1(obj.eg.time,1:length(obj.eg.time),corrected_twtt);
end

% ======================================================================
%% We plot the whole data matrix and then use xlim and ylim to control the
% limits of what is displayed
set(obj.eg.h_image,{'XData','YData'},obj.eg.image_xaxis,obj.eg.image_yaxis);
obj.left_panel.imagewin.set_cdata(obj, obj.eg.image_data);
zoom reset;

%% Add labels, y-direction of axis, title, and set the colormap
xlabel(obj.right_panel.axes.handle,obj.eg.x_label);
ylabel(obj.right_panel.axes.handle,obj.eg.y_label);
set(obj.right_panel.axes.handle,'YDir',obj.eg.y_order);
title_str = sprintf('%s %d to %d',obj.eg.cur_sel.day_seg, ...
  obj.eg.frame_idxs(1), obj.eg.frame_idxs(end));
title(title_str,'Interpreter','none','Parent',obj.right_panel.axes.handle);

%% Set axis limits
x_lims = interp1(obj.eg.image_gps_time,obj.eg.image_xaxis,[x_min x_max],'linear','extrap');
x_lim_low = max([min(x_lims) obj.eg.image_xaxis(1)]);
x_lim_high = min([max(x_lims) obj.eg.image_xaxis(end)]);
xlim(obj.right_panel.axes.handle,[x_lim_low,x_lim_high]);
if y_min == -inf
  y_min = min(obj.eg.image_yaxis([1 end]));
end
if y_max == inf
  y_max = max(obj.eg.image_yaxis([1 end]));
end
ylim(obj.right_panel.axes.handle,sort([y_min y_max]))

%% Plot layers
%% WARNING: DO NOT IMPLEMEMNT WITH SCATTER... TOO SLOW RENDERING
layer_data_x = obj.eg.layer.x;
for idx = 1:length(layer_data_x)
  % Convert x-axis units
  layer_x_curUnit = interp1(obj.eg.image_gps_time,...
    obj.eg.image_xaxis,...
    layer_data_x{idx},'linear');
  
  % get manual/auto pts (use them for layer handles)
  layer_manual = obj.eg.layer.type{idx} == 1;
  
  % Manual points (plot this way to handle empty XData or YData
  obj.layer_h(2*(idx-1)+1) = plot(obj.right_panel.axes.handle,1,1,'bx');
  set(obj.layer_h(2*(idx-1)+1),{'XData','YData'}, ...
    {layer_x_curUnit(layer_manual),layer_y_curUnit{idx}(layer_manual)});
  % Auto and manual points
  obj.layer_h(2*(idx-1)+2) = plot(obj.right_panel.axes.handle, ...
    layer_x_curUnit,layer_y_curUnit{idx},'b--');
  
  layer_good = obj.eg.layer.qual{idx} == 1;
  layer_moderate = obj.eg.layer.qual{idx} == 2;
  layer_derived = obj.eg.layer.qual{idx} == 3;
  layer_y_curUnit_good = layer_y_curUnit{idx};
  layer_y_curUnit_good(~layer_good) = NaN;
  layer_y_curUnit_moderate= layer_y_curUnit{idx};
  layer_y_curUnit_moderate(~layer_moderate) = NaN;
  layer_y_curUnit_derived= layer_y_curUnit{idx};
  layer_y_curUnit_derived(~layer_derived) = NaN;

  % Good manual points (plot this way to handle empty XData or YData
  obj.quality_h(6*(idx-1)+1) = plot(obj.right_panel.axes.handle,1,1,'gx');
  set(obj.quality_h(6*(idx-1)+1),{'XData','YData'}, ...
    {layer_x_curUnit(layer_good&layer_manual),layer_y_curUnit{idx}(layer_good&layer_manual)});
  
  % Good auto points (plot this way to handle empty XData or YData
  obj.quality_h(6*(idx-1)+2) = plot(obj.right_panel.axes.handle,1,1,'g--');
  set(obj.quality_h(6*(idx-1)+2),{'XData','YData'}, ...
    {layer_x_curUnit,layer_y_curUnit_good});

  % Moderate manual points (plot this way to handle empty XData or YData
  obj.quality_h(6*(idx-1)+3) = plot(obj.right_panel.axes.handle,1,1,'yx');
  set(obj.quality_h(6*(idx-1)+3),{'XData','YData'}, ...
    {layer_x_curUnit(layer_moderate&layer_manual),layer_y_curUnit{idx}(layer_moderate&layer_manual)});
  
  % Moderate auto points (plot this way to handle empty XData or YData
  obj.quality_h(6*(idx-1)+4) = plot(obj.right_panel.axes.handle,1,1,'y--');
  set(obj.quality_h(6*(idx-1)+4),{'XData','YData'}, ...
    {layer_x_curUnit,layer_y_curUnit_moderate});

  % Derived manual points (plot this way to handle empty XData or YData
  obj.quality_h(6*(idx-1)+5) = plot(obj.right_panel.axes.handle,1,1,'rx');
  set(obj.quality_h(6*(idx-1)+5),{'XData','YData'}, ...
    {layer_x_curUnit(layer_derived&layer_manual),layer_y_curUnit{idx}(layer_derived&layer_manual)});
  
  % Derived auto points (plot this way to handle empty XData or YData
  obj.quality_h(6*(idx-1)+6) = plot(obj.right_panel.axes.handle,1,1,'r--');
  set(obj.quality_h(6*(idx-1)+6),{'XData','YData'}, ...
    {layer_x_curUnit,layer_y_curUnit_derived});

  %% Update obj.eg layers with current units
  obj.eg.layer.x_curUnit{idx} = layer_x_curUnit;
  obj.eg.layer.y_curUnit{idx} = layer_y_curUnit{idx};
end

%% Plot crossover data
obj.eg.crossovers.h = [];
for idx = 1:length(obj.eg.crossovers.x_curUnit)
  obj.eg.crossovers.h(2*(idx-1)+1) = plot(obj.right_panel.axes.handle, ...
    obj.eg.crossovers.x_curUnit(idx),obj.eg.crossovers.y_curUnit(idx),'o', ...
    'markersize',25);
  obj.eg.crossovers.h(2*(idx-1)+2) = plot(obj.right_panel.axes.handle, ...
    obj.eg.crossovers.x_curUnit(idx),obj.eg.crossovers.y_curUnit(idx),'x', ...
    'markersize',25);
end

%% Set layer and cross over visibility
obj.set_visibility();

%% Plot cursor (if valid)
xlimits = xlim(obj.right_panel.axes.handle);
start_gps = obj.eg.image_gps_time(find(obj.eg.image_xaxis>=xlimits(1),1));
stop_gps = obj.eg.image_gps_time(find(obj.eg.image_xaxis<=xlimits(2),1,'last'));

if ~isempty(obj.cursor.gps_time) ...
    && obj.cursor.gps_time >= start_gps && obj.cursor.gps_time <= stop_gps
  obj.echo_update_cursors();
  str = obj.status_text_cursor();
  obj.status_text_set(str,'replace');
end

%% Plot frame indicator (color lines)
% obj.frame_indicator_h = [];
% picker_pick_frame_indicator(curHandle,obj.right_panel.axes.handle);

%% Apply the display mode
obj.left_panel.imagewin.set_cdata(obj.eg.image_data);

%% Plot new selection on flight path
notify(obj,'update_echowin_flightline');

return;

