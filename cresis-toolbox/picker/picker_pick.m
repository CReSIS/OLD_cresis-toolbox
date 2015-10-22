function picker_pick(cmd,param)
% picker_pick(cmd,param)
%
% Runs commands related to the picker picking window.
%
% Called from picker.m and support functions.

global gCtrl;
global hui;

% =================================================================
% =================================================================
% Load a new dataset into the pick window
% =================================================================
% =================================================================
if cmd == 1
  % Create filename based on source directory and current selection
  fn = gCtrl.source.src_fns{gCtrl.source.cur_sel}{gCtrl.source.cur_src};
  [tmp name] = fileparts(fn);

  if ~exist(fn,'file')
    fprintf('  File does not exist\n');
    return;
  end
  
  % Update layer information from file if this current instance of the
  % picker has not modified this frame before
  if gCtrl.source.modified(gCtrl.source.cur_sel) == ' '
    source_fn = gCtrl.source.fns{gCtrl.source.cur_sel};
    load(source_fn,'layerData');
    if length(gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{1}.quality) == length(layerData{1}.quality)
      gCtrl.source.layers{gCtrl.source.cur_sel}.layerData = layerData;
    end
  end
  
  % Load data and metadata
  tmp = load(fn);
  if isfield(tmp,'Truncate_Bins')
    % This is a FMCW compressed data file.  See kuband_readme.doc or
    % snow_readme.doc for details.
    tmp.Time = tmp.Time(tmp.Truncate_Bins);
    tmp.Depth = tmp.Depth(tmp.Truncate_Bins);
  end
  
  % If only the loaded pick file is the same frame, keep axis, cursor
  % etc. the same.
  if gCtrl.source.cur_pick ~= gCtrl.source.cur_sel || ~ishandle(hui.pickfig.handle)
    redo_figure = true;
    gCtrl.undo_stack_ptr = 0;
    gCtrl.undo_stack = [];
  else
    redo_figure = false;
  end
  if ~ishandle(hui.pickfig.handle)
    new_plot = true;
  else
    new_plot = false;
  end
  
  % If window is already open, save current layers to memory
  if ishandle(hui.pickfig.handle)
    gCtrl.source.layers{gCtrl.source.cur_pick}.layerData{1}.value{1}.data ...
      = get(hui.pickfig.layer_h(1),'YData') / 1e6;
    gCtrl.source.layers{gCtrl.source.cur_pick}.layerData{1}.value{2}.data ...
      = get(hui.pickfig.layer_h(2),'YData') / 1e6;
    gCtrl.source.layers{gCtrl.source.cur_pick}.layerData{2}.value{1}.data ...
      = get(hui.pickfig.layer_h(3),'YData') / 1e6;
    gCtrl.source.layers{gCtrl.source.cur_pick}.layerData{2}.value{2}.data ...
      = get(hui.pickfig.layer_h(4),'YData') / 1e6;
    
    for cur_layer = 1:2
      good = get(hui.pickfig.quality_h(3*(cur_layer-1)+1),'YData');
      moderate = get(hui.pickfig.quality_h(3*(cur_layer-1)+2),'YData');
      derived = get(hui.pickfig.quality_h(3*(cur_layer-1)+3),'YData');
      qual = 1*isfinite(good) + 2*isfinite(moderate) + 3*isfinite(derived);
      qual(qual < 1 | qual > 3) = 1;
      gCtrl.source.layers{gCtrl.source.cur_pick}.layerData{cur_layer}.quality = qual;
    end
  end
  
  % Select the figure for plotting
  figure(hui.pickfig.handle);
  
  if ~redo_figure
    % Store critical figure information
    cur_axis = axis;
  else
    % Redoing figure, reset information
    gCtrl.pick.cur_idx = 1;
    set(hui.pickfig.handle,'WindowButtonUpFcn',@picker_pick_button);
    set(hui.pickfig.handle,'WindowKeyReleaseFcn',@picker_pick_key);
  end
  
  clf;
  set(hui.pickfig.handle,'NumberTitle','off');
  set(hui.pickfig.handle,'Name','3: Pick');
  
  % Keep track of important variables used for cursor, interpolation, etc.
  gCtrl.pick.time = tmp.Time;
  gCtrl.pick.gps_time = tmp.GPS_time;
  gCtrl.pick.data = double(tmp.Data);
  
  gCtrl.pick.xaxis = 1:length(gCtrl.source.geo{gCtrl.source.cur_sel}.GPS_time);
  
  % -----------------------------------------------------------------------
  % Interpolate data so it matches layer file
  % 1. GPS times may be different so interpolate across range lines
  %    to line these up
  if length(gCtrl.pick.gps_time) ~= length(gCtrl.source.geo{gCtrl.source.cur_sel}.GPS_time) ...
    || any(gCtrl.pick.gps_time ~= gCtrl.source.geo{gCtrl.source.cur_sel}.GPS_time)
    gCtrl.pick.data = interp1(gCtrl.pick.gps_time,gCtrl.pick.data.', ...
      gCtrl.source.geo{gCtrl.source.cur_sel}.GPS_time).';
    % nan_idxs = find(isnan(gCtrl.pick.data(1,:))==1);
    % if ~isempty(nan_idxs)
    %   gCtrl.pick.data(:,nan_idxs) = [];
    %   gCtrl.pick.xaxis(nan_idxs) = [];
    % end
  end
  % 2. Motion compensation may be different so adjust fast-time
  if isfield(gCtrl.source.geo{gCtrl.source.cur_sel},'Elevation')
    oldElevation = interp1(gCtrl.pick.gps_time,tmp.Elevation, ...
      gCtrl.source.geo{gCtrl.source.cur_sel}.GPS_time);
    physical_constants;
    dr = (tmp.Time(2)-tmp.Time(1)) * c/2;
    for rline = 1:size(gCtrl.pick.data,2)
      range_corr = oldElevation(rline)-gCtrl.source.geo{gCtrl.source.cur_sel}.Elevation(rline);
      if isfinite(range_corr) && abs(range_corr) > dr/2
        gCtrl.pick.data(:,rline) = circshift(gCtrl.pick.data(:,rline),-round(range_corr/dr));
      end
    end
  end
  
  % Load up the image and the cursor
  display_mode = get(hui.fig.ctrl_panel.display_modePM,'Value');
  if display_mode == 1
    hui.pickfig.image.h = imagesc(gCtrl.pick.xaxis, gCtrl.pick.time*1e6, lp(gCtrl.pick.data,1));
  elseif display_mode == 2
    hui.pickfig.image.h = imagesc(gCtrl.pick.xaxis, gCtrl.pick.time*1e6, lp(local_detrend(gCtrl.pick.data,[1 1],[1 10],5),1));
  elseif display_mode == 3
    hui.pickfig.image.h = imagesc(gCtrl.pick.xaxis, gCtrl.pick.time*1e6, lp(local_detrend(gCtrl.pick.data,[1 1],[3 10],3),1));
  elseif display_mode == 4
    hui.pickfig.image.h = imagesc(gCtrl.pick.xaxis, gCtrl.pick.time*1e6, lp(local_detrend(gCtrl.pick.data,[31 101],[3 10],4),1));
  elseif display_mode == 5
    hui.pickfig.image.h = imagesc(gCtrl.pick.xaxis, gCtrl.pick.time*1e6, detrending(gCtrl.pick.data)); 
  end
  xlabel('Range line');
  ylabel('Two-way Propagation (us)');
  picker_colormap(hui.pickfig.handle);
  hui.pickfig.axes.handle = gca;
  title(sprintf('Data Frame ID: %s, Source: %s %s', gCtrl.source.frm_id(gCtrl.source.cur_sel,:), ...
    gCtrl.source.src_disp{gCtrl.source.cur_sel}{gCtrl.source.cur_src}, char(gCtrl.source.modified(gCtrl.source.cur_sel))), ...
    'Interpreter','none');
  hold on;
  cursor_xpos = interp1(gCtrl.source.geo{gCtrl.source.cur_sel}.GPS_time, ...
    1:length(gCtrl.source.geo{gCtrl.source.cur_sel}.GPS_time), gCtrl.pick.cursor, 'linear','extrap');
  hui.pickfig.cursor.h = plot([cursor_xpos cursor_xpos], ...
    [gCtrl.pick.time([1 end])*1e6], 'k--');
  hold off;
  
  if ~redo_figure
    axis(cur_axis);
  else
    new_ylim = picker_ylimits(get(hui.fig.ctrl_panel.ylim_TE,'String'), ylim);
    ylim(new_ylim);
  end
  
  hold on;
  % --------------------------------------------------------------------
  % Plot layers onto data if they exist
  hui.pickfig.layer_h(1) = plot(1e6*gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{1}.value{1}.data,'m.');
  hui.pickfig.layer_h(2) = plot(1e6*gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{1}.value{2}.data,'m--');
  hui.pickfig.layer_h(3) = plot(1e6*gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{2}.value{1}.data,'r.');
  hui.pickfig.layer_h(4) = plot(1e6*gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{2}.value{2}.data,'r--');
  % --------------------------------------------------------------------
  % Plot quality level onto data if they exist
  % Quality for surface (3 plots)
  manual = gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{1}.value{1}.data;
  derived = gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{1}.value{2}.data;
  qual_plot = derived;
  qual_plot(isfinite(manual)) = manual(isfinite(manual));
  qual = gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{1}.quality;
  if length(qual) ~= length(qual_plot)
    fprintf('Quality length mismatch with layer data (ERROR!)\n');
    keyboard;
  end
  qual_tmp = qual_plot;
  qual_tmp(qual~=1) = inf;
  hui.pickfig.quality_h(1) = plot(1e6*qual_tmp,'go');
  qual_tmp = qual_plot;
  qual_tmp(qual~=2) = inf;
  hui.pickfig.quality_h(2) = plot(1e6*qual_tmp,'yo');
  qual_tmp = qual_plot;
  qual_tmp(qual~=3) = inf;
  hui.pickfig.quality_h(3) = plot(1e6*qual_tmp,'ro');
  % Quality for bottom (3 plots)
  manual = gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{2}.value{1}.data;
  derived = gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{2}.value{2}.data;
  qual_plot = derived;
  qual_plot(isfinite(manual)) = manual(isfinite(manual));
  qual = gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{2}.quality;
  qual_tmp = qual_plot;
  qual_tmp(qual~=1) = inf;
  hui.pickfig.quality_h(4) = plot(1e6*qual_tmp,'go');
  qual_tmp = qual_plot;
  qual_tmp(qual~=2) = inf;
  hui.pickfig.quality_h(5) = plot(1e6*qual_tmp,'yo');
  qual_tmp = qual_plot;
  qual_tmp(qual~=3) = inf;
  hui.pickfig.quality_h(6) = plot(1e6*qual_tmp,'ro');
  hold off;
  
  picker_pick_layer_visible;
  
  
  % Successful, so update cur_pick and map
  picker_map(3,gCtrl.source.cur_sel);
  
  if redo_figure
    xlim([1 length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time)]);
  end
  
  % --------------------------------------------------------------------
  % Plot landmarks on figure if using landmark tool
  if get(hui.fig.ctrl_panel.toolPM,'Value') == 10
    picker_pick_landmarks;
  end
  
elseif cmd == 2 && ishandle(hui.pickfig.handle)
  % Change the image processing or display mode
  display_mode = get(hui.fig.ctrl_panel.display_modePM,'Value');
  if display_mode == 1
    set(hui.pickfig.image.h,'CData',lp(gCtrl.pick.data,1));
  elseif display_mode == 2
    set(hui.pickfig.image.h,'CData',lp(local_detrend(gCtrl.pick.data,[1 1],[1 10],5),1));
  elseif display_mode == 3
    set(hui.pickfig.image.h,'CData',lp(local_detrend(gCtrl.pick.data,[1 1],[3 10],3),1));
  elseif display_mode == 4
    set(hui.pickfig.image.h,'CData',lp(local_detrend(gCtrl.pick.data,[31 101],[3 10],4),1));
  elseif display_mode == 5
    set(hui.pickfig.image.h,'CData',detrending(gCtrl.pick.data));
  end
  
elseif cmd == 3 && ishandle(hui.pickfig.handle)
  % Set y-limits
  figure(hui.pickfig.handle);
  new_ylim = picker_ylimits(get(hui.fig.ctrl_panel.ylim_TE,'String'), ylim);
  ylim(new_ylim);
  
end

return;




