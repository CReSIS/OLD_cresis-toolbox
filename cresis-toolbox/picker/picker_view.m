function picker_view(cmd,param)
% picker_view(cmd,param)
%
% Runs commands related to the picker viewing window.
%
% Called from picker.m and support functions.

global gCtrl;
global hui;

% =================================================================
% =================================================================
% Load a new dataset into the view window
% =================================================================
% =================================================================
if cmd == 1
  % Create filename based on source directory and current selection
  fn = gCtrl.source.src_fns{gCtrl.source.cur_sel}{gCtrl.source.cur_src};
  [tmp name] = fileparts(fn);
  
  if ~exist(fn,'file')
    return;
  end
  
  % If only the loaded view file is the same frame, keep axis, cursor
  % etc. the same.
  if gCtrl.source.cur_view ~= gCtrl.source.cur_sel || ~ishandle(hui.viewfig.handle)
    redo_figure = true;
  else
    redo_figure = false;
  end
  
  % Select the figure for plotting
  figure(hui.viewfig.handle);
  
  if ~redo_figure
    % Store critical figure information
    cur_axis = axis;
  else
    % Redoing figure, reset information
    gCtrl.view.cur_idx = 1;
    set(hui.viewfig.handle,'WindowButtonUpFcn',@picker_view_button);
    set(hui.viewfig.handle,'WindowKeyReleaseFcn',@picker_view_key);
  end
  
  clf;
  set(hui.viewfig.handle,'NumberTitle','off');
  set(hui.viewfig.handle,'Name','4: View');
  
  % Load layers
  tmp = load(fn);
  if isfield(tmp,'Truncate_Bins')
    % This is a FMCW compressed data file.  See kuband_readme.doc or
    % snow_readme.doc for details.
    tmp.Time = tmp.Time(tmp.Truncate_Bins);
    tmp.Depth = tmp.Depth(tmp.Truncate_Bins);
  end
  
  % Keep track of important variables used for cursor, interpolation, etc.
  gCtrl.view.time = tmp.Time;
  gCtrl.view.gps_time = tmp.GPS_time;
  gCtrl.view.data = double(tmp.Data);
  
  gCtrl.view.xaxis = 1:length(gCtrl.source.geo{gCtrl.source.cur_sel}.GPS_time);
  
  % -----------------------------------------------------------------------
  % Interpolate data so it matches layer file
  % 1. GPS times may be different so interpolate across range lines
  %    to line these up
  if length(gCtrl.view.gps_time) ~= length(gCtrl.source.geo{gCtrl.source.cur_sel}.GPS_time) ...
    || any(gCtrl.view.gps_time ~= gCtrl.source.geo{gCtrl.source.cur_sel}.GPS_time)
    gCtrl.view.data = interp1(gCtrl.view.gps_time,gCtrl.view.data.', ...
      gCtrl.source.geo{gCtrl.source.cur_sel}.GPS_time).';
    % nan_idxs = find(isnan(gCtrl.view.data(1,:))==1);
    % if ~isempty(nan_idxs)
    %   gCtrl.view.data(:,nan_idxs) = [];
    %   gCtrl.view.xaxis(nan_idxs) = [];
    % end
  end
  % 2. Motion compensation may be different so adjust fast-time
  if isfield(gCtrl.source.geo{gCtrl.source.cur_sel},'Elevation')
    oldElevation = interp1(gCtrl.view.gps_time,tmp.Elevation, ...
      gCtrl.source.geo{gCtrl.source.cur_sel}.GPS_time);
    physical_constants;
    dr = (tmp.Time(2)-tmp.Time(1)) * c/2;
    for rline = 1:size(gCtrl.view.data,2)
      range_corr = oldElevation(rline)-gCtrl.source.geo{gCtrl.source.cur_sel}.Elevation(rline);
      if isfinite(range_corr) && abs(range_corr) > dr/2
        gCtrl.view.data(:,rline) = circshift(gCtrl.view.data(:,rline),-round(range_corr/dr));
      end
    end
  end
  
  % -----------------------------------------------------------------------
  % Load up the image and the cursor
  display_mode = get(hui.fig.ctrl_panel.display_modePM,'Value');
  if display_mode == 1
    hui.viewfig.image.h = imagesc(gCtrl.view.xaxis, gCtrl.view.time*1e6, lp(gCtrl.view.data,1));
  elseif display_mode == 2
    hui.viewfig.image.h = imagesc(gCtrl.view.xaxis, gCtrl.view.time*1e6, lp(local_detrend(gCtrl.view.data,[1 1],[1 10],5),1));
  elseif display_mode == 3
    hui.viewfig.image.h = imagesc(gCtrl.view.xaxis, gCtrl.view.time*1e6, lp(local_detrend(gCtrl.view.data,[1 1],[3 10],3),1));
  elseif display_mode == 4
    hui.viewfig.image.h = imagesc(gCtrl.view.xaxis, gCtrl.view.time*1e6, lp(local_detrend(gCtrl.view.data,[31 101],[3 10],4),1));
  elseif display_mode == 5
    hui.pickfig.image.h = imagesc(gCtrl.view.xaxis, gCtrl.view.time*1e6, detrending(gCtrl.view.data)); 
  end
  xlabel('Range line');
  ylabel('Two-way Propagation (us)');
  picker_colormap(hui.viewfig.handle);
  hui.viewfig.axes.handle = gca;
  title(sprintf('Data Frame ID: %s, Source: %s', gCtrl.source.frm_id(gCtrl.source.cur_sel,:), ...
    gCtrl.source.src_disp{gCtrl.source.cur_sel}{gCtrl.source.cur_src}), ...
    'Interpreter','none');
  hold on;
  hui.viewfig.cursor.h = plot([gCtrl.view.cursor gCtrl.view.cursor], ...
    [gCtrl.view.time([1 end])*1e6], 'k--');
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
  hui.viewfig.layer_h(1) = plot(1e6*gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{1}.value{1}.data,'m.');
  hui.viewfig.layer_h(2) = plot(1e6*gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{1}.value{2}.data,'m--');
  hui.viewfig.layer_h(3) = plot(1e6*gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{2}.value{1}.data,'r.');
  hui.viewfig.layer_h(4) = plot(1e6*gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{2}.value{2}.data,'r--');
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
  hui.viewfig.quality_h(1) = plot(1e6*qual_tmp,'go');
  qual_tmp = qual_plot;
  qual_tmp(qual~=2) = inf;
  hui.viewfig.quality_h(2) = plot(1e6*qual_tmp,'yo');
  qual_tmp = qual_plot;
  qual_tmp(qual~=3) = inf;
  hui.viewfig.quality_h(3) = plot(1e6*qual_tmp,'ro');
  % Quality for bottom (3 plots)
  manual = gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{2}.value{1}.data;
  derived = gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{2}.value{2}.data;
  qual_plot = derived;
  qual_plot(isfinite(manual)) = manual(isfinite(manual));
  qual = gCtrl.source.layers{gCtrl.source.cur_sel}.layerData{2}.quality;
  qual_tmp = qual_plot;
  qual_tmp(qual~=1) = inf;
  hui.viewfig.quality_h(4) = plot(1e6*qual_tmp,'go');
  qual_tmp = qual_plot;
  qual_tmp(qual~=2) = inf;
  hui.viewfig.quality_h(5) = plot(1e6*qual_tmp,'yo');
  qual_tmp = qual_plot;
  qual_tmp(qual~=3) = inf;
  hui.viewfig.quality_h(6) = plot(1e6*qual_tmp,'ro');
  hold off;
  
  picker_view_layer_visible;
  
  % Successful, so update cur_view and map
  picker_map(4,gCtrl.source.cur_sel);
  
  if redo_figure
    xlim([1 length(gCtrl.source.geo{gCtrl.source.cur_view}.GPS_time)]);
  end
  
elseif cmd == 2 && ishandle(hui.viewfig.handle)
  % Load up the image and the cursor
  display_mode = get(hui.fig.ctrl_panel.display_modePM,'Value');
  if display_mode == 1
    set(hui.viewfig.image.h,'CData',lp(gCtrl.view.data,1));
  elseif display_mode == 2
    set(hui.viewfig.image.h,'CData',lp(local_detrend(gCtrl.view.data,[1 1],[1 10],5),1));
  elseif display_mode == 3
    set(hui.viewfig.image.h,'CData',lp(local_detrend(gCtrl.view.data,[1 1],[3 10],3),1));
  elseif display_mode == 4
    set(hui.viewfig.image.h,'CData',lp(local_detrend(gCtrl.view.data,[31 101],[3 10],3),1));
  elseif display_mode == 5
    set(hui.viewfig.image.h,'CData',detrending(gCtrl.view.data));
  end
  
elseif cmd == 3 && ishandle(hui.viewfig.handle)
  % Set y-limits
  figure(hui.viewfig.handle);
  new_ylim = picker_ylimits(get(hui.fig.ctrl_panel.ylim_TE,'String'), ylim);
  ylim(new_ylim);
  
end

return;

