function create_movie (params, options)
% tomo.create_movie (params, options)
%
% Generate 3D "flyover" movies from existing DEM files.
% Usually called from tomo.run_create_movie
%
% Authors: Victor Berger, John Paden
%
% See also: tomo.run_create_movie.m, tomo.plot_surf_geotiff.m,
%  and tomo.freezeColors.m

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s (%s)\n', dbstack_info(1).name, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

physical_constants;

if options.showMAP
  fprintf('\nLoading GeoTIFF, this may take a minute.. (%s)\n', datestr(now));
  mapper.proj = geotiffinfo(options.pathMAP);
  [mapper.RGB_bg, mapper.R_bg, ~] = geotiffread(options.pathMAP);
  fprintf('  Done (%s)\n', datestr(now));
end

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) ...
      || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  options.videoPathTMP = ct_filename_out(param, options.videoPath);
  
  if ~exist(options.videoPathTMP, 'dir')
    mkdir(options.videoPathTMP);
  end
  
  % Load frames file
  load(ct_filename_support(param,'','frames'));
  
  if isempty(param.cmd.frms)
    param.cmd.frms = 1:length(frames.frame_idxs);
  end
  
  for frm_idx = param.cmd.frms
    
    drawCO = options.showCrossovers;
    if options.createVideo
      options.videoName = sprintf('Movie_%s_%03d',param.day_seg, frm_idx);
      if exist(strcat(fullfile(options.videoPathTMP, options.videoName), '.mp4'), 'file') || ...
          exist(strcat(fullfile(options.videoPathTMP, options.videoName), '.avi'), 'file')
        fprintf('\nFile for %s already exists. Continuing to next...', sprintf('%s_%03d',param.day_seg, frm_idx));
        continue
      end
    end
    
    fprintf('\nLoading DEM, this may take a minute.. (%s)\n', datestr(now));
    try
      param.topDEM    = load(fullfile(ct_filename_out(param, options.DEM_source, ''), sprintf('%s_%03d_top',param.day_seg, frm_idx)));
      param.bottomDEM = load(fullfile(ct_filename_out(param, options.DEM_source, ''), sprintf('%s_%03d_bottom',param.day_seg, frm_idx)));
      records_fn      = ct_filename_support(param,'','records');
      records = load(records_fn,'gps_time');
    catch ME
      warning('Attention: DEM for reference frame %s not found. Skipping this frame.', sprintf('%s_%03d',param.day_seg, frm_idx));
      continue;
    end
    
    fprintf('  Done (%s)', datestr(now));
    
    f1 = figure;
    ax1 = gca;
    
    % Draw reference TOP surface on figure 1
    filt_elev = wiener2(param.topDEM.points.elev(:, 1:length(param.topDEM.points.x)),[1 5]);
    vert_scale = 2;
    s1 = surf(param.topDEM.points.x(:, 1:length(param.topDEM.points.x))/vert_scale, param.topDEM.points.y(:, 1:length(param.topDEM.points.x))/vert_scale, filt_elev, ...
      'FaceAlpha', .5, 'EdgeColor', 'none', 'Parent', ax1);
    
    hold all
    
    % Draw cross-over TOP surface(s) on figure 1
    if drawCO && exist('data', 'var') && isfield(data, 'frame')
      for surf_idx = 1:length(crossover_frms)
        if ~crossover_badmask(surf_idx) || strcmp(crossover_frms{surf_idx}{1},sprintf('%s_%03d',param.day_seg, frm_idx))
          continue;
        end
        surf(crossover_DEM{surf_idx}.topDEM.points.x(:, 1:length(crossover_DEM{surf_idx}.topDEM.points.x))/vert_scale, ...
          crossover_DEM{surf_idx}.topDEM.points.y(:, 1:length(crossover_DEM{surf_idx}.topDEM.points.x))/vert_scale, ...
          crossover_DEM{surf_idx}.topDEM.points.elev(:, 1:length(crossover_DEM{surf_idx}.topDEM.points.x)), ...
          'FaceAlpha', .5, 'EdgeColor', 'none', 'Parent', ax1);
        hold all
      end
    end
    
    warning off
    colormap gray;
    camlight
    shading interp
    hold all
    tomo.freezeColors
    
    % Draw reference BOTTOM surface on figure 1
    filt_elev = wiener2(param.bottomDEM.points.elev(:, 1:length(param.bottomDEM.points.x)),[1 5]);
    s2 = surfl(param.bottomDEM.points.x(:, 1:length(param.bottomDEM.points.x))/vert_scale, ...
      param.bottomDEM.points.y(:, 1:length(param.bottomDEM.points.x))/vert_scale, ...
      filt_elev(:, 1:length(param.bottomDEM.points.x)), 'light', 'Parent', ax1);
    
    % Draw cross-over BOTTOM surface(s) on figure 1
    if drawCO && exist('data', 'var') && isfield(data, 'frame')
      for surf_idx = 1:length(crossover_frms)
        if ~crossover_badmask(surf_idx) || strcmp(crossover_frms{surf_idx}{1},sprintf('%s_%03d',param.day_seg, frm_idx))
          continue;
        end
        
        surfl(crossover_DEM{surf_idx}.bottomDEM.points.x(:, 1:length(crossover_DEM{surf_idx}.bottomDEM.points.x))/vert_scale, ...
          crossover_DEM{surf_idx}.bottomDEM.points.y(:, 1:length(crossover_DEM{surf_idx}.bottomDEM.points.x))/vert_scale, ...
          crossover_DEM{surf_idx}.bottomDEM.points.elev(:, 1:length(crossover_DEM{surf_idx}.bottomDEM.points.x)), ...
          'light','Parent', ax1);
        shading(ax1, 'interp');
      end
    end
    
    camlight
    shading interp
    colormap(demcmap(param.bottomDEM.points.elev, 24));
    axis(ax1, 'off');
    
    [B,A] = butter(2, .002);
    top_elev = filtfilt(B, A, interp_finite(param.topDEM.points.elev(32, :)));
    bot_elev = filtfilt(B, A, interp_finite(param.bottomDEM.points.elev(32, :)));
    
    % Function call to plot_DEM generates map and DEM
    param.frm = frm_idx;
    param = merge_structs(param, options);
    param.radar_name = options.sys;
    [f2, ax2, cax2] = tomo.plot_DEM(param, 'gt', mapper, 'title', 'off');
    caxis(ax1, cax2);
    
    min_x = min( inf, min(param.bottomDEM.points.x(:)));
    max_x = max(-inf, max(param.bottomDEM.points.x(:)));
    min_y = min( inf, min(param.bottomDEM.points.y(:)));
    max_y = max(-inf, max(param.bottomDEM.points.y(:)));
    
    figure_dots_per_km = 20;
    map_buffer = 2e3;
    
    map_min_x = min_x-map_buffer;
    map_max_x = max_x+map_buffer;
    map_min_y = min_y-map_buffer;
    map_max_y = max_y+map_buffer;
    
    set(ax2,'Units','pixels');
    map_axes = get(ax2,'Position');
    set(f2,'Units','pixels');
    map_pos = get(f2,'Position');
    
    map_new_axes = map_axes;
    map_new_axes(3) = round(figure_dots_per_km*(map_max_x-map_min_x)/1e3);
    map_new_axes(4) = round(figure_dots_per_km*(map_max_y-map_min_y)/1e3);
    map_pos(3) = map_new_axes(3) + map_pos(3)-map_axes(3);
    map_pos(4) = map_new_axes(4) + map_pos(4)-map_axes(4);
    
    set(f2,'Position',map_pos);
    set(ax2,'Position',map_new_axes);
    set(f2,'PaperPositionMode','auto');
    
    camzoom(ax1, options.zoom)
    oset       = 300;
    axis_equal = 0;
    slice_idx  = 1;
    
    %% Pad figures to obtain total size
    titleFrameWidth = 40;
    
    f1.Position = [200+options.videoSize(1)/2, 200, options.videoSize(1)/2,...
      options.videoSize(2)-titleFrameWidth];
    f2.Position = [200, 200, options.videoSize(1)/2, ...
      options.videoSize(2)-titleFrameWidth];
    
    set(ax2, 'Units', 'normalized')
    set(ax2, 'OuterPosition', [.1, .1, .75, .75]);
    
    drawnow;
    
    % Generate line on 2D map figure
    h = line(param.bottomDEM.points.x(:, find(~isnan(param.bottomDEM.points.x(32,:)),1))/1e3,...
      param.bottomDEM.points.y(:, find(~isnan(param.bottomDEM.points.x(32,:)),1))/1e3, ...
      ones(size(param.bottomDEM.points.elev(:, slice_idx))),...
      'LineWidth', 2, 'Color','r','Parent', ax2);
    
    f3 = figure(3);
    f3.Position(3) = options.videoSize(1);
    f3.Position(4) = titleFrameWidth;
    ax3            = gca;
    axis(ax3, 'off');
    
    seg_year = str2double(param.day_seg(1:4));
    seg_month = str2double(param.day_seg(5:6));
    seg_day = str2double(param.day_seg(7:8));
    
    % Start time (used in title)
    ST = records.gps_time(frames.frame_idxs(frm_idx));
    % Stop time (used in title)
    if frm_idx == length(frames.frame_idxs)
      last_idx = length(records.gps_time);
    else
      last_idx = frames.frame_idxs(frm_idx+1)-1;
    end
    ET = records.gps_time(last_idx);
    
    [year,month,day,hour,minute,sec] = datevec(epoch_to_datenum(ST));
    days_offset = datenum(year,month,day) - datenum(seg_year,seg_month,seg_day);
    hour = hour + 24*days_offset;
    start_time_str = sprintf('%02d:%02d:%04.1f', hour, minute, sec);
    
    [year,month,day,hour,minute,sec] = datevec(epoch_to_datenum(ET));
    days_offset = datenum(year,month,day) - datenum(seg_year,seg_month,seg_day);
    hour = hour + 24*days_offset;
    stop_time_str = sprintf('%02d:%02d:%04.1f', hour, minute, sec);
    
    txt = title(sprintf('%s %s: "%s"  %s: %s to %s GPS', options.sys, param.season_name, ...
      param.cmd.mission_name, sprintf('%s_%03d',param.day_seg, frm_idx), start_time_str, ...
      stop_time_str),'Interpreter','none','FontWeight','normal','FontSize', 15);
    
    txt.Units = 'normalized';
    txt.HorizontalAlignment = 'center';
    txt.Position = [0.5000 -1.2000 0];
    
    fr3 = getframe(f3);
    close(f3);
    
    warning on
    
    if options.createVideo
      try
        v = VideoWriter(fullfile(options.videoPathTMP, options.videoName), options.videoFormat);
      catch ME
        warning('The selected video format is not supported. Using .avi instead.');
        warning('MPEG-4: only systems with Windows 7 and later, or Mac OS X 10.7 and later.');
        v = VideoWriter(fullfile(options.videoPathTMP, options.videoName));
      end
      v.FrameRate = options.frameRate;
      fprintf('\nWriting video to \n%s  (%s)\n', fullfile(options.videoPathTMP, options.videoName), datestr(now));
      open(v);
    end
    
    %% Move along-track through frame to create animation
    while slice_idx <= length(param.bottomDEM.points.x)-oset
      if isnan(param.topDEM.points.x(32, slice_idx)) || isnan(param.bottomDEM.points.x(32, slice_idx)) ...
          || isnan(param.bottomDEM.points.x(32, slice_idx+oset))
        slice_idx = slice_idx + 1;
        continue
      end
      
      % Update line location on 2D map figure
      h.XData = param.bottomDEM.points.x(:, slice_idx)/1e3;
      h.YData = param.bottomDEM.points.y(:, slice_idx)/1e3;
      h.ZData = ones(size(param.bottomDEM.points.elev(:, slice_idx)));
      
      hold on
      
      % Update camera position and target
      campos(ax1, [param.bottomDEM.points.x(32,slice_idx)/vert_scale, param.bottomDEM.points.y(32, slice_idx)/vert_scale, top_elev(slice_idx)+100]);
      camtarget(ax1, [param.bottomDEM.points.x(32,slice_idx+oset)/vert_scale, param.bottomDEM.points.y(32, slice_idx+oset)/vert_scale, bot_elev(slice_idx+oset)]);
      
      if ~axis_equal
        axis(ax1, 'equal')
        axis_equal = 1;
      end
      
      % Capture images, concatenate frames, and write to video
      if options.createVideo
        fr1 = getframe(f1);
        fr2 = getframe(f2);
        frame.colormap = [];
        frame.cdata = horzcat(fr2.cdata, fr1.cdata);
        frame.cdata = vertcat(fr3.cdata, frame.cdata);
        writeVideo(v, frame);
      else
        pause(0.002);
      end
      
      slice_idx = slice_idx + options.frameSkip;
    end
    
    if options.createVideo
      fprintf('  Done (%s)\n', datestr(now));
      close(v)
    end
    
    close all;
    clear f1 f2 fr1 fr2 frame flowline OPSparams WKT_polygon ST ET
    fprintf('\n=====================================================================\n');
  end
end

fprintf('\n\nALL DONE ::: %s\n\n', datetime('now'));
