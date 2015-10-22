function picker_pick_button(src,event)
% picker_pick_button(src,event)
%
% Support function for picker_pick (the pick figure window)
%
% Author: John Paden, Aric Beaver

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

[x,y,but] = get_mouse_info(hui.pickfig.handle,hui.pickfig.axes.handle);
y = y / gCtrl.tool.layer_multiple;

% Find the gps time of click
x_time = interp1(1:length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time), ...
  gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time, x, 'linear','extrap');
fprintf('Pick: x = %.2f (GPS %s), y = %.2f us, but = %d\n', x, ...
  datestr(epoch_to_datenum(x_time),'yyyymmdd HH:MM:SS.FFF'), y, but);

if but == 1
  % ===================================================================
  % Left mouse button: activate tool
  
  tool = get(hui.fig.ctrl_panel.toolPM,'Value');
  cur_layers = get(hui.fig.ctrl_panel.layerPM,'Value');
  if cur_layers == 3
    cur_layers = [1 2];
  end
  begin_command = true;
  frame_modified = false;
  
  % menuString{1} = '(b)rowse';
  % menuString{2} = '(e)nter pnt';
  % menuString{3} = '(m)ax pnt';
  % menuString{4} = '(d)elete';
  % menuString{5} = '(q)uality';
  % menuString{6} = '(i)nterp';
  % menuString{7} = '(s)nake';
  % menuString{8} = '(M)ax';
  if tool == 1
    % Browse tool (does not do anything)
    
  elseif tool == 2
    % ================================================================
    % Enter Point tool
    % ================================================================
    for cur_layer = cur_layers
      layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
      xlims = xlim(hui.pickfig.axes.handle);
      if x < xlims(1)
        x = xlims(1);
      elseif x > xlims(2)
        x = xlims(2);
      end
      x = round(x);
      if x < 1
        x = 1;
      elseif x > length(layer_data)
        x = length(layer_data);
      end
      rlines = round(x);
      % Update quality
      picker_update_layer(cur_layer,1,rlines,y,begin_command);
      begin_command = false;
      frame_modified = true;
    end
    
  elseif tool == 3
    % ================================================================
    % Max Point tool
    % ================================================================
 
    for cur_layer = cur_layers
      layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
      xlims = xlim(hui.pickfig.axes.handle);
      if x < xlims(1)
        x = xlims(1);
      elseif x > xlims(2)
        x = xlims(2);
      end
        x = round(x);
      if x < 1
        x = 1;
      elseif x > length(layer_data)
        x = length(layer_data);
      end
      rlines = round(x);
      tool_param1_str = get(hui.fig.ctrl_panel.tool_param1_TE,'String');
      search_range = [];
      try
      % Assumes the ylimit_str is a matlab expression that can be evaluated
        search_range = eval(sprintf('[%s]', tool_param1_str));
        if length(search_range) == 1
          search_range = -search_range:search_range;
        else
          search_range = round(min(search_range)):round(max(search_range));
        end
      end
      if isempty(search_range)
        search_range = -5:5;
      end

      % A is the area of interest 
      A = get(hui.pickfig.image.h,'CData');
      A = interp1(gCtrl.pick.xaxis,A.',1:length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time)).';
      center = find(gCtrl.pick.time >= y.*1e-6,1,'first');
      first_bin = center + search_range(1);
      last_bin = center + search_range(end);
      A = A(first_bin:last_bin,rlines);
      % Find the Max value  
      [tmp, idx_max] = max(A,[],1);
      idx_max = idx_max + first_bin - 1;
      val = interp1(1:length(gCtrl.pick.time),gCtrl.pick.time*1e6,idx_max,'linear','extrap');    
      picker_update_layer(cur_layer,1,rlines,val,begin_command);
      begin_command = false;
      frame_modified = true;
    end
    
  elseif tool == 4
    % ================================================================
    % Delete region tool
    % ================================================================
    if isnan(gCtrl.tool.selection.x)
      gCtrl.tool.selection.x = x;
      gCtrl.tool.selection.y = y;
    else
      for cur_layer = cur_layers
        layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
        [xbnds ybnds outside] = picker_tool_box_sel(gCtrl.tool.selection.x, x, ...
          gCtrl.tool.selection.y, y, axis(hui.pickfig.axes.handle), ...
          length(layer_data));
        rlines = round(xbnds(1)):round(xbnds(2));
        % Determine which points are in the rubberband/box
        ROI = layer_data(rlines);
        if outside
          rlines_changed = rlines;
        else
          rlines_changed = rlines(ROI>ybnds(1) & ROI<ybnds(2));
        end
        picker_update_layer(cur_layer,1,rlines_changed,inf*rlines_changed,begin_command);
        begin_command = false;
        % Remove derived points
        layer_data = get(hui.pickfig.layer_h(2*cur_layer),'YData');
        % Determine which points are in the rubberband/box
        ROI = layer_data(rlines);
        if outside
          rlines_changed = rlines;
        else
          rlines_changed = rlines(ROI>ybnds(1) & ROI<ybnds(2));
        end
        picker_update_layer(cur_layer,2,rlines_changed,inf*rlines_changed,begin_command);
      end
      gCtrl.tool.selection.x = NaN;
      frame_modified = true;
    end
    
  elseif tool == 5
    % ================================================================
    % Change quality region tool
    % ================================================================
    if isnan(gCtrl.tool.selection.x)
      gCtrl.tool.selection.x = x;
      gCtrl.tool.selection.y = y;
    else
      for cur_layer = cur_layers
        layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
        [xbnds ybnds outside] = picker_tool_box_sel(gCtrl.tool.selection.x, x, ...
          gCtrl.tool.selection.y, y, axis(hui.pickfig.axes.handle), ...
          length(layer_data));
        rlines = round(xbnds(1)):round(xbnds(2));
        % Determine which points are in the rubberband/box
        ROI = layer_data(rlines);
        if outside
          rlines_changed = rlines;
        else
          rlines_changed = rlines(ROI>ybnds(1) & ROI<ybnds(2));
        end
        picker_update_layer(cur_layer,1,rlines_changed,layer_data(rlines_changed),begin_command);
        begin_command = false;
        % Change quality of derived points
        layer_data = get(hui.pickfig.layer_h(2*cur_layer),'YData');
        ROI = layer_data(rlines);
        if outside
          rlines_changed = rlines;
        else
          rlines_changed = rlines(ROI>ybnds(1) & ROI<ybnds(2));
        end
        picker_update_layer(cur_layer,2,rlines_changed,layer_data(rlines_changed),begin_command);
      end
      gCtrl.tool.selection.x = NaN;
      frame_modified = true;
    end
    
  elseif tool == 6
    % ================================================================
    % Interpolate region tool
    % ================================================================
    if isnan(gCtrl.tool.selection.x)
      gCtrl.tool.selection.x = x;
      gCtrl.tool.selection.y = y;
    else
      for cur_layer = cur_layers
        layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
        [xbnds ybnds outside] = picker_tool_box_sel(gCtrl.tool.selection.x, x, ...
          gCtrl.tool.selection.y, y, axis(hui.pickfig.axes.handle), ...
          length(layer_data));
        % Get manual points that we will interpolate from
        rlines = round(xbnds(1)):round(xbnds(2));
        ROI = layer_data(rlines);
        if outside
          rlines_valid = rlines(isfinite(ROI));
        else
          rlines_valid = rlines(ROI>ybnds(1) & ROI<ybnds(2));
        end
        
        if length(rlines_valid) < 2
          fprintf('Less than 2 manual points to interpolate in specified region.');
        else
          values = layer_data(rlines_valid);
          rlines = rlines_valid(1):rlines_valid(end);

          % Interpolate onto derived points
          picker_update_layer(cur_layer,2,rlines, ...
            interp1(rlines_valid,values,rlines),begin_command);
          begin_command = false;
        end
      end
      gCtrl.tool.selection.x = NaN;
      frame_modified = true;
    end
    
  elseif tool == 7
    % ================================================================
    % Snake region tool
    % ================================================================
    if isnan(gCtrl.tool.selection.x)
      gCtrl.tool.selection.x = x;
      gCtrl.tool.selection.y = y;
    else
      for cur_layer = cur_layers
        layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
        [xbnds ybnds outside] = picker_tool_box_sel(gCtrl.tool.selection.x, x, ...
          gCtrl.tool.selection.y, y, axis(hui.pickfig.axes.handle), ...
          length(layer_data));
        % Get manual points that we will snake from
        rlines = round(xbnds(1)):round(xbnds(2));
        ROI = layer_data(rlines);
        if outside
          rlines_valid = rlines(isfinite(ROI));
        else
          rlines_valid = rlines(ROI>ybnds(1) & ROI<ybnds(2));
        end
        
        if length(rlines_valid) < 2
          fprintf('Less than 2 manual points to interpolate in specified region.');
        else
          values = layer_data(rlines_valid);
          rlines = rlines_valid(1):rlines_valid(end);

          % Run snake on values
          vals = picker_snake(rlines_valid,values);
  
          picker_update_layer(cur_layer,2,rlines,vals,begin_command);
          begin_command = false;
        end
      end
      gCtrl.tool.selection.x = NaN;
      frame_modified = true;
    end
    
  elseif tool == 8
    % ================================================================
    % Max region tool
    % ================================================================
    if isnan(gCtrl.tool.selection.x)
      gCtrl.tool.selection.x = x;
      gCtrl.tool.selection.y = y;
    else
      for cur_layer = cur_layers
        layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
        [xbnds ybnds outside] = picker_tool_box_sel(gCtrl.tool.selection.x, x, ...
          gCtrl.tool.selection.y, y, axis(hui.pickfig.axes.handle), ...
          length(layer_data));
        rlines = round(xbnds(1)):round(xbnds(2));
        % Obtain data that was selected
        A = get(hui.pickfig.image.h,'CData');
        A = interp1(gCtrl.pick.xaxis,A.',1:length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time)).';
        if outside
            first_bin = 1;
            last_bin = size(A,2);
        else
            first_bin = find(gCtrl.pick.time >= ybnds(1).*1e-6,1,'first');
            last_bin = find(gCtrl.pick.time >= ybnds(2).*1e-6,1,'first');
        end
        A = A(first_bin:last_bin,rlines);
        % Find the Max value  
        [tmp, idx_max] = max(A,[],1);
        idx_max = idx_max + first_bin - 1;

        vals = interp1(1:length(gCtrl.pick.time),gCtrl.pick.time*1e6,idx_max,'linear','extrap');
        picker_update_layer(cur_layer,2,rlines,vals,begin_command);
        begin_command = false;
      end
      gCtrl.tool.selection.x = NaN;        
      frame_modified = true;
    end
    
  elseif tool == 9
    % ================================================================
    % Leading edge tool
    % ================================================================
    if isnan(gCtrl.tool.selection.x)
      gCtrl.tool.selection.x = x;
      gCtrl.tool.selection.y = y;
    else
      for cur_layer = cur_layers
        layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
        [xbnds ybnds outside] = picker_tool_box_sel(gCtrl.tool.selection.x, x, ...
          gCtrl.tool.selection.y, y, axis(hui.pickfig.axes.handle), ...
          length(layer_data));
        % Get manual points that we will snake from
        rlines = round(xbnds(1)):round(xbnds(2));
        
        % Obtain data that was selected
        A = get(hui.pickfig.image.h,'CData');
        if outside
            first_bin = 1;
            last_bin = size(A,2);
        else
            first_bin = find(gCtrl.pick.time >= ybnds(1).*1e-6,1,'first');
            last_bin = find(gCtrl.pick.time >= ybnds(2).*1e-6,1,'first');
        end
        first_bin = min(first_bin,size(A,1)-5);
        last_bin = max(last_bin,first_bin+5);
        
        ROI = layer_data(rlines);
        if outside
          rlines_valid = rlines(isfinite(ROI));
        else
          rlines_valid = rlines(ROI>ybnds(1) & ROI<ybnds(2));
        end
        
        if length(rlines_valid) < 2
          fprintf('Less than 2 manual points to interpolate in specified region.\n');
        else
          values = layer_data(rlines_valid);
          rlines = rlines_valid(1):rlines_valid(end);

          % Run snake on values
          vals = picker_threshold(rlines_valid,values,first_bin,last_bin);
  
          picker_update_layer(cur_layer,2,rlines,vals,begin_command);
          begin_command = false;
        end
      end
      gCtrl.tool.selection.x = NaN;
      frame_modified = true;
    end
    
%   elseif tool == 10
%     % ================================================================
%     % Leading edge tool
%     % ================================================================
%     if isnan(gCtrl.tool.selection.x)
%       gCtrl.tool.selection.x = x;
%       gCtrl.tool.selection.y = y;
%     else
%       for cur_layer = cur_layers
%         layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
%         [xbnds ybnds outside] = picker_tool_box_sel(gCtrl.tool.selection.x, x, ...
%           gCtrl.tool.selection.y, y, axis(hui.pickfig.axes.handle), ...
%           length(layer_data));
%         rlines = round(xbnds(1)):round(xbnds(2));
%         % Obtain data that was selected
%         A = get(hui.pickfig.image.h,'CData');
%         %A = interp1(gCtrl.pick.xaxis,A.',1:length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time)).';
%         if outside
%             first_bin = 1;
%             last_bin = size(A,2);
%         else
%             first_bin = find(gCtrl.pick.time >= ybnds(1).*1e-6,1,'first');
%             last_bin = find(gCtrl.pick.time >= ybnds(2).*1e-6,1,'first');
%         end
%         first_bin = min(first_bin,size(A,1)-5);
%         last_bin = max(last_bin,first_bin+5);
%         A = A(first_bin:last_bin,rlines);
%         
%         tool_param1_str = get(hui.fig.ctrl_panel.tool_param1_TE,'String');
%         noise_pow_thresh = [];
%         num_lower_bins = [];
%         try
%           % Assumes the ylimit_str is a matlab expression that can be evaluated
%           eval_result = eval(sprintf('[%s]', tool_param1_str));
%           noise_pow_thresh = eval_result(1);
%           num_lower_bins = eval_result(2);
%         end
%         if isempty(noise_pow_thresh)
%           noise_pow_thresh = 10;
%         end
%         if isempty(num_lower_bins)
%           num_lower_bins = 3;
%         end
%         
%         % Find the background noise power
%         noise_pow = mean(mean(A(1:5,:)));
%         idx_max = inf*zeros(1,size(A,2));
%         for rline = 1:size(A,2)
%           % Find the leading edge
%           leading_edge_state = 0;
%           for bin = 1:size(A,1)-2
%             if leading_edge_state == 1
%               val_lower = true;
%               for bin_low = bin:bin+num_lower_bins-1
%                 if A(bin_low,rline) >= A(bin-1,rline)
%                   val_lower = false;
%                   break;
%                 end
%               end
%               if val_lower
%                 idx_max(rline) = bin + first_bin - 2;
%                 break;
%               end
%             end
%             if leading_edge_state == 0 && A(bin,rline) > noise_pow+noise_pow_thresh
%               leading_edge_state = 1;
%             end
%           end
%         end
%           vals = interp1(1:length(gCtrl.pick.time),gCtrl.pick.time*1e6,idx_max,'linear','extrap');
%           picker_update_layer(cur_layer,2,rlines,vals,begin_command);
%           begin_command = false;
%       end
%       gCtrl.tool.selection.x = NaN;        
%       frame_modified = true;
%     end
    
    elseif tool == 10    
    % ================================================================
    % Landmark requisition tool
    % ================================================================
    if isnan(gCtrl.tool.selection.x)
      gCtrl.tool.selection.x = x;
      gCtrl.tool.selection.y = y;
    else
      cur_layer = cur_layers(1);
      
      if isempty(gCtrl.source.landmarks.fn)
        [fn fn_dir] = uiputfile('*.mat','Select a landmark file','landmarks.mat');
        if isnumeric(fn)
          return;
        end
        gCtrl.source.landmarks.fn = fullfile(fn_dir,fn);
        if exist(gCtrl.source.landmarks.fn,'file')
          load(gCtrl.source.landmarks.fn)
        else
          landmarks = [];
          save(gCtrl.source.landmarks.fn,'landmarks');
        end
      end
      
      % Enable notes box
      set(hui.landfig.ctrl_panel.notesBox_TE,'Enable','on');
      set(hui.landfig.ctrl_panel.notesBox_TE,'ForegroundColor',[0 0 0])
      
      % Called to obtain landmark info
      layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
      [xbnds ybnds outside] = picker_tool_box_sel(gCtrl.tool.selection.x, x, ...
        gCtrl.tool.selection.y, y, axis(hui.pickfig.axes.handle), ...
        length(layer_data));
      
      rlines = round(xbnds(1)):round(xbnds(2));
      rline_start = rlines(1);
      rline_stop = rlines(end);
      rbin_start = ybnds(1);
      rbin_stop = ybnds(2);
      
      gpstime_start = interp1(1:length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time), ...
        gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time, xbnds(1), 'linear','extrap');
      gpstime_stop = interp1(1:length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time), ...
        gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time, xbnds(2), 'linear','extrap');
      
      frm_id = gCtrl.source.frm_id(gCtrl.source.cur_pick,:);
%       radar_name = 'radar_name';
%       season_name = 'season_name';
      radar_name = gCtrl.radar_name;
      season_name = gCtrl.season_name;
      notes = '';
      lm_types = get(hui.landfig.ctrl_panel.typePM,'String');
      lm_type = char(lm_types(get(hui.landfig.ctrl_panel.typePM,'Value')));
      
      
      if exist(gCtrl.source.landmarks.fn,'file')
        load(gCtrl.source.landmarks.fn)
      else
        landmarks = [];
      end
      
      new_landmark = struct('frm_id',{frm_id},'rbin_start',{rbin_start},'rbin_stop',{rbin_stop},...
        'rline_start',{rline_start},'rline_stop',{rline_stop},'gpstime_start',{gpstime_start},'gpstime_stop',{gpstime_stop},...
        'radar_name',{radar_name},'season_name',{season_name},'notes',{notes},'type',{lm_type});
      if ~isempty(landmarks)
        [tmp insert_idx] = sortrows(cat(1,landmarks.frm_id,frm_id));
        insert_idx = find(insert_idx == length(insert_idx));
        landmarks = [landmarks(1:insert_idx-1) new_landmark landmarks(insert_idx:end)];  
        save(gCtrl.source.landmarks.fn,'landmarks')
      else
        insert_idx = 1;
        landmarks = new_landmark;
        [path,name,ext] = fileparts(gCtrl.source.landmarks.fn);
        warning off; mkdir(path); warning on;
        save(gCtrl.source.landmarks.fn,'landmarks')
      end
      
      menuString = {};
      for landmark_idx = 1:length(landmarks)
        menuString{landmark_idx} = horzcat(landmarks(landmark_idx).frm_id,' - ',landmarks(landmark_idx).type);
      end
      set(hui.landfig.ctrl_panel.landmarksLB,'String',menuString);
      if ~length(menuString) == 1
        set(hui.landfig.ctrl_panel.landmarksLB,'Value',get(hui.landfig.ctrl_panel.landmarksLB,'Value'));
      else
        set(hui.landfig.ctrl_panel.landmarksLB,'Value',1)
      end

      % Create and plot new landmark handle
      x1 = landmarks(insert_idx).rline_start;
      x2 = landmarks(insert_idx).rline_stop;
      y1 = landmarks(insert_idx).rbin_start;
      y2 = landmarks(insert_idx).rbin_stop;
      
      figure(hui.pickfig.handle)
      hold on;
      if ~isempty(hui.pickfig.landmarks_h)
        hui.pickfig.landmarks_h(insert_idx) = plot([x1 x1 x2 x2 x1],[y1 y2 y2 y1 y1],...
          'Color','y','LineStyle','- -','Visible','on');
      else
        hui.pickfig.landmarks_h(insert_idx) = plot([x1 x1 x2 x2 x1],[y1 y2 y2 y1 y1],...
          'Color','y','LineStyle','- -','Visible','on');
      end
      hold off;
      
      debug = 1;
      if debug == 0
        fprintf('Landmark acquisition\nFrame_ID: %s\nRange Bins: %fus to %fus\nRange Lines: %u to %u\nGPS Time Range: %f to %f\nRadar Name: %s\nMission Name: %s\nNotes: %s\nLandmark Type: %s\n',...
          frm_id,rbin_start,rbin_stop,rline_start,rline_stop,gpstime_start,gpstime_stop,radar_name,season_name,notes,lm_type)
      elseif debug == 1
        fprintf('Landmark acquired\n')
      end
      gCtrl.tool.selection.x = NaN;
    end
    
    elseif tool == 11
    % ================================================================
    % Convert layers region tool
    % ================================================================
    if isnan(gCtrl.tool.selection.x)
      gCtrl.tool.selection.x = x;
      gCtrl.tool.selection.y = y;
    else
      cur_layer = get(hui.fig.ctrl_panel.layerPM,'Value');
      % Get layer to be converted
      if cur_layer == 3
        fprintf('  The (c)onvert layer tool only works with layers 1 and 2\n')
        return;
      elseif cur_layer == 1
        convert_layer = 2;
      elseif cur_layer == 2
        convert_layer = 1;  
      end

      % Get selection box region
      [xbnds ybnds outside] = picker_tool_box_sel(gCtrl.tool.selection.x, x, ...
        gCtrl.tool.selection.y, y, axis(hui.pickfig.axes.handle), ...
        length(gCtrl.pick.data));

      % Get valid rangelines where manual or auto points exist that we will convert
      if outside  
        tmp_manual_layer_data = get(hui.pickfig.layer_h(2*convert_layer-1),'YData');
        tmp_auto_layer_data = get(hui.pickfig.layer_h(2*convert_layer),'YData');
        rlines = ceil(xbnds(1)):floor(xbnds(2));
        tmp_manual_ROI = tmp_manual_layer_data(rlines);
        tmp_auto_ROI = tmp_auto_layer_data(rlines);
        rlines_valid_auto = rlines(isfinite(tmp_auto_ROI));
        rlines_valid_manual = rlines(isfinite(tmp_manual_ROI));
      else
        tmp_manual_layer_data = get(hui.pickfig.layer_h(2*convert_layer-1),'YData');
        tmp_auto_layer_data = get(hui.pickfig.layer_h(2*convert_layer),'YData');
        rlines = ceil(xbnds(1)):floor(xbnds(2));
        tmp_manual_ROI = tmp_manual_layer_data(rlines);
        tmp_auto_ROI = tmp_auto_layer_data(rlines);
        rlines_valid_manual = rlines(tmp_manual_ROI>ybnds(1)...
          & tmp_manual_ROI<ybnds(2));
        rlines_valid_auto = rlines(tmp_auto_ROI>ybnds(1)...
          & tmp_auto_ROI<ybnds(2));
      end
      
     if ~isempty(rlines_valid_manual) && isempty(rlines_valid_auto)
        % Convert only manual points then return because no auto layer to
        % convert, then delete old auto points from current layer
        picker_update_layer(cur_layer,1,rlines_valid_manual,...
          tmp_manual_layer_data(rlines_valid_manual),begin_command);
        begin_command = false;
        picker_update_layer(cur_layer,2,rlines,...
          Inf(1,floor(xbnds(2))-ceil(xbnds(1))+1),begin_command);
        if gCtrl.tool.layer_switch 
          % Replace manual points from convert layer with current layer
          % using Shift+c
          begin_command = false;
          picker_update_layer(convert_layer,1,rlines_valid_manual,...
            Inf(1,length(rlines_valid_manual)),begin_command);
          gCtrl.tool.layer_switch = false;
        end

      elseif isempty(rlines_valid_manual) && ~isempty(rlines_valid_auto)
        % Convert only auto points then return because no manual layer to
        % convert, then delete old manual points current layer 
        picker_update_layer(cur_layer,2,rlines_valid_auto,...
          tmp_auto_layer_data(rlines_valid_auto),begin_command);
        begin_command = false;
        picker_update_layer(cur_layer,1,rlines,...
          Inf(1,floor(xbnds(2))-ceil(xbnds(1))+1),begin_command);
        if gCtrl.tool.layer_switch
          % Replace auto points from convert layer with current layer
          % using Shift+c
          begin_command = false;
          picker_update_layer(convert_layer,2,rlines,...
            Inf(1,floor(floor(2))-ceil(xbnds(1))+1),begin_command);
          gCtrl.tool.layer_switch = false;
        end

      else  
        % Both auto and manual points exist, iterate through them both
        picker_update_layer(cur_layer,1,rlines_valid_manual,...
          tmp_manual_layer_data(rlines_valid_manual),begin_command);
        begin_command = false;
        picker_update_layer(cur_layer,2,rlines_valid_auto,...
          tmp_auto_layer_data(rlines_valid_auto),begin_command);
        if gCtrl.tool.layer_switch
          % Replace manual and auto points from convert layer with current layer
          begin_command = false;
          picker_update_layer(convert_layer,1,rlines,...
            Inf(1,floor(xbnds(2))-ceil(xbnds(1))+1),begin_command);
          begin_command = false;
          picker_update_layer(convert_layer,2,rlines,...
            Inf(1,floor(xbnds(2))-ceil(xbnds(1))+1),begin_command);
          gCtrl.tool.layer_switch = false;
        end
      end


      gCtrl.tool.selection.x = NaN;
      frame_modified = true;
    end
  end


  if frame_modified
    % Update frames list and picker title if this frame had not been
    % modified since it was last saved
    if gCtrl.source.modified(gCtrl.source.cur_pick) == ' '
      gCtrl.source.modified(gCtrl.source.cur_pick) = '*';
      menuString = get(hui.fig.ctrl_panel.framesLB,'String');
      menuString{gCtrl.source.cur_pick}(end) = gCtrl.source.modified(gCtrl.source.cur_pick);
      set(hui.fig.ctrl_panel.framesLB,'String',menuString);
      title_str = get(get(hui.pickfig.axes.handle,'Title'),'String');
      title_str(end) = gCtrl.source.modified(gCtrl.source.cur_pick);
      set(get(hui.pickfig.axes.handle,'Title'),'String',title_str);
    end
  end
  
  
elseif but == 3
  % ===================================================================
  % Right mouse button: move cursors
  
  % Find the gps time of click
  cursor_time = interp1(1:length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time), ...
    gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time, x, 'linear','extrap');
  
  % ------------------------------------------------------------------
  % Find the closest index
  [min_dist min_idx] = min(abs(cursor_time - gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time));
  
  % Update the cursor data structures
  gCtrl.pick.cur_idx = min_idx;
  gCtrl.pick.cursor = gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time(gCtrl.pick.cur_idx);
  
  % Update the map cursor
  set(hui.mapfig.pick_cursor.h,'XData', ...
    gCtrl.source.geo{gCtrl.source.cur_pick}.X(gCtrl.pick.cur_idx));
  set(hui.mapfig.pick_cursor.h,'YData', ...
    gCtrl.source.geo{gCtrl.source.cur_pick}.Y(gCtrl.pick.cur_idx));

  % Update the map axis to be centered on this point
  xlims = get(hui.mapfig.axes.handle,'XLim');
  ylims = get(hui.mapfig.axes.handle,'YLim');
  set(hui.mapfig.axes.handle,'XLim',xlims - mean(xlims) + gCtrl.source.geo{gCtrl.source.cur_pick}.X(gCtrl.pick.cur_idx));
  set(hui.mapfig.axes.handle,'YLim',ylims - mean(ylims) + gCtrl.source.geo{gCtrl.source.cur_pick}.Y(gCtrl.pick.cur_idx));

  % Update the pick figure cursor
  if ishandle(hui.pickfig.handle)
    % If the pick window is open...
    cursor_xpos = interp1(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time, ...
      1:length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time), gCtrl.pick.cursor, 'linear','extrap');
    set(hui.pickfig.cursor.h,'XData',[cursor_xpos cursor_xpos]);
    set(hui.pickfig.cursor.h,'YData',[gCtrl.pick.time([1 end])*1e6]);
  end
  
  % Print surface, bed, and thickness in us
  physical_constants;
  rline_idx = min_idx;
  surface_time = gCtrl.source.layers{gCtrl.source.cur_pick}.layerData{1}.value{2}.data(rline_idx);
  bottom_time = gCtrl.source.layers{gCtrl.source.cur_pick}.layerData{2}.value{2}.data(rline_idx);
  thickness_time = bottom_time-surface_time;
  if ~isfinite(surface_time) && ~isfinite(bottom_time)
    % Print nothing    
  elseif isfinite(surface_time) && isfinite(bottom_time)
    fprintf('  Surface: %.2fus, Bottom: %.2fus, Thickness: %.2fus (%.2fm)\n',...
      surface_time*1e6,bottom_time*1e6,thickness_time*1e6,thickness_time*c/2/sqrt(er_ice));
  elseif isfinite(surface_time) && ~isfinite(bottom_time)
    fprintf('  Surface: %.2fus, Bottom: n/a, Thickness: n/a\n',surface_time*1e6);
  elseif ~isfinite(surface_time) && isfinite(bottom_time)
    fprintf('  Surface: n/a, Bottom: %.2fus, Thickness: n/a\n',bottom_time*1e6);
  else
    % Do nothing?
  end
  
end

return
