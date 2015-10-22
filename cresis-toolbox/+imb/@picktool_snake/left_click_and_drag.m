function cmds = left_click_and_drag(obj,param)
% cmds = left_click_and_drag(obj,param)
%
% Snake tool

image_x = param.image_x;
image_y = param.image_y;
image_c = param.image_c;
cur_layers = param.cur_layers;
x = param.x;
y = param.y;
cmds = [];

fprintf('Snake points %f to %f, %f to %f\n', x, y);

param.x_bounds = 3;
param.y_bounds = 1;

%% Get snake method
tool_idx = get(obj.top_panel.tool_PM,'Value');
if tool_idx == 1
  
  %=========================================================================
  %% Basic Snake
  for layer_idx = 1:length(cur_layers)
    cur_layer = cur_layers(layer_idx);
    
    [manual_idxs,auto_idxs,point_idxs] = find_matching_pnts(obj,param,cur_layer);
    
    if length(manual_idxs) < 1
      warning('Insufficient points to snake');
      continue;
    elseif ~isempty(auto_idxs)
      % Run snake on values
      [y_new] = obj.snake(image_c,image_x,image_y,param.layer.x(manual_idxs), ...
        param.layer.y{cur_layer}(manual_idxs),param.layer.x(auto_idxs));
      
      cmds(end+1).undo_cmd = 'insert';
      cmds(end).undo_args = {cur_layer, auto_idxs, ...
        param.layer.y{cur_layer}(auto_idxs), ...
        param.layer.type{cur_layer}(auto_idxs), ...
        param.layer.qual{cur_layer}(auto_idxs)};
      cmds(end).redo_cmd = 'insert';
      cmds(end).redo_args = {cur_layer, auto_idxs, ...
        y_new, ...
        2*ones(size(auto_idxs)), param.cur_quality*ones(size(auto_idxs))};
    end
  end
  
elseif tool_idx == 2
  %% ==========================================================================
  % Crandall
  
  % setup params
  top_smooth = 10^(get(obj.bottom_panel.topSmoothSlider, 'Value')*6);
  top_peak = get(obj.bottom_panel.topPeakRatioSlider, 'Value');
  bottom_smooth = 10^(get(obj.bottom_panel.bottomSmoothSlider, 'Value')*6);
  bottom_peak = get(obj.bottom_panel.bottomPeakRatioSlider, 'Value');
  repulse = 10^(get(obj.bottom_panel.repulseSlider, 'Value')*2);
  opts = [top_smooth bottom_smooth top_peak bottom_peak repulse];
  
  for layer_idx = 1:length(cur_layers)
    cur_layer = cur_layers(layer_idx);
    x_auto = param.layer_x_auto{layer_idx};
    y_auto = param.layer_y_auto{layer_idx};
    x_manual = param.layer_x_manu{layer_idx};
    y_manual = param.layer_y_manu{layer_idx};
    manual_auto = layer_type{cur_layer};
    
    % combine manual and auto points
    [layer_data_x,I] = sort([x_manual,x_auto]);
    y_cat = [y_manual,y_auto];
    layer_data_y = y_cat(I);
    % Determine which points are in the rubberband/box
    dx = x_data(2) - x_data(1);
    x_idx_1 = find(layer_data_x >= x(1)-dx/2,1,'first');
    x_idx_2 = find(layer_data_x <= x(2)+dx/2,1,'last');
    old_idx = x_idx_1:x_idx_2;
    manual_auto = manual_auto(old_idx);
    old_idx_manual = old_idx(logical(manual_auto-2));
    old_idx_auto = old_idx(logical(manual_auto-1));
    y_old_manual = layer_data_y(old_idx_manual); % in current unit
    x_old_manual = layer_data_x(old_idx_manual); % in current unit
    A_indices = x_data >= (x(1)-dx/2) & x_data <= (x(2)+dx/2);
    A = A(:,A_indices);
    if ~isfield(param,'point') || isempty(param.point.x)
      pnts_y = interp1(y_data,1:length(y_data),...
        y_old_manual,'linear','extrap');
      pnts_x = interp1(x_data,1:length(x_data),...
        x_old_manual,'linear','extrap');
      pnts = [round(pnts_x);round(pnts_y)];
      if length(pnts_x) < 1  % no manual points
        [AA path]=stereo(1, double(A), opts);
      elseif cur_layer == 1  % have manual points && surface
        [AA path]=stereo(1, double(A), opts,double(pnts),[]);
      elseif cur_layer == 2  % have manual points && bottom
        [AA path]=stereo(1, double(A), opts,[],double(pnts));
      end
      y_new = interp1(1:length(y_data),y_data,path(cur_layer,:),'linear','extrap');
      x_new = x_data(A_indices);
      % Run snake on values
      cmds(end+1).cmd = 'snake crandall';
      cmds(end).begin_cmd = true;
      cmds(end).cur_layer = cur_layer;
      cmds(end).type = manual_auto; % auto
      cmds(end).qual = param.cur_qual;
      cmds(end).old_idx = old_idx;
      cmds(end).new_pnts.x = x_new;
      cmds(end).new_pnts.y = y_new;
    else
      x_old_manual(end+1) = param.point.x;
      y_old_manual(end+1) = param.point.y;
      [x_old_manual sort_idx] = sort(x_old_manual);
      y_old_manual = y_old_manual(sort_idx);
      pnts_y = interp1(y_data,1:length(y_data),...
        y_old_manual,'linear','extrap');
      pnts_x = interp1(x_data,1:length(x_data),...
        x_old_manual,'linear','extrap');
      pnts = [round(pnts_x);round(pnts_y)];
      if length(pnts_x) < 1  % no manual points
        [AA path]=stereo(1, double(A), opts);
      elseif cur_layer == 1  % have manual points && surface
        [AA path]=stereo(1, double(A), opts,double(pnts),[]);
      elseif cur_layer == 2  % have manual points && bottom
        [AA path]=stereo(1, double(A), opts,[],double(pnts));
      end
      y_new = interp1(1:length(y_data),y_data,path(cur_layer,:),'linear','extrap');
      x_new = x_data(A_indices);
      %-------first delete auto points
      cmds(end+1).cmd = 'delete auto pnts';
      cmds(end).begin_cmd = true;
      cmds(end).cur_layer = cur_layer;
      cmds(end).type = 2; % auto
      cmds(end).qual = param.cur_qual;
      cmds(end).old_idx = old_idx;
      cmds(end).new_pnts.x = [];
      cmds(end).new_pnts.y = [];
      %-------then add the single manual point
      cmds(end+1).cmd = 'insert';
      cmds(end).begin_cmd = false;
      cmds(end).cur_layer = cur_layer;
      cmds(end).type = 1; % manual
      cmds(end).qual = param.cur_qual;
      cmds(end).old_idx = [];
      cmds(end).new_pnts.x = param.point.x;
      cmds(end).new_pnts.y = param.point.y;
      %-------last add derived points
      cmds(end+1).cmd = 'snake Crandall';
      cmds(end).begin_cmd = false;
      cmds(end).cur_layer = cur_layer;
      cmds(end).type = 2; % auto
      cmds(end).qual = param.cur_qual;
      cmds(end).old_idx = [];
      cmds(end).new_pnts.x = x_new;
      cmds(end).new_pnts.y = y_new;
    end
    if get(obj.top_panel.reinterp_mode_cbox,'Value') == 1
      
      obj.last_tool = tool_idx;
      obj.last_layers = cur_layers;
      obj.last_range_gps = interp1(x_data,param.image_gps_time,...
        x,'linear','extrap');
    end
    
  end
  
  
  
  
  
  
  
elseif tool_idx == 3;
  %%==========================================================================
  % Panton
  
  
end

return

