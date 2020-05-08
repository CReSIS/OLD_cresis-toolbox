function cmds = left_click_and_drag(obj,param)
% cmds = left_click_and_drag(obj,param)
%
% Detect tool
%
% Compile with
%   mex -largeArrayDims viterbi.cpp

image_x = param.image_x;
image_y = param.image_y;
image_c = param.image_c;
cur_layers = param.cur_layers;
x = param.x;
y = param.y;
cmds = [];

fprintf('Performing viterbi tracking on points %f to %f, %f to %f\n', x, y);

param.x_bounds = 3;
param.y_bounds = 1;

tool_idx = get(obj.top_panel.tool_PM,'Value');
if tool_idx == 1
  
  %=========================================================================
  
  for layer_idx = 1:length(cur_layers)
    cur_layer = cur_layers(layer_idx);
    
    [selected_manual_idxs,~,~] = find_matching_pnts(obj,param,cur_layer);
    param.x_bounds = 1;
    [all_manual_idxs,~,~] = find_matching_pnts(obj,param,cur_layer);
    
    if true
      
      % Nx: number of along track records/range lines
      Nx = length(image_x);
      
      % Match GT points with axis coordinates
      gt = [interp1(image_x, 1:length(image_x),param.layer.x(selected_manual_idxs), 'nearest', 'extrap');
        interp1(image_y, 1:length(image_y),param.layer.y{cur_layer}(selected_manual_idxs), 'nearest', 'extrap')];
      upper_bound = interp1(image_y, 1:length(image_y),param.y(1), 'nearest', 'extrap');
      lower_bound = interp1(image_y, 1:length(image_y),param.y(2), 'nearest', 'extrap');
      x_points       = gt(1, :);
      y_points       = gt(2, :);
      gt             = nan(1, Nx);
      gt(x_points)   = y_points;
      
      % Echogram Parameters
      viterbi_data   = image_c;
      mask           = inf * ones([1 Nx]);
      mask_dist      = round(bwdist(mask == 0));
      
      %% Detrending
      if 1
        % viterbi_data(end-140:end,:) = 0;
        viterbi_data = echo_norm(viterbi_data,struct('scale',[-40 90]));
        % % Along track filtering
        % viterbi_data = fir_dec(viterbi_data,ones(1,5)/5,1);
        % % Estimate noise level
        % noise_value = mean(mean(viterbi_data(end-80:end-60,:)));
        % % Estimate trend
        % trend = mean(viterbi_data,2);
        % trend(trend<noise_value) = noise_value;
        % % Subtract trend
        % viterbi_data = bsxfun(@minus,viterbi_data,trend);
        % % Remove bad circular convolution wrap around at end of record
        % viterbi_data(end-70:end,:) = 0;
      end
      
      dt = param.echo_time(2) - param.echo_time(1);
      zero_bin = floor(-param.echo_time(1)/dt + 1);
      
      %% Distance-to-Ice-Margin model
      clear DIM DIM_costmatrix;
      global gRadar
      DIM = load(fullfile(gRadar.path, '+tomo', 'Layer_tracking_2D_parameters_Matrix.mat'));
      DIM_costmatrix = DIM.Layer_tracking_2D_parameters;
      DIM_costmatrix = DIM_costmatrix .* (200 ./ max(DIM_costmatrix(:)));
      
      % Surface and multiple suppression weights
      try
        layers = eval(obj.top_panel.layers_TE.String);
      catch
        layers = [1];
      end
      try
        layers_weight = eval(obj.top_panel.layers_weight_TE.String);
      catch
        layers_weight = [1000];
      end
      
      try
        mult_weight = eval(obj.top_panel.mult_weight_TE.String);
      catch ME
        mult_weight = 100;
      end
      try
        mult_weight_decay = eval(obj.top_panel.mult_weight_decay_TE.String);
      catch ME
        mult_weight_decay = 0;
      end
      try
        mult_weight_local_decay = eval(obj.top_panel.mult_weight_local_decay_TE.String);
      catch ME
        mult_weight_local_decay = 0.8;
      end
      
      try
        max_slope = eval(obj.top_panel.max_slope_TE.String);
      catch ME
        max_slope = -1;
      end
      try
        along_track_weight = eval(obj.top_panel.transition_weight_TE.String);
      catch ME
        along_track_weight = 1;
      end
      try
        image_mag_weight = eval(obj.top_panel.image_mag_weight_TE.String);
      catch ME
        image_mag_weight = 1;
      end
      try
        gt_weight = -eval(obj.top_panel.ground_truth_weight_TE.String);
      catch ME
        gt_weight = -1;
      end
      try
        gt_cutoff = eval(obj.top_panel.ground_truth_cutoff_TE.String);
      catch ME
        gt_cutoff = 5;
      end
      
      if obj.top_panel.r_echo.Value
        hori_bounds = [1 Nx];
      end
      if obj.top_panel.r_sel.Value
        hori_bounds = round(param.x([1 2]));
      end
      if obj.top_panel.r_extr.Value
        selected_manual_idxs = find(~isnan(gt)&gt>=upper_bound&gt<=lower_bound);
        if length(selected_manual_idxs) < 2
          error('Less than 2 gt points are selected but horizontal bounding option is set to extreme gt points.');
        end
        hori_bounds = selected_manual_idxs([1 end]);
      end
      bound_gt_idxs = find(~isnan(gt));
      bound_gt_idxs = bound_gt_idxs(bound_gt_idxs >= hori_bounds(1) & bound_gt_idxs <= hori_bounds(2));


      layers_bins = zeros(length(layers),length(image_x));
      layer_costs = zeros(length(layers),length(image_x));
      for layers_idx = 1:length(layers)
        % Interpolate surface layer to match image x-axis coordinates
        new_layer_bins = interp1(param.layer.x,param.layer.y{layers(layers_idx)},image_x);
        % Interpolate surface layer y-axis units to image pixels
        new_layer_bins = interp1(image_y, 1:length(image_y),new_layer_bins);
        % Interpolate all non-finite values using surrounding data
        new_layer_bins = interp_finite(new_layer_bins, 0);
        layers_bins(layers_idx,:) = new_layer_bins;
        layer_costs(layers_idx,:) = layers_weight(layers_idx);
      end
      
      along_track_slope = diff(layers_bins(1,:));
      
      % TODO[reece]: Scale with method Prof. Paden suggested, not based on axis resolutions -- ask for refresher
      transition_weights = ones(1, size(viterbi_data, 2) - 1) * along_track_weight;

      if ~obj.top_panel.surf_slope_cbox.Value
        along_track_slope(:) = 0;
      end

      layers_bins = [layers_bins; gt];
      
      gt_costs = nan(1, size(viterbi_data, 2));
      gt_costs(~isnan(gt)) = gt_weight;
      layer_costs = [
        layer_costs;
        gt_costs
      ];
      
      gt_cutoffs = nan(1, size(viterbi_data, 2));
      gt_cutoffs(~isnan(gt)) = gt_cutoff;
      layer_cutoffs = [
        nan(length(layers), size(viterbi_data, 2));
        gt_cutoffs
      ];

      upper_bounds = ones(1, Nx)*upper_bound;
      lower_bounds = ones(1, Nx)*lower_bound;

      gt_idxs = find(~isnan(gt(1, :)));

      upper_bounds(gt_idxs) = max([gt(gt_idxs) - gt_cutoff; upper_bounds(gt_idxs)]);
      lower_bounds(gt_idxs) = min([gt(gt_idxs) + gt_cutoff; lower_bounds(gt_idxs)]);

      upper_bounds = max([upper_bounds; layers_bins(1, :)]);

      viterbi_timer = tic;
      [y_new, debug] = tomo.viterbi2(double(viterbi_data), along_track_slope, along_track_weight, upper_bounds, lower_bounds);
      fprintf('Viterbi call took %.2f sec.\n', toc(viterbi_timer));
      
      figure('NumberTitle', 'off', 'Name', 'Suppression On');
      imagesc(debug);
      hold on;
      plot(y_new, 'r');
      hold off;

      bounding_idxs = hori_bounds(1):hori_bounds(2);
      y_new = y_new(bounding_idxs);
      % Resample and interpolate y_new to match layer axes
      y_new = interp1(1:length(image_y), image_y,y_new,'linear','extrap');
      y_new = interp1(image_x(bounding_idxs),y_new,param.layer.x(bounding_idxs),'linear', 'extrap');
      
      cmds(end+1).undo_cmd = 'insert';
      % Quality measurement from Viterbi algorithm result
      if false % obj.top_panel.quality_output_cbox.Value
        % quality calculation not implemented
        try
          thrs = str2double(obj.top_panel.quality_threshold_TE.String);
        catch ME
          thrs = -20;
        end
        quality = ones(size(cost));
        quality(cost < thrs) = 3;
        cmds(end).undo_args = {cur_layer, bounding_idxs, ...
          param.layer.y{cur_layer}(bounding_idxs), ...
          param.layer.type{cur_layer}(bounding_idxs), quality};
      else
        cmds(end).undo_args = {cur_layer, bounding_idxs, ...
          param.layer.y{cur_layer}(bounding_idxs), ...
          param.layer.type{cur_layer}(bounding_idxs), ...
          param.layer.qual{cur_layer}(bounding_idxs)};
      end
      
      cmds(end).redo_cmd = 'insert';
      if false % obj.top_panel.quality_output_cbox.Value
        cmds(end).redo_args = {cur_layer, bounding_idxs, y_new, ...
          2*ones(size(bounding_idxs)), quality};
      else
        type = 2*ones(size(bounding_idxs));
        if obj.top_panel.r_echo.Value && ~isempty(all_manual_idxs)
            type(all_manual_idxs - bounding_idxs(1) + 1) = 1;
        elseif ~isempty(bound_gt_idxs) % Only affects manual points within selection for other options
            type(bound_gt_idxs - bounding_idxs(1) + 1) = 1;
        end
        cmds(end).redo_args = {cur_layer, bounding_idxs, y_new, ...
        type, param.cur_quality*ones(size(bounding_idxs))};
      end
    end
  end
else
  %%% Do nothing for now
end

return
