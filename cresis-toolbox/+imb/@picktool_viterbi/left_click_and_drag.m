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
    
    [manual_idxs,auto_idxs_initial,~] = find_matching_pnts(obj,param,cur_layer);
    
    if ~isempty(auto_idxs_initial)
      
      % Nx: number of along track records/range lines
      Nx = length(image_x);
      
      % Match GT points with axis coordinates
      gt = [interp1(image_x, 1:length(image_x),param.layer.x(manual_idxs), 'nearest', 'extrap');
        interp1(image_y, 1:length(image_y),param.layer.y{cur_layer}(manual_idxs), 'nearest', 'extrap')];
      gt_bound_idxs  = gt(1, 1):gt(1, end);
      x_points       = gt(1, :);
      y_points       = gt(2, :);
      gt             = nan(1, Nx);
      gt(x_points)   = y_points;
      
      % Echogram Parameters
      viterbi_data   = image_c;
      mask           = inf * ones([1 Nx]);
      hori_bounds    = gt_bound_idxs([1 end]);
      vert_bounds    = [];
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
        transition_weight = eval(obj.top_panel.transition_weight_TE.String);
      catch ME
        transition_weight = 1;
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
      
      slope = diff(layers_bins(1,:));
      
      % TODO[reece]: Scale with method Prof. Paden suggested, not based on axis resolutions -- ask for refresher
      transition_weights = ones(1, size(viterbi_data, 2) - 1) * transition_weight;

      if ~obj.top_panel.surf_slope_cbox.Value
        slope(:) = 0;
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
      
      viterbi_timer = tic;
      y_new = tomo.viterbi(double(viterbi_data), double(layers_bins), double(layer_costs), ...
        double(layer_cutoffs), double(mask), double(image_mag_weight), double(slope), ...
        double(max_slope), int64(hori_bounds), int64(vert_bounds), double(mask_dist), double(DIM_costmatrix), ...
        double(transition_weights), double(mult_weight), double(mult_weight_decay), ...
        double(mult_weight_local_decay), int64(zero_bin));
      
      fprintf('Viterbi call took %.2f sec.\n', toc(viterbi_timer));
      
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
        type(manual_idxs - manual_idxs(1) + 1) = 1;
        cmds(end).redo_args = {cur_layer, bounding_idxs, y_new, ...
          type, param.cur_quality*ones(size(bounding_idxs))};
      end
    end
  end
else
  %%% Do nothing for now
end

return
