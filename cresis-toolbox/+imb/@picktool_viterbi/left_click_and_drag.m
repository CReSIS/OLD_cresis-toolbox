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
    
    [manual_idxs,auto_idxs,point_idxs] = find_matching_pnts(obj,param,cur_layer);
    
    if length(manual_idxs) < 1
      warning('Insufficient points to track');
      continue;
    elseif ~isempty(auto_idxs)
      
      % Nx: number of along track records/range lines
      Nx = length(image_x);
      custom_data.mu = [11.2575 11.3748 11.4393 11.4555 11.4323   11.3666   11.2668   11.1332   10.9900 10.8484   10.6916];
      custom_data.sigma = [5.4171    5.2945    5.2187    5.1939    5.2174    5.3247    5.4643    5.6571    5.8428 6.0477    6.2935];
      
      % Interpolate surface layer to match image x-axis coordinates
      surf_bins = interp1(param.layer.x,param.layer.y{1},image_x);
      % Interpolate surface layer y-axis units to image pixels
      surf_bins = interp1(image_y, 1:length(image_y),surf_bins);
      % Interpolate all non-finite values using surrounding data
      surf_bins = interp_finite(surf_bins);
      
      % Match GT points with axis coordinates
      gt = [param.layer.x(manual_idxs); interp1(image_y, 1:length(image_y),param.layer.y{cur_layer}(manual_idxs))];
      
      % Echogram Parameters
      viterbi_data   = image_c;
      bottom_bin     = -1;
      mask           = inf * ones([1 Nx]);
      egt_weight     = -1;
      slope          = round(diff(surf_bins));
      bounds         = [];
      viterbi_weight = ones([1 Nx]);
      mu_size        = 31;
      mu             = log10(exp(-(-(mu_size-1)/2 : (mu_size-1)/2).^4/1));
      mu(mu<-30)     = -30;
      mu             = mu - mean(mu);
      sigma          = sum(abs(mu))/10*ones(1,mu_size);
      mc             = -1 * ones(1, Nx);
      mc_weight      = 0;
      CF.sensory_distance = 200;
      CF.max_cost = 50;
      CF.lambda = 0.075;
      
      try
        smooth_weight = str2double(obj.top_panel.smoothness_weight_TE.String);
      catch ME
        smooth_weight = 3;
      end
      
      try
        smooth_var = str2double(obj.top_panel.smoothness_variance_TE.String);
      catch ME
        smooth_var = Inf;
      end
      
      try
        repulsion = str2double(obj.top_panel.repulsion_TE.String);
      catch ME
        repulsion = 150000;
      end
      
      try
        ice_bin_thr = str2double(obj.top_panel.icebinthr_TE.String);
      catch ME
        ice_bin_thr = 10;
      end
      
      %% Detrending
      if 1
        % Along track filtering
        viterbi_data = fir_dec(viterbi_data,ones(1,5)/5,1);
        % Estimate noise level
        noise_value = mean(mean(viterbi_data(end-80:end-60,:)));
        % Estimate trend
        trend = mean(viterbi_data,2);
        trend(trend<noise_value) = noise_value;
        % Subtract trend
        viterbi_data = bsxfun(@minus,viterbi_data,trend);
        % Remove bad circular convolution wrap around at end of record
        viterbi_data(end-70:end,:) = 0;
      end
      
      %% Column restriction between first and last selected GT points
      if obj.top_panel.column_restriction_cbox.Value
        viterbi_data   = viterbi_data(:, auto_idxs);
        surf_bins      = surf_bins(:, auto_idxs);
        mask           = mask(:, auto_idxs);
        viterbi_weight = viterbi_weight(:, auto_idxs);
        mc             = mc(:, auto_idxs);
        slope          = round(diff(surf_bins));
      end
      
      %% Top suppression
      if obj.top_panel.top_sup_cbox.Value
        tic
        topbuffer = 10;
        botbuffer = 30;
        filtvalue = 50;
        for rline = 1 : size(viterbi_data, 2)
          column_chunk = viterbi_data(round(surf_bins(rline) - topbuffer) : ...
            round(surf_bins(rline) + botbuffer), rline);
          viterbi_data(round(surf_bins(rline) - topbuffer) : ...
            round(surf_bins(rline) + botbuffer), rline) = imgaussfilt(column_chunk, filtvalue);
        end
        fprintf('Top suppression took %.2f sec.\n', toc);
      end
      
      %% Multiple suppression
      if obj.top_panel.mult_sup_cbox.Value
        tic
        topbuffer = 10;
        botbuffer = 5;
        filtvalue = 50;
        for rline = 1 : size(viterbi_data, 2)
          column_chunk = viterbi_data(round(2*surf_bins(rline) - topbuffer) : ...
            round(2*surf_bins(rline) + botbuffer), rline);
          viterbi_data(round(2*surf_bins(rline) - topbuffer) : ...
            round(2*surf_bins(rline) + botbuffer), rline) = imgaussfilt(column_chunk, filtvalue);
        end
        fprintf('Multiple suppression took %.2f sec.\n', toc);
      end
      
      %% Call viterbi.cpp
      tic
      [y_new, cost] = tomo.viterbi(double(viterbi_data), double(surf_bins), ...
        double(bottom_bin), double(gt), double(mask), double(mu), ...
        double(sigma), double(egt_weight), double(smooth_weight), ...
        double(smooth_var), double(slope), int64(bounds), ...
        double(viterbi_weight), double(repulsion), double(ice_bin_thr), ...
        double(mc), double(mc_weight), ...
        double(CF.sensory_distance), double(CF.max_cost), double(CF.lambda));
      fprintf('Viterbi call took %.2f sec.\n', toc);

      if obj.top_panel.column_restriction_cbox.Value
        y_new(end) = y_new(end-1);
      else
        y_new = y_new(auto_idxs);
      end
      
      % Interpolate layer to match image y-axis
      y_new  = interp1(1:length(image_y), image_y, y_new);
      cmds(end+1).undo_cmd = 'insert';
      % Quality measurement from Viterbi algorithm result
      if obj.top_panel.quality_output_cbox.Value
        try
          thrs = str2double(obj.top_panel.quality_threshold_TE.String);
        catch ME
          thrs = -20;
        end        
        quality = ones(size(cost));
        quality(cost < thrs) = 3;
        cmds(end).undo_args = {cur_layer, auto_idxs, ...
        param.layer.y{cur_layer}(auto_idxs), ...
        param.layer.type{cur_layer}(auto_idxs), quality};
      else
        cmds(end).undo_args = {cur_layer, auto_idxs, ...
        param.layer.y{cur_layer}(auto_idxs), ...
        param.layer.type{cur_layer}(auto_idxs), ...
        param.layer.qual{cur_layer}(auto_idxs)};
      end

      cmds(end).redo_cmd = 'insert';
      if obj.top_panel.quality_output_cbox.Value
        cmds(end).redo_args = {cur_layer, auto_idxs, y_new, ...
          2*ones(size(auto_idxs)), quality};        
      else
        cmds(end).redo_args = {cur_layer, auto_idxs, y_new, ...
          2*ones(size(auto_idxs)), param.cur_quality*ones(size(auto_idxs))};
      end
    end
  end
else
  %%% Do nothing for now
end

return
