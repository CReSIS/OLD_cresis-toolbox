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
      mask_dist      = round(bwdist(mask == 0));
      
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
        mask_dist      = round(bwdist(mask == 0));
        viterbi_weight = viterbi_weight(:, auto_idxs);
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
      
      %% Distance-to-Ice-Margin model
      % Obtained from geostatistical analysis of 2014 Greenland P3
      prob.DIM_means = [    16.0582   27.3381   34.0492   40.5714...
        46.9463   52.3826   58.4664   63.7750   70.5143   76.3140...
        81.3519   86.9523   92.2285   98.1430   102.8310  107.0558...
        112.4538  116.3923  121.0109  125.1486  128.7611  133.3286...
        136.3500  139.5058  142.3812  146.1565  148.3781  151.1398...
        153.5642  155.2606  157.4824  159.8529  161.0239  163.1799...
        164.2849  166.1013  166.0728  167.3132  168.1448  169.1323...
        169.7559  170.3869  171.4264  171.8926  171.5201  171.8870...
        171.6407  172.3505  171.3533  171.6161];
      
      prob.DIM_vars  = [     20.5619   24.9082   28.0037   32.0840...
        35.6021   38.9544    42.4034   45.3588   48.9714   52.1360...
        55.0826   57.5144    60.3847   63.1485   65.1199   67.3734...
        69.3662   71.2849    72.9471   74.3759   75.5521   76.9737...
        77.9961   79.3596    79.9999   81.0342   81.6340   82.2424...
        82.9658   83.4794    84.1632   84.4168   85.0014   85.3065...
        85.7757   86.1880    86.3563   86.5577   86.7289   87.0748...
        87.1360   87.2473    87.2828   86.6350   86.5453   86.2989...
        86.4736   86.7318   87.1606   87.6966];
      
      costm = ones(1201, 101);
      fd = 18 * fir_dec(1:101, ones(1,5)/3.7.');
      
      for t = 1:1200
        for i = 1:101
          costm(t,i) = 200 * exp(-0.075 .* t) - 200 * exp(-0.075 * 50);
          if t > fd(i)
            costm(t,i) = 200;
          end
        end
      end
      
      prob.DIM_means = [prob.DIM_means prob.DIM_means(end)* ones(1,50)];
      prob.DIM_vars  = [prob.DIM_vars  prob.DIM_vars(end) * ones(1,50)];
      
      DIM_costmatrix = ones(1200, 100);
      for DIM = 1 : 100
        for T = 1 : 1200
          var = DIM * 0.05 * prob.DIM_vars(end);
          cost = 0.1 ./ normpdf(T, prob.DIM_means(DIM), var);
          cost(cost > 800) = 800;
          DIM_costmatrix(T, DIM) = cost;
        end
      end
      
      for k = 1:100
        DIM_costmatrix(1:50, k) = DIM_costmatrix(1:50, k) + 0.06 * k * costm(1:50, k);
      end
      
      DIM_costmatrix(DIM_costmatrix < 0) = 0;
      DIM_costmatrix(DIM_costmatrix > 800) = 800;
      
      DIM_costmatrix = DIM_costmatrix(3:end, :);
      DIM_costmatrix(end+1, :) = DIM_costmatrix(end,:);
      DIM_costmatrix(end+1, :) = DIM_costmatrix(end,:);
      DIM_costmatrix = DIM_costmatrix ./ 4;
      
      % Visualization of DIM_costmatrix
      if 0
        figure; (imagesc(DIM_costmatrix)); colorbar; hold on;
        xlabel('Distance to nearest ice-margin [m]');
        ylabel('Ice thickness distribution [range-bin]')
        title('Added cost');
      end

      transition_mu = -1 * ones(1, size(viterbi_data, 2));
      transition_sigma = -1 * ones(1, size(viterbi_data, 2));

      tic
      y_new = tomo.viterbi(double(viterbi_data), double(surf_bins), ...
        double(bottom_bin), double(gt), double(mask), double(mu), ...
        double(sigma), double(egt_weight), double(smooth_weight), ...
        double(smooth_var), double(slope), int64(bounds), ...
        double(viterbi_weight), double(repulsion), double(ice_bin_thr), ...
        double(mask_dist), double(DIM_costmatrix), ...
        double(transition_mu), double(transition_sigma));
      toc
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
