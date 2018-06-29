function [ stats_obj ] = cluster_kernel_TRWS_3D( sources, refs, detect_param_obj, stats_obj )
%TRWS_RUN Summary of this function goes here
%   Detailed explanation goes here
fprintf('Run TRWS --- Smooth1: %d, Smooth2: %d, GT: %d, smooth3: %d \n', detect_param_obj.smooth_1,...
      detect_param_obj.smooth_2, detect_param_obj.smooth_3);

   len = length(refs);

   for i = 1:len
     clear loaded;
     clear img_data_obj;
     clear ref_data_obj
     fn_1 = sources{i};
     fn_3 = refs{i};

    try
      loaded = load(fn_1);
    catch
      str = 'error1';
      continue;
    end

    img_data_obj.img = loaded.Topography.img;
    img_data_obj.x_size = size(img_data_obj.img, 2);
    img_data_obj.middle_x =  img_data_obj.x_size / 2;
    img_data_obj.middle_x_adjusted = ceil(img_data_obj.middle_x) + 1;
    img_data_obj.range_of_images =  1:size(img_data_obj.img, 3);
    img_data_obj.number_of_images = length(img_data_obj.range_of_images);
%     img_data_obj.mu =  mean(loaded.Topography.mu);
%     img_data_obj.sigma = mean(loaded.Topography.sigma);

    bounds_relative = [3 2 0 0];

    try
      loaded = load(fn_3);
    catch
      str = 'error3';
      continue;
    end
    
    ref_data_obj.this_layer = loaded.surf(2);    % bottom   
    ref_data_obj.active_layer = loaded.surf(ref_data_obj.this_layer.active);
    ref_data_obj.control_layer = loaded.surf(ref_data_obj.this_layer.gt);
    ref_data_obj.surf_layer = loaded.surf(ref_data_obj.this_layer.top);
    ref_data_obj.mask_layer = loaded.surf(ref_data_obj.this_layer.mask);
    start_idx = 1;
 
  labels = [];
  slices_total = 1:size(img_data_obj.img,3);
  
%% TWRS setup

      slice_range = slices_total;
      num_loops = 10;
      threshold = [detect_param_obj.threshold inf]; % check this later
      
      range = Inf;
      polynomial = [];
      mu_size = detect_param_obj.correlation;
      smooth = [detect_param_obj.smooth_1 detect_param_obj.smooth_2 detect_param_obj.smooth_3];
      cols = 1:size(img_data_obj.img,2); % check this later
      
      left_edge_en = 1;
      right_edge_en = 1;
      top_edge_en = 1;
      bottom_edge_en = 1;
      
      
      slice_range = min(slice_range):max(slice_range);
      slices = slice_range;            
      slices = intersect(slices,1:size(img_data_obj.img,3));
      
            
      if ~left_edge_en
        start_slice_idx = 1;
      else
        start_slice_idx = 2;
      end
      
      if ~right_edge_en
        end_slice_idx = length(slices);
      else
        end_slice_idx = length(slices)-1;
      end
 
      gt = [];
      
      if ~isempty(ref_data_obj.control_layer)
        % Create ground truth input
        % 1. Each column is one ground truth input
        % 2. Row 1: relative slice/range-line, Row 2: x, Row 3: y
        for idx = 1:length(slices)
          slice = slices(idx);
          
          mask = isfinite(ref_data_obj.control_layer.x(:,slice)) ...
            & isfinite(ref_data_obj.control_layer.y(:,slice));
          mask(1:bounds_relative(1)) = 0;
          mask(end-bounds_relative(2)+1:end) = 0;
          
          gt = cat(2,gt,[(idx-1)*ones(1,sum(mask)); ...
            ref_data_obj.control_layer.x(mask,slice).'-1; ...
           ref_data_obj.control_layer.y(mask,slice).'+0.5]);
          bottom_bin = ref_data_obj.control_layer.y(33,slices);
        end
      else
        bottom_bin = NaN*zeros(1,length(slices));
      end
      
      if isempty(ref_data_obj.surf_layer)
        %error('trws cannot be run without a surface surface');
        surf_bins = NaN*ref_data_obj.active_layer.y(:,slices);
      else
        surf_bins = ref_data_obj.surf_layer.y(:,slices);
      end
      surf_bins(isnan(surf_bins)) = -1;      
      bottom_bin(isnan(bottom_bin)) = -1;
      
      if isempty(ref_data_obj.mask_layer)
        mask = ones(img_data_obj.x_size,length(slices));
      else
        mask = ref_data_obj.mask_layer.y(:,slices);
      end
      
      begin_slice = max(1, min(slices)-1);
      end_slice = min(size(img_data_obj.img,3), max(slices)+1);
      edge = [ref_data_obj.active_layer.y(:,begin_slice), ref_data_obj.active_layer.y(:,end_slice)];
      edge(mask(:,1) == 0,1) = surf_bins(mask(:,1) == 0,1);
      edge(mask(:,end) == 0,2) = surf_bins(mask(:,end) == 0,end);
      
      trws_data = img_data_obj.img(:,:,slices);
      trws_data(trws_data>threshold(2)) = threshold(2);

     for idx = 1:length(slices)
        for col = 1:size(img_data_obj.img, 2)
          tmp = trws_data(1:min(end,round(surf_bins(col,idx))+70),col,idx);
          tmp(tmp>threshold(1)) = threshold(1);
          trws_data(1:min(end,round(surf_bins(col,idx))+70),col,idx) = tmp;
        end
     end
     
      if ~left_edge_en
        edge(:,1) = -1;
      end
      if ~right_edge_en
        edge(:,end) = -1;
      end
      if ~isempty(polynomial)
        smooth_slope = polyval(polynomial, linspace(-1,1,size(sb.data,2)-1));
      else
        smooth_slope = [];
      end

      smooth_weight = smooth(1:2);
      smooth_var = smooth(3);
      mu = sinc(linspace(-1.5,1.5,mu_size));
      sigma = sum(mu)/20*ones(1,mu_size);
      mask = 90*fir_dec(fir_dec(double(shrink(mask,2)),ones(1,5)/3.7).',ones(1,5)/3.7).';
      mask(mask>=90) = inf;
      bounds = [bounds_relative(1) size(trws_data,2)-bounds_relative(2)-1 -1 -1];
      
      if top_edge_en && ~isempty(cols) && cols(1) > bounds(1)+1
        % Add ground truth from top edge as long as top edge is bounded
        gt = cat(2,gt,[slices-slices(1); ...
        ref_data_obj.active_layer.x(cols(1)-1,slices)-1; ...
        ref_data_obj.active_layer.y(cols(1)-1,slices)+0.5]);
      end
      
      if bottom_edge_en && ~isempty(cols) && cols(end) < bounds(2)
        % Add ground truth from top edge as long as top edge is bounded
        gt = cat(2,gt,[slices-slices(1); ...
          ref_data_obj.active_layer.x(cols(end)+1,slices)-1; ...
          ref_data_obj.active_layer.y(cols(end)+1,slices)+0.5]);
      end
       
      rows = [];
      if isfinite(range) && size(gt,2) > 0
        % Restrict search to range of rows around ground truth
        rows = max(1,min(round(gt(3,:)-range))) : min(size(trws_data,1),max(round(gt(3,:)+range)));
        if length(rows) < length(mu)+1
          error('Error: Range restriction leaves too few rows of data. Increase "Row range" option.');
        end
        trws_data = trws_data(rows,:,:);
        surf_bins = surf_bins - rows(1) + 1;
        bottom_bin = bottom_bin - rows(1) + 1;
        gt(3,:) = gt(3,:) - rows(1) + 1;
      end
      
    
    
     labels = tomo.trws(double(trws_data), ...
        double(surf_bins), double(bottom_bin), double(gt), double(mask), ...
        double(mu), double(sigma), smooth_weight, smooth_var, ...
        double(smooth_slope), double(edge), double(num_loops), int64(bounds));
        
%         figure; imagesc(lp(detect_data)); hold on; plot(labels);        
%      keyboard
     
      bound_front = bounds(1);
      bound_back = bounds(2);
      %%%
%       ref1=ref_data_obj.this_layer.y(:, 1);
%       a =  labels(:,1);
%       abs(a(bound_front:bound_back) - ref1(bound_front:bound_back))'      
      %%%
      
      
      for idx_chk = 1:size(img_data_obj.img,3)
        ref = ref_data_obj.this_layer.y(:, idx_chk);
        label = labels(:,idx_chk);
        layer_diff  = abs(label(bound_front:bound_back) - ref(bound_front:bound_back));
        rmse        = sqrt(mean(abs(layer_diff(:)).^2));
        mean_diff   = nanmean(layer_diff(:));
        median_diff = nanmedian(layer_diff(:));
        min_diff    = nanmin(layer_diff(:));
        max_diff    = nanmax(layer_diff(:));

        stats_obj.rmse_f(stats_obj.counter_f) = rmse;
        stats_obj.diff_f(stats_obj.counter_f) = mean_diff;
        stats_obj.med_f(stats_obj.counter_f)  = median_diff;
        stats_obj.min_df(stats_obj.counter_f) = min_diff;
        stats_obj.max_df(stats_obj.counter_f) = max_diff;
        stats_obj.total_diff  = cat(2, stats_obj.total_diff, layer_diff(:)');
        stats_obj.counter_f   = stats_obj.counter_f + 1;
      end
  
 end
end

