function  [ stats_obj ]  = cluster_kernel_Viterbi_3D(sources, refs, detect_param_obj, stats_obj)

  %% load the frame and ground truth (use a try catch block later)  
 fprintf('Run Viterbi --- SmoothWeight: %d, SmoothVariance: %d, GT: %d, Repulsion: %d \n', detect_param_obj.smooth_weight,...
      detect_param_obj.smooth_var, detect_param_obj.egt_weight, detect_param_obj.repulsion);
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
%   img_data_obj.mu =  mean(loaded.Topography.mu);
%   img_data_obj.sigma = mean(loaded.Topography.sigma);
  
  bounds_relative = [3 2 0 0]; % only for 2014 Greenland P3

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
  
  %% processing
%   if detect_param_obj.previous % if previous flag is on
%     start_idx = 2;
%     detect_param_obj.viterbi_weight = 2;
%   else
%     start_idx = 1;
%     detect_param_obj.viterbi_weight = 1;
%   end

  start_idx = 1;
  labels = [];
  slices = 1:size(img_data_obj.img,3);
  
  for idx = start_idx : length(slices)
    % do set_ground_truth
    slice = slices(idx);      
    slice_range = 3;
    [gt, bottom_bin, surf_bins, mask, viterbi_weight] = set_ground_truth(slice, bounds_relative, ...
      detect_param_obj.previous, img_data_obj, ref_data_obj, idx, labels, slices, slice_range);
    
    % check if viterbi weight is the same
    if(1 ~= unique(viterbi_weight))
      disp(detect_param_obj.viterbi_weight);
      disp(viterbi_weight);
      error('different viterbi weight\n');
    end
    
    % setup detect_data
    detect_data = set_detect_data(img_data_obj, detect_param_obj, idx);
    bounds = [bounds_relative(1) size(detect_data,2)-bounds_relative(2)-1];
    mc            = -1 * ones(1, size(detect_data,2));
    mc_weight     = 0;
    mu_size = 11;
    mu = sinc(linspace(-1.5, 1.5, mu_size));
    sigma         = sum(mu)/20*ones(1,mu_size);
    
    mask           = 90*fir_dec(fir_dec(double(shrink(mask,2)),ones(1,5)/3.7).',ones(1,5)/3.7).';
    mask(mask>=90) = inf;
    
    if(slice_range+1 > size(mask,2))      
       mask = mask(:,slice_range);
    else
       mask = mask(:,slice_range+1);
    end
         
    
     labels = tomo.viterbi(double(detect_data), ...
          double(surf_bins), double(bottom_bin), ...
          double(gt), double(mask), ...
          double(mu), ...
          double(sigma), ...
          double(detect_param_obj.egt_weight), ...
          double(detect_param_obj.smooth_weight), ... 
          double(detect_param_obj.smooth_var), ...
          double(detect_param_obj.slope), ...
          int64(bounds), double(viterbi_weight), ...
          double(detect_param_obj.repulsion), ...
          double(detect_param_obj.ice_b_threshold), ...
          double(mc), double(mc_weight)); 
     
      %    keyboard   
      %    figure; imagesc(lp(detect_data)); hold on; plot(labels);

      bound_front = bounds(1);
      bound_back = bounds(2);
      ref = ref_data_obj.this_layer.y(:, idx);
      layer_diff  = abs(labels(bound_front:bound_back).' - ref(bound_front:bound_back));
      rmse        = sqrt(mean(layer_diff(:).^2));
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
      stats_obj.counter_f = stats_obj.counter_f + 1;

  end
 end
 
end

