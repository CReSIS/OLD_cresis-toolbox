function [ gt, bottom_bin, surf_bins, mask, viterbi_weight] = set_ground_truth( slice, bounds_relative, previous, img, ref, idx, labels, slices, slice_range)
%SET_GROUND_TRUTH 
%
% Input:
%  previous: logical; 1 or 0
%  img: img_Data object
%  ref : ref_Data object
%  idx: current index
%  labels: previous detect result

  % If previous is turned on, then ground truth is set to the detected
  % result of the previous slide
  if previous
    slice_prev = idx-1;
    if idx == 2
      % if index == 2; we set it to the hand picked ones (in ref_Data object)
      gt = [ref.active_layer.x(:,slice_prev).'-1; ...
        ref.active_layer.y(:,slice_prev).'+0.5];
    else    
      gt = [ref.active_layer.x(:,slice_prev).'-1; ...
                labels(:).'+0.5];
    end
  else
    gt = [];
  end
  
  % if hand picked data in this slide is presented, we use them
  if ~isempty(ref.control_layer)
     mask = isfinite(ref.control_layer.x(:,slice)) ...
      & isfinite(ref.control_layer.y(:,slice));
     % from mask, we got 
     gt_this_slide = [ref.control_layer.x(mask, slice).'-1; ...
       ref.control_layer.y(mask, slice).'+0.5];
     % concatenate the result (will be 65 columns)
     gt = cat(2, gt, gt_this_slide);

     % replace columns in gt that has index overlaps with gt_this_slide 
    [~,unique_indices] = unique(gt(1,:),'last','legacy');
    gt = gt(:,unique_indices);
    
    viterbi_weight = ones([1 (size(img.img,2))]);        
    

    % Victor said changed to 2 
    viterbi_weight(1+gt(1,:)) = 2;
 
    [~,sort_indices] = sort(gt(1,:)); % check if this is redundant
    gt = gt(:,sort_indices);

    % setup bottom bin    
    bottom_bin = ref.control_layer.y(img.middle_x_adjusted, slice);%.'+0.5;
  else
    bottom_bin = NaN; % why?
    viterbi_weight = ones([1 length(gt)]);
  end
  
   
  if ~isempty(ref.control_layer.gt)
    % Create ground truth input
    % 1. Each column is one ground truth input
    % 2. Row 1: relative slice/range-line, Row 2: x, Row 3: y
    slice_range = 3;
    slices = slice - slice_range:slice + slice_range;
    slices = intersect(slices, 1:size(img.img,3));
    
    
    for idx = 1:length(slices)
      slice = slices(idx);

      mask = isfinite(ref.control_layer.x(:,slice)) ...
        & isfinite(ref.control_layer.y(:,slice));
      mask(1:bounds_relative(1)) = 0;
      mask(end-bounds_relative(2)+1:end) = 0;
    end    
  end
  
  
  if isempty(ref.control_layer.mask)
    mask = ones(size(img.img,2), length(slices));
  else
    mask = ref.mask_layer.y(:,slices);
  end

  
  if isempty(ref.surf_layer)
    surf_bins = NaN * ref.get_active_y_at(slice);
  else
    surf_bins = ref.surf_layer.y(:,slice);
  end  
  
  surf_bins(isnan(surf_bins)) = -1;
  bottom_bin(isnan(bottom_bin)) = -1;
end

