function idxs_box = boxing_1D(idx, box_bounds_each_side, box_bounds_max)

idx_box_min = max([idx - box_bounds_each_side, 1]);
idx_box_max = min([idx + box_bounds_each_side, box_bounds_max]);
idxs_box = idx_box_min:idx_box_max;
