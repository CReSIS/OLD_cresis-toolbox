function labels = tracker_viterbi(data,track)
% labels = tracker_viterbi(data,track)
%

gt = track.crossovers;
physical_constants;

try
  transition_weight = track.viterbi.transition_weight;
catch ME
  transition_weight = 1;
  warning('transition_weight not specified or invalid. Defaulting to %.2f.', ...
          transition_weight);
end
try
  gt_cutoff = track.viterbi.gt_cutoff;
catch ME
  gt_cutoff = 50;
  warning('gt_cutoff not specified or invalid. Defaulting to %.0f.', ...
          gt_cutoff);
end

along_track_slope = zeros(1, size(data, 2) - 1);

upper_bounds = track.min_bin;
lower_bounds = track.max_bin;
upper_bounds(gt(1, :)) = gt(2, :) - gt_cutoff;
lower_bounds(gt(1, :)) = gt(2, :) + gt_cutoff;

upper_bounds(upper_bounds < 1) = 1;
lower_bounds(lower_bounds > size(data, 1)) = size(data, 1);

labels = tomo.viterbi2(single(data), along_track_slope, transition_weight, ...
                       upper_bounds, lower_bounds);
