function labels = tracker_viterbi(data,track)
% labels = tracker_viterbi(data,track)
%

gt = track.crossovers;
physical_constants;

try
  transition_weight = track.viterbi.transition_weight;
catch ME
  transition_weight = 1;
  warning("transition_weight not specified or invalid. Using default of %d.", transition_weight);
end
try
  gt_cutoff = track.viterbi.gt_cutoff;
catch ME
  gt_cutoff = 50;
  warning("gt_cutoff not specified or invalid. Using default of %d.", gt_cutoff);
end

along_track_slope = zeros(1, size(data, 2) - 1);

upper_bounds = nan(1, size(data, 2));
lower_bounds = nan(1, size(data, 2));
upper_bounds(gt(1, :)) = gt(2, :) - gt_cutoff;
lower_bounds(gt(1, :)) = gt(2, :) + gt_cutoff;

labels = tomo.viterbi2(double(data), along_track_slope, transition_weight, upper_bounds, lower_bounds);
