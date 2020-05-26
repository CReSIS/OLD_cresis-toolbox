function labels = tracker_viterbi(data,track)
% labels = tracker_viterbi(data,track)
%

gt = track.crossovers;

try
  transition_weight = track.viterbi.transition_weight;
catch ME
  transition_weight = 1;
end
try
  gt_cutoff = track.viterbi.gt_cutoff;
catch ME
  gt_cutoff = 50;
end

elevation = param.echowin.eg.elev;
vel_air = c/2;
vel_ice = c/(sqrt(er_ice)*2);
dt = param.echo_time(2) - param.echo_time(1);
along_track_slope = diff(elevation);

yaxis_choice = get(param.echowin.left_panel.yaxisPM,'Value');
if yaxis_choice == 1 % TWTT
  drange = dt * vel_air;
elseif yaxis_choice == 2 % WGS_84 Elevation
  drange = dt * vel_ice;
elseif yaxis_choice == 3 % Range
  drange = dt * vel_ice;
elseif yaxis_choice == 4 % Range bin
  drange = dt * vel_air;
elseif yaxis_choice == 5 % Surface flat
  drange = dt * vel_ice;
end
along_track_slope = round(along_track_slope / drange);

%% Call viterbi.cpp
upper_bounds = nan(1, size(data, 2));
lower_bounds = nan(1, size(data, 2));
upper_bounds(gt(1, :)) = gt(2, :) - gt_cutoff;
lower_bounds(gt(1, :)) = gt(2, :) + gt_cutoff;

labels = tomo.viterbi(double(data), along_track_slope, transition_weight, upper_bounds, lower_bounds);
