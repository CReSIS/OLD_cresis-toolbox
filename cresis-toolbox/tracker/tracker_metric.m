% tracker_metric
%
% Function to compare the labeled layer data to the tracker output and
% provide a metric or score.

% For binary comparison

% 1. Do fir_dec in fast time
% Points for getting within N pixels of layer
%   [0 0 0 1 2 3 2 1 0 0 0 0 1 2 3 2 1 0 0 0 ... ]

% 2. Threshold values
% Support threshold
%   [0 0 0 1 1 1 1 1 0 0 0 0 1 1 1 1 1 0 0 0 ... ]

% 3. Apply secondary operation
% And any other functions that are vectorized easily such as square
%   [0 0 0 1 4 9 4 1 0 0 0 0 1 4 9 4 1 0 0 ... ]

% Score is a distance measure between these two
% Positive true
%   True positive: NN predicts layer and there really is a layer
%   True negative: NN predicts no layer and there really is no layer
%   False negative: NN predicts no layer and there really is a layer
%   False positive: NN predicts layer and there really is no layer
% positive/negative correlates to NN prediction
% true/false correlates to whether or not the NN predicted it correctly

% /cresis/snfs1/scratch/paden/ct_user_tmp/tracker_sim/masoud/layer/
% /cresis/snfs1/scratch/paden/ct_user_tmp/tracker_sim/masoud/layer_binary/
% /cresis/snfs1/scratch/paden/ct_user_tmp/tracker_sim/masoud/image/

% For vector comparison
% Find the closest layer on average (various )

