function labels = tracker_viterbi(data,track)
% labels = tracker_viterbi(data,track)
%

gt = track.crossovers;
surf_bins = track.dem;
ice_mask.mask = track.ice_mask;

%% Distance-to-Ice-Margin model
if ~isfield(track.viterbi, 'DIM_matrix') || isempty(track.viterbi.DIM_matrix)
  prob.DIM_means = [6.908407 14.603709 22.077745 29.333980 36.375879 43.206908 49.830531 56.250214 62.469423 68.491622 74.320277 79.958853 85.410815 90.679629 95.768760 100.681673 105.421833 109.992706 114.397756 118.640450 122.724252 126.652627 130.429042 134.056960 137.539848 140.881171 144.084393 147.152981 150.090399 152.900113 155.585588 158.150289 160.597681 162.931230 165.154401 167.270659 169.283470 171.196299 173.012610 174.735870 176.369543 177.917095 179.381991 180.767696 182.077676 183.315395 184.484320 185.587915 186.629645 187.612977 188.541374 189.418303 190.247228 191.031615 191.774929 192.480636 193.152200 193.793087 194.406763 194.996691 195.566338 196.119170 196.658650 197.188245 197.711419 198.231639 198.752368 199.277073 199.809219 200.352271 200.909693 201.484953 202.081514 202.702842 203.352402 204.033660 204.750081 205.505130 206.302272 207.144972 208.036696 208.980910 209.981077 211.040664 212.163136 213.351958 214.610596 215.942514 217.351178 218.840053 220.412604 222.072297 223.822597 225.666969 227.608878 229.651790 231.799170 234.054483 236.421194 238.902770];
  prob.DIM_vars  = [19.032807 25.055615 28.640839 32.753438 36.144273 39.329838 42.688930 45.511082 49.036818 52.177650 55.028380 57.516627 60.519048 63.217206 65.211203 67.459337 69.609678 71.543557 73.182822 74.615772 75.628159 77.127086 78.155483 79.447090 80.011376 81.108576 81.618789 82.287856 82.979740 83.561585 84.281769 84.648076 85.290095 85.566969 86.052342 86.487424 86.675812 86.959733 87.181337 87.641261 87.674246 87.947628 87.895269 87.286380 87.202972 86.878606 87.151259 87.477659 88.049960 88.587946 88.515276 89.070799 88.756636 88.345201 87.754785 87.689382 87.240118 86.800999 86.164340 86.085916 85.803664 85.356194 85.831974 85.264038 85.222428 84.898093 84.652262 84.332790 84.249144 83.871931 83.552786 83.233334 82.842279 82.658637 82.008042 81.694151 81.421515 80.901673 80.885452 81.070003 80.524210 80.776716 80.320438 80.445820 80.085639 79.751146 79.557559 78.923447 78.522063 77.525973 77.426494 76.624448 76.855826 77.277564 76.777165 76.716292 75.970217 77.149291 76.900846 76.890210];
  gauss = @(x, mean, var)((1 / (sqrt(2 * pi * (var.^2)))) * (exp(-((x - mean).^2)/(2 * var.^ 2))));
  costm = ones(50, 100);
  for t = 1:50
    for i = 1:100
      costm(t,i) = 10 * exp(-0.01 .* t) - 10 * exp(-0.01 * 50);
    end
  end
  alpha = 0.025;
  beta = 0.05;
  
  DIM_costmatrix = ones(1000, 100);
  for DIM = 1 : 100;
    for T = 1 : 1000
      if alpha / gauss(T, prob.DIM_means(DIM), DIM * beta * prob.DIM_vars(DIM)) >= 200
        DIM_costmatrix(T, DIM) = 200;
      else
        DIM_costmatrix(T, DIM) = alpha/ gauss(T, prob.DIM_means(DIM), DIM * beta * prob.DIM_vars(DIM));
      end
    end
  end
  
  DIM_costmatrix = DIM_costmatrix(15:end,:);
  
  for k = 1:100
    DIM_costmatrix(1:50, k) = DIM_costmatrix(1:50, k) + k * costm(1:50, k);
  end
else
  global gRadar;
  DIM = load(fullfile(gRadar.path, track.viterbi.DIM_matrix));
  DIM_costmatrix = DIM.Layer_tracking_2D_parameters;
  DIM_costmatrix = DIM_costmatrix .* (200 ./ max(DIM_costmatrix(:)));
end

Nx = size(data, 2);
slope = diff(surf_bins);

try
  surf_weight = track.viterbi.surf_weight;
catch ME
  surf_weight = 1000;
end
try
  mult_weight = track.viterbi.mult_weight;
catch ME
  mult_weight = 100;
end
try
  mult_weight_decay = track.viterbi.mult_weight_decay;
catch ME
  mult_weight_decay = 0;
end
try
  mult_weight_local_decay = track.viterbi.mult_weight_local_decay;
catch ME
  mult_weight_local_decay = .8;
end
try
  manual_slope = track.viterbi.manual_slope;
catch ME
  manual_slope = 0;
end
try
  max_slope = track.viterbi.max_slope;
catch ME
  max_slope = -1;
end
try
  transition_weight = track.viterbi.transition_weight;
catch ME
  transition_weight = 1;
end
try
  image_mag_weight = track.viterbi.image_mag_weight;
catch ME
  image_mag_weight = 1;
end
try
  gt_weight = -track.viterbi.gt_weight;
catch ME
  gt_weight = -1;
end
try
  gt_cutoff = track.viterbi.gt_cutoff;
catch ME
  gt_cutoff = 5;
end

transition_weights = ones(1, Nx-1) * transition_weight;

if ~track.viterbi.use_surf_for_slope
  slope(:) = 0;
end

zero_bin = track.zero_bin;

%% Set variable echogram tracking parameters
bounds = [0 Nx];
ice_mask.mask_dist = round(bwdist(ice_mask.mask == 0));

ice_mask.mask_dist = round(bwdist(ice_mask.mask == 0));
ice_mask.mask_dist = round(ice_mask.mask_dist ./ 45);

DIM_costmatrix = 18 .* DIM_costmatrix;

%% Call viterbi.cpp

gt_layer = ones(1, size(surf_bins, 2)) * NaN;
gt_layer(gt(1, :)) = gt(2, :);
layers = [surf_bins; gt_layer];

surf_costs = ones(1, size(surf_bins, 2)) * surf_weight;
gt_costs = ones(1, size(surf_bins, 2)) * gt_weight;
layer_costs = [surf_costs; gt_costs];

surf_cutoffs = ones(1, size(surf_bins, 2)) * -1;
gt_cutoffs = ones(1, size(surf_bins, 2)) * gt_cutoff;
layer_cutoffs = [surf_cutoffs; gt_cutoffs];

labels = tomo.viterbi(double(data), double(layers), ...
  double(layer_costs), double(layer_cutoffs), double(ice_mask.mask), ...
  double(image_mag_weight), double(slope), double(max_slope), ...
  int64(bounds), double(ice_mask.mask_dist), double(DIM_costmatrix), ...
  double(transition_weights), double(mult_weight), ...
  double(mult_weight_decay), double(mult_weight_local_decay), int64(zero_bin));
