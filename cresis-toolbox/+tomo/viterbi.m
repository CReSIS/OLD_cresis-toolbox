% Perform ice-bottom tracking on 2D MUSIC imagery with the Viterbi algorithm
% Can be used on 2D echograms or 2D slices from 3D imagery
%
% Compile with mex -v -largeArrayDims viterbi.cpp
%    -v: verbose
%    -largeArrayDims: required for 64 bit
% 
% Usage: 
% [labels] = viterbi(input_img, surface_gt, bottom_gt, extra_gt, ice_mask, 
%   mean, var, egt_weight, smooth_weight, smooth_var, smooth_slope, 
%   bounds, viterbi_weight, repulsion, ice_bin_thr, mc, mc_weight)
%
% No input arguments may be missing 
%
% All ground truth (surface, bottom, and extra) must have values which are
% less than Nt-length(mu) or they will have no effect.
%
% input_img: double, 2D image 
% surface_gt: double, 2D surface, 1 by Nx vector, contains surface ground
%   truth
% bottom_gt: double, scalar, contains center ground truth bin
% extra_gt: double, 2 by N matrix where N is the number of extra ground
%   truth points, contains extra ground truth
% ice_mask: double, 2D surface mask, 1 by Nx vector, contains surface mask
%   A value of zero means no ice and the surface and bottom are forced to
%   be the same.
% mean: double, Nmean size vector, image magnitude template in Nt dimension
% var: double, Nmean size vector, image magnitude weights
%   corresponding to mean
% egt_weight: double, scalar, contains extra_gt weight 
% smooth_weight: double, scalar, contains smoothness weight
% smooth_var: double, scalar, contains smoothness variance
% smooth_slope: double, 1 by Nx-1 vector, contains expected slope 
% bounds: double, 1 by 2 vector, default is [0 Nx]. It should contain
%  the start and stop x-axis indexes. All indexes are zero-indexed. Setting
%  any element to a negative value will cause the default to be used for 
%  that particular element.
% viterbi_weight: double, 1 by Nx vector. Contains weighting factor for all
%  points. Default is ones(1 Nx). 
% repulsion: double, scalar. Contains scaling factor for surface repulsion.
% ice_bin_thr: double, scalar. Contains scanning range for proximity to
%  location of zero ice_mask (no ice). Allows for a more smooth transition
%  between icy and non-icy regions.
% mc: double, 1 by Nx vector. Additional bottom-layer evidence provided by
%  mass-conservation methods. Set values to -1 to ignore mass-conservation. 
% mc_weight: double, scalar. Contains scaling factor for mass-conservation
%  evidence. Set value to 0 to ignore mass-conservation.
%
% labels: 1 by Nx, contains tracked ice bottom layer
%
% surf: double, 2D surface, Ndoa by Nx matrix, contains ice bottom surface
%
% Authors:  Victor Berger, Mingze Xu, John Paden
