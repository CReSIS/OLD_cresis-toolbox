%[] = detect(input_image, surface_gt, bottom_gt, extra_gt, ice_mask, 
%   mean, var, mid, egt_weight, smooth_weight, smooth_var, smooth_slope)
%
%HMM Viterbi detection function for ice-bed layer in MUSIC slices
%
%Compile with mex -v -largeArrayDims detect.cpp
%    -v: verbose
%    -largeArrayDims: required for 64 bit
% 
%Inputs
% input_img: radar echogram to be analyzed (double matrix Nt by N)
% surface_gt: surface bins matrix (double matrix 1 by N)
% bottom_gt: bottom bin (scalar for bottom at column 32)
% extra_gt: manual bottom points (double matrix 2 by Ngt)
%   extra_gt(1,:): indicates the column in input_img of the ground truth point
%   extra_gt(2,:): indicates the row in input_img  of the ground truth point
% ice_mask: ice(true) or rock(false) mask (double matrix 1 by N)
% mean: mean of image peak template (double matrix 1 by Ntemplate)
% var: variance of image peak template (double matrix 1 by Ntemplate)
% mid: -1
% egt_weight: weight attributed to manual ground truth points (e.g. 10).
% smooth_weight: -1
% smooth_var: -1
% smooth_slope: to assume topography is flat (no slope): zeros(1, N - 1)
%
%Outputs
% labels: location of bottom layer for each column (double 1 by N)
%
%Examples: 
%  /+imb/slicetool_detect.m 
%  /+imb/@picktool_detect/left_click_and_drag.m
%
% Authors:  Mingze Xu
