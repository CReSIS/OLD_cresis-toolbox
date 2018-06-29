%% data setup

% Here, we select the source data and the reference (i.e. ground truth/manually
% picked data) for the grid or random search we'd like to test on

param.radar_name = 'rds';
param.season_name = '2014_Greenland_P3';
param.out_type = 'CSARP_post/music3D';
param.surfdata_source = 'CSARP_post/surfData';
bounds_relative = [3 2 0 0];

% Put the segment and frame in a "hash map"
% The key is the segment
% The values with respect to the key are the frames under that segment

segments_and_frame = containers.Map;
segments_and_frame('20140325_05') = [1];
% segments_and_frame('20140325_06') = [];
% segments_and_frame('20140325_07') = [];
% segments_and_frame('20140401_03') = [];
% segments_and_frame('20140506_01') = [46];
% segments_and_frame('20140325_06') = [1];
% segments_and_frame('20140325_07') = [5];
% segments_and_frame('20140401_03') = [33];
% segments_and_frame('20140506_01') = [45 46];
%% number of trials
% Here we set up how many random vectors we'd like to generate
% This is equivalent to how many tests we want to do in a random search

number_of_trials = 1;

%% range setup
% Minimum and maximum value we want to setup for each component in a
% random vector; a random vector here is simply a combination of input parameters
% for Viterbi Algorithm

Viterbi_param.smooth_weight_min = 0;
Viterbi_param.smooth_weight_max = 200;

Viterbi_param.smooth_var_min = 0;
Viterbi_param.smooth_var_max = 200;
  
Viterbi_param.repulsion_min = 0;
Viterbi_param.repulsion_max = 3000;

Viterbi_param.egt_weight_min = 2;
Viterbi_param.egt_weight_max = 30;

% it will be random integer
Viterbi_param.ice_b_thr_min = 0;
Viterbi_param.ice_b_thr_max = 5; 

