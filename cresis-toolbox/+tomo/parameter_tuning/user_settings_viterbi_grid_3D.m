
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
% segments_and_frame('20140325_05') = [2]; % doesn't work
% segments_and_frame('20140325_06') = [1];
% segments_and_frame('20140325_07') = [5];
% segments_and_frame('20140401_03') = [33];
segments_and_frame('20140506_01') = [46];


%% grid setup
Viterbi_param.smooth_weight = [120 125 ];
Viterbi_param.smooth_var = [30 35];
Viterbi_param.egt_weight = [0];
Viterbi_param.repulsion = [3000];
Viterbi_param.ice_b_thr = [0];
