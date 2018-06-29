function [ result ] = run_entire_season_Viterbi_3D( param_vec )
%RUN_ENTIRE_SEASON Summary of this function goes here
%   This function will run in parallel

%% the entire data 
 
param.radar_name = 'rds';
param.season_name = '2014_Greenland_P3';
param.out_type = 'CSARP_post/music3D';
param.surfdata_source = 'CSARP_post/surfData';
%bounds_relative = [3 2 0 0];

% the entire 2014 Greenland P3 data
segments_and_frame = containers.Map;
segments_and_frame('20140325_05') = 1:2; % doesn't work
segments_and_frame('20140325_06') = 1;
segments_and_frame('20140325_07') = 1:5;
segments_and_frame('20140401_03') = 1:48;
segments_and_frame('20140506_01') = 1:46;

sources = {};
references = {};

for segment_name = keys(segments_and_frame)
  param.day_seg = segment_name{:};                                            % get the string
  
  for frame = segments_and_frame(segment_name{:})
    source_fn = fullfile(ct_filename_out(param,param.out_type,''), ...
      sprintf('Data_%s_%03d.mat',param.day_seg,frame));                       % get the file path as a string from the image data
    ref_fn = fullfile(ct_filename_out(param,param.surfdata_source, ...
      'CSARP_surfData'),sprintf('Data_%s_%03d.mat',param.day_seg,frame));     % get the reference path as a string
    sources = [sources {source_fn}];
    references = [references {ref_fn}];
  end
end

num_slices = 330712; % pre-calculated before

%% parameters setup
% default param
Viterbi_param.threshold = 13.5;
Viterbi_param.slope = zeros(1,63);
Viterbi_param.num_columns = [];
Viterbi_param.previous = false;

% tuned param
Viterbi_param.smooth_weight = param_vec(1);
Viterbi_param.smooth_var = param_vec(2);
Viterbi_param.repulsion = param_vec(3);
Viterbi_param.ice_b_thr = param_vec(4);
Viterbi_param.egt_weight = param_vec(5);

%% create parameter structures
% 'grid' is fine here since we're running only one combination
[detect_params_array, stats_array, ~]...
  = setup_parameters(Viterbi_param, num_slices, 'grid', 'Viterbi', '3D', []);  % setup the data structures for cluster processing

result_stat_struct = cluster_kernel_Viterbi_3D(sources, references, detect_params_array(1), stats_array(1));
result_stat_struct = compute_hit_ratios(result_stat_struct);
result_stat_struct = compute_errors(result_stat_struct);

result = [result_stat_struct.comb_vector , ...
  result_stat_struct.hit_ratios, ...
  result_stat_struct.rmse, ...
  result_stat_struct.mean_difference, ...
  result_stat_struct.median_difference];  

end

