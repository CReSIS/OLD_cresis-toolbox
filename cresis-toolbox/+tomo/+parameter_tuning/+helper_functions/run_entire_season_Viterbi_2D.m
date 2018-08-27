function [ result ] = run_entire_season_Viterbi_2D( params, param_override, OPS_Surface, OPS_Bottom, OPS_data, ...
  OPS_crossover_data, param_vec )
%RUN_ENTIRE_SEASON_VITERBI_2D Summary of this function goes here
%   Run the entire season on Antartica

%% setup  (refer to run_tracker.m)

options.name       = 'CSARP_post/mvdr';
options.debug      = false;
geotiff_fn         = 'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_icemask_grounded_and_shelves.tif';
geotiff2_fn        = 'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_rockmask.tif';
% geotiff_fn         = 'greenland/IceMask/GimpIceMask_90m_v1.1.tif';
geotiff_fn         = ct_filename_gis([],geotiff_fn);
geotiff2_fn         = ct_filename_gis([],geotiff2_fn);
% geotiff2_fn        = '';


%% experiment starts (refer to tracker.m)

physical_constants;

warning('off','all');

Viterbi_param.smooth_weight = param_vec(1);
Viterbi_param.smooth_var = param_vec(2);
Viterbi_param.repulsion = param_vec(3);
Viterbi_param.ice_b_thr = param_vec(4);

[detect_params_array, stats_array, ~] = ...
setup_parameters(Viterbi_param, [], 'grid', 'Viterbi', '2D', []);  

                            

%% do the job
result_stat_struct = cluster_kernel_viterbi_2D(params, param_override, options,...
  geotiff_fn, geotiff2_fn, detect_params_array(1), stats_array(1), OPS_Surface, OPS_Bottom, OPS_data, OPS_crossover_data );

result_stat_struct = compute_hit_ratios(result_stat_struct);
result_stat_struct = compute_errors(result_stat_struct);

result = [result_stat_struct.comb_vector , ...
  result_stat_struct.hit_ratios, ...
  result_stat_struct.rmse, ...
  result_stat_struct.mean_difference, ...
  result_stat_struct.median_difference];  

end

