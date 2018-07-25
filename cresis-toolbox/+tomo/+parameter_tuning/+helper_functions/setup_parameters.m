function [ detect_params_array, stats_array, num_combinations ] = setup_parameters( algo_param, num_images, method, algorithm, data_type, number_of_trials)
%SETUP_3D_PARAMETERS Summary of this function goes here
% Setup the combinations of parameters to feed in to the tests later on
if strcmp(data_type, '3D')
  if strcmp(method, 'grid')
    if strcmp(algorithm, 'TRWS')
      %% Create parameter combinations fpr grid search of parameter setup on 3D data
      combinations = combvec(algo_param.correlation, algo_param.smooth_1, algo_param.smooth_2, algo_param.smooth_3)';
      [detect_params_array, stats_array, num_combinations] = set_data_structures(combinations, algo_param, num_images, 'TRWS', '3D');
    elseif  strcmp(algorithm, 'Viterbi')
      %% Viterbi grid search combinations of parameter setup on 3D data
      combinations = combvec(algo_param.smooth_weight, algo_param.smooth_var, algo_param.repulsion, algo_param.ice_b_thr, algo_param.egt_weight)';
      [detect_params_array, stats_array, num_combinations] = set_data_structures(combinations, algo_param, num_images, 'Viterbi', '3D');
    end
  elseif strcmp(method, 'random')
    if strcmp(algorithm, 'Viterbi')
      %% Viterbi random search combinations of parameter setup on 3D data
      random_vectors = make_random_vectors(algo_param, number_of_trials, 'Viterbi', '3D');
      [detect_params_array, stats_array, num_combinations] = set_data_structures(random_vectors, algo_param, num_images, 'Viterbi', '3D');
    elseif strcmp(algorithm, 'TRWS')
      %% TRWS random search combinations of parameter setup on 3D data
      random_vectors = make_random_vectors(algo_param, number_of_trials, 'TRWS', '3D');
      [detect_params_array, stats_array, num_combinations] = set_data_structures(random_vectors, algo_param, num_images, 'TRWS', '3D');
    end
  end
elseif strcmp(data_type, '2D')
  if strcmp(method, 'grid') && strcmp(algorithm, 'Viterbi')
    %% Viterbi -- grid search combinations of parameter setup on 2D data
    combinations = combvec(algo_param.smooth_weight, algo_param.smooth_var, algo_param.repulsion, algo_param.ice_b_thr)';
    [detect_params_array, stats_array, num_combinations] = set_data_structures(combinations, algo_param, [], 'Viterbi', '2D');
  elseif strcmp(method, 'random') && strcmp(algorithm, 'Viterbi')
    %% Viterbi  -- random search combinations of parameter setup on 2D data
    
    random_vectors = make_random_vectors(algo_param, number_of_trials, 'Viterbi', '2D');
    
    [detect_params_array, stats_array, num_combinations] = set_data_structures(random_vectors, algo_param, [], 'Viterbi', '2D');
  end
end
end

function [detect_array, statistic_array, num_comb] = set_data_structures(combinations, algorithm_param, number_of_images, algorithm_name, data_type)

num_comb = size(combinations, 1);
detect_array = struct.empty(0,num_comb);
statistic_array = struct.empty(0,num_comb);

for idx_i = 1:num_comb
  if strcmp(data_type, '2D')
    detect_array(idx_i).sw = combinations(idx_i,1);
    detect_array(idx_i).sv = combinations(idx_i,2);
    detect_array(idx_i).repls = combinations(idx_i,3);
    detect_array(idx_i).ibt = combinations(idx_i,4);
    
    statistic_array(idx_i).counter_f = 1;
    statistic_array(idx_i).total_diff = [];
    statistic_array(idx_i).hit_ratio_flag = true;
    statistic_array(idx_i).hit_ratio = [];
    statistic_array(idx_i).comb_vector = combinations(idx_i,:);
    
  elseif strcmp(data_type, '3D')
    if strcmp(algorithm_name, 'TRWS')
      detect_array(idx_i).correlation = combinations(idx_i,1);
      detect_array(idx_i).smooth_1 = combinations(idx_i,2);
      detect_array(idx_i).smooth_2 = combinations(idx_i,3);
      detect_array(idx_i).smooth_3 = combinations(idx_i,4);
      detect_array(idx_i).slope = algorithm_param.slope;
      detect_array(idx_i).cols = algorithm_param.num_columns;
      detect_array(idx_i).threshold = algorithm_param.threshold;
    elseif strcmp(algorithm_name, 'Viterbi')
      detect_array(idx_i).previous = algorithm_param.previous;
      detect_array(idx_i).threshold = algorithm_param.threshold;
      detect_array(idx_i).slope = algorithm_param.slope;
      detect_array(idx_i).cols = algorithm_param.num_columns;
      detect_array(idx_i).viterbi_weight = [];
      detect_array(idx_i).smooth_weight = combinations(idx_i,1);
      detect_array(idx_i).smooth_var = combinations(idx_i,2);
      detect_array(idx_i).repulsion = combinations(idx_i,3);
      detect_array(idx_i).ice_b_threshold = combinations(idx_i,4);
      detect_array(idx_i).egt_weight = combinations(idx_i,5);
      
      info_v = sprintf('smooth weight: %d; smooth var: %d; repulsion: %d, ice distance threshold: %d egt_weight: %d', ...
        detect_array(idx_i).smooth_weight, ...
        detect_array(idx_i).smooth_var, ...
        detect_array(idx_i).repulsion, ...
        detect_array(idx_i).ice_b_threshold, ...
        detect_array(idx_i).egt_weight );
      statistic_array(idx_i).title = info_v;
    end
    
    statistic_array(idx_i).rmse_f = double.empty(0, number_of_images);
    statistic_array(idx_i).diff_f = double.empty(0, number_of_images);
    statistic_array(idx_i).med_f = double.empty(0, number_of_images);
    statistic_array(idx_i).max_df = double.empty(0, number_of_images);
    statistic_array(idx_i).min_df = double.empty(0, number_of_images);
    statistic_array(idx_i).counter_f = 1;
    statistic_array(idx_i).total_diff = [];
    statistic_array(idx_i).comb_vector = combinations(idx_i,:);
    statistic_array(idx_i).hit_ratio = [];
    statistic_array(idx_i).hit_ratio_flag = 1;
  end
end
end


function [random_vecs] = make_random_vectors(algo_param, num_trials, algo_name, data_type)
rng('shuffle');               % seeds the random number generator based on the current time
random_vecs = [];
for index = 1:num_trials
  if strcmp(algo_name, 'TRWS')
    correlation = round((algo_param.correlation_max-algo_param.correlation_min).*rand(1,1) + algo_param.correlation_min);
    smooth_1 = (algo_param.smooth_1_max-algo_param.smooth_1_min).*rand(1,1) + algo_param.smooth_1_min;
    smooth_2 = (algo_param.smooth_2_max-algo_param.smooth_2_min).*rand(1,1) + algo_param.smooth_2_min;
    smooth_3 = (algo_param.smooth_3_max-algo_param.smooth_3_min).*rand(1,1) + algo_param.smooth_3_min;    
    random_vecs = vertcat(random_vecs, [correlation smooth_1 smooth_2 smooth_3]);
    
  elseif strcmp(algo_name, 'Viterbi')
    repulsion = (algo_param.repulsion_max-algo_param.repulsion_min).*rand(1,1) + algo_param.repulsion_min;
    smooth_weight = (algo_param.smooth_weight_max-algo_param.smooth_weight_min).*rand(1,1) + algo_param.smooth_weight_min;
    smooth_var = (algo_param.smooth_var_max-algo_param.smooth_var_min).*rand(1,1) + algo_param.smooth_var_min;
    ice_b_thr = randi([algo_param.ice_b_thr_min algo_param.ice_b_thr_max],1,1);
    
    if strcmp(data_type, '3D')
      egt_weight = (algo_param.egt_weight_max-algo_param.egt_weight_min).*rand(1,1) + algo_param.egt_weight_min;
      
      random_vecs = vertcat(random_vecs, [smooth_weight smooth_var repulsion ice_b_thr egt_weight]);
    elseif  strcmp(data_type, '2D')
      random_vecs = vertcat(random_vecs, [smooth_weight smooth_var repulsion ice_b_thr]);
    end
  end
end
end