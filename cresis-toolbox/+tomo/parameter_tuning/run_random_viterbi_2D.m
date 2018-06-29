%% This is a script to run parameter tunning
% WHAT IT DOES: RANDOM SEARCH on 2D data with Viterbi Algorithm
% WARNING: DO NOT RUN THIS ON STRANGERSON; it'll run out of memory

tic

%% user settings
user_settings_viterbi_random_2D;                                  % (FILE TO EDIT FOR EACH TEST) set up the image sources, references, and the parameter combinations              

%% Pre-load OPS files
[OPS_Surface, OPS_Bottom, OPS_data, ...
  OPS_crossover_data, total_bins] = setup_OPS(params, param_override);
                                                                  % preload OPS data since OPS server does not support parallel file access  

%%
cpu_time = 500 + .0032*total_bins;                                % estimated cpu_time for cluster processing
memory_requirement = (50)*100000;                                 % estimated memory usage for 2D images

%% create parameter structures
[detect_params_array, stats_array, ~] = setup_parameters(Viterbi_param, [], 'random', 'Viterbi', '2D', num_trials);
                                                                  % setup the data structures for cluster processing 

%% Cluster
                                                                                                             
result = cluster_Viterbi_2D( cpu_time, memory_requirement, detect_params_array, ...
  stats_array, num_trials, params, param_override, options, geotiff_fn, geotiff2_fn,...
  OPS_Surface, OPS_Bottom, OPS_data, OPS_crossover_data);         % run the test on different combination of parameters on the cluster   

%% data reorganization
data_matrix = result_make_data_matrix(result);                    % extract data from the cluster result and make a data matrix

%% data visualization
result_visualize(data_matrix, 'random', ...
  'Viterbi', '2D');                                               % show a table of the ranked result of rmse from low to high

%% find optimal set of parameters using fmincon

best_parameter_set = result_find_optimal_param_2D(data_matrix, ...
  params, param_override, options, geotiff_fn, geotiff2_fn, ...
  OPS_Surface, OPS_Bottom, OPS_data, OPS_crossover_data);         % run fmincon with the best point we've found (best set of parameters) on the response surface as the initial point 
                                                                  % this optimization use the same dataset
%% save the result
save_result(data_matrix, '2D','Viterbi', detect_params_array,...
  'random search', params, best_parameter_set);                   % save the test result in the RESULTS folder

fprintf('Parameter tuning complete. Please check the saved file in the tuning_results folder under /script \n');
toc
