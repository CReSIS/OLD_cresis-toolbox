%% This is a script to run parameter tunning 
% WHAT IT DOES: RANDOM SEARCH on 3D data with Viterbi Algorithm
tic
%% User Settings
user_settings_viterbi_random_3D;                            % (FILE TO EDIT FOR EACH TEST) set up the image sources, references, and the parameter combinations
[sources, references, num_slices] = ...
  setup_3Ddata(param, segments_and_frame);                  % setup the file paths for 3D images and references for cluster processing 

%% default parameters
% default parameters are parameters that we choose NOT to 
% include in the random search; they are fixed as constants throughout the
% tests
Viterbi_param.threshold = 13.5;
Viterbi_param.slope = zeros(1,63);
Viterbi_param.num_columns = [];
Viterbi_param.previous = false;

%% time for computation
% a heuristic to estimate the computational resource requred for each
% individual task
cpu_time = .0164*num_slices + 300;                          % estimated cpu_time (via one variable regression)
memory_requirement = (.3*num_slices)*1000000;               % estimated memory usage

%% create parameter structures

[detect_params_array, stats_array, num_combinations] = ...
  setup_parameters(Viterbi_param, num_slices, 'random', 'Viterbi', '3D', number_of_trials);
                                                            % setup the data structures for cluster processing      
                                                        
%% Initialize Cluster and start tasks
result = cluster_viterbi_3D(cpu_time, memory_requirement, detect_params_array, ...
  stats_array, num_combinations, sources, references);      % run the test on different combination of parameters on the cluster

%% data reorganization 
data_matrix = result_make_data_matrix(result);              % extract data from the cluster result and make a data matrix 

%% data visualization
result_visualize(data_matrix, 'random',...
  'Viterbi', '3D');                                         % show a table of the ranked result of rmse from low to high
                                                            % and show a plot for each parameter vs. the error      
%% find optimal set of parameters using fmincon              

best_parameter_set = result_find_optimal_param_3D(data_matrix, Viterbi_param,...
  sources, references, num_slices, 'Viterbi', '3D');        % run fmincon with the best point we've found (best set of parameters) on the response surface as the initial point 
                                                            % this optimization use the same dataset
                                                        

%% save the result
save_result(data_matrix, '3D', 'Viterbi', Viterbi_param,...
  'random search', segments_and_frame, best_parameter_set); % save the test result in the RESULTS folder

fprintf('Parameter tuning complete. Please check the saved file in the tuning_results folder under /script \n');
toc

                                                        
                                                        
