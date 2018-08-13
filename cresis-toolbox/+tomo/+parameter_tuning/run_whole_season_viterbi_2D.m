
% 3D data with Viterbi
% Retrieve data from the test we're interested

v_grid_results_path = '../../../../tuning_results/grid_search_viterbi_2D/';
v_random_results_path = '../../../../tuning_results/random_search_viterbi_2D/';
v_whole_results_path = '../../../../tuning_results/Entire_season_Viterbi_2D.mat';

listings_grid = dir(v_grid_results_path);
listings_random = dir(v_random_results_path);

best_parameters = containers.Map;
best_parameters('grid') = [];
best_parameters('random') = [];


% load the current saved file
parameter_tested_before = [];
parameter_tested_results_before = [];

try
  saved_file = load(v_whole_results_path);
  parameter_tested_results_before = table2array(saved_file.result_Viterbi_2D);
  parameter_tested_before = table2array(saved_file.result_Viterbi_2D(:,1:4));
  % this are the parameters that had run the whole season before;
  % i.e. we know their results -- we don't want to run them again unless
  % the whole season is set to be something else
catch
  fprintf('No saved file on the whole season result');
end


for idx = 3:length(listings_grid)
  file_struct = listings_grid(idx);
  file = load(strcat(v_grid_results_path, file_struct.name));
  best_result = table2array(file.best_parameter_table);
  best_result = best_result(1:4);
  if exist('Entire_season_Viterbi_2D.mat','file') && ismember(best_result, parameter_tested_before, 'rows')
    continue;
  else
    best_parameters('grid') = [best_parameters('grid'); best_result];
  end
end

for idx = 3:length(listings_random)
  file_struct = listings_random(idx);
  file = load(strcat(v_random_results_path, file_struct.name));
  best_result = table2array(file.best_parameter_table);
  best_result = best_result(1:4);
  if exist(v_whole_results_path,'file') && ismember(best_result, parameter_tested_before, 'rows')
    continue;
  else
    best_parameters('random') = [best_parameters('random'); best_result];
  end
end

all_parameters = [best_parameters('grid'); best_parameters('random')];
%% run the entire season

% get the total number of best combination
num_combs = size(all_parameters,1);
results = double.empty(0,13); % 13 for Viterbi 2D

% some setup before the run
clear('param_override');
global gRadar;
params = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_DC8_all.xls'),'','post');
% params = read_param_xls(ct_filename_param('rds_param_2017_Greenland_P3.xls'),'','post');

% Input checking
if ~exist('params','var')
  error('Use run_tracker: A struct array of parameters must be passed in\n');
end

if exist('param_override','var')
  param_override = merge_structs(gRadar, param_override);
else
  param_override = gRadar;
end


for i = 1:41
  if params(i).cmd.generic ==1
    fprintf('[%d] ', i);
  end
end

[OPS_Surface, OPS_Bottom, OPS_data, ...
  OPS_crossover_data, ~] = setup_OPS(params, param_override);


parfor idx = 1:num_combs
  results(idx,:) = run_entire_season_Viterbi_2D(params, param_override, ...
    OPS_Surface, OPS_Bottom, OPS_data, OPS_crossover_data, all_parameters(idx,:));
end

 
results = [parameter_tested_results_before; results]; % concatenate the results we've had
result_Viterbi_2D = rank_result(results, 'Viterbi', '2D');  % rank everything via the mean error
save(v_whole_results_path, 'result_Viterbi_2D');