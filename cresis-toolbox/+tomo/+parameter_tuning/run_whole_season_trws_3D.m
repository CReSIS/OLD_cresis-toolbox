
% 3D data with TRWS
% Retrieve data from the test we're interested
trws_grid_results_path = '../../../../tuning_results/grid_search_trws_3D/';
trws_random_results_path = '../../../../tuning_results/random_search_trws_3D/';
trws_whole_results_path = '../../../../tuning_results/Entire_season_TRWS_3D.mat';

listings_grid = dir(trws_grid_results_path);
listings_random = dir(trws_random_results_path);

best_parameters = containers.Map;
best_parameters('grid') = []; 
best_parameters('random') = [];


% load the current saved file
parameter_tested_before = [];
parameter_tested_results_before = [];
try
  saved_file = load(trws_whole_results_path);  
  parameter_tested_results_before = table2array(saved_file.result_Viterbi_3D);
  parameter_tested_before = table2array(saved_file.result_Viterbi_3D(:,1:5));
  % this are the parameters that had run the whole season before;
  % i.e. we know their results -- we don't want to run them again unless
  % the whole season is set to be something else
catch
  fprintf('No saved file on the whole season result for TRWS \n');
end


for idx = 3:length(listings_grid)
  file_struct = listings_grid(idx);
  file = load(strcat(trws_grid_results_path, file_struct.name));
  best_result = table2array(file.best_parameter_table);
  best_result = best_result(1:4);
    
  if exist(trws_whole_results_path,'file') && ismember(best_result, parameter_tested_before, 'rows')
    continue;
  else
    best_parameters('grid') = [best_parameters('grid'); best_result];
  end    
end

for idx = 3:length(listings_random)
  file_struct = listings_random(idx);
  file = load(strcat(trws_random_results_path, file_struct.name));
  best_result = table2array(file.best_parameter_table);
  best_result = best_result(1:4);
  if exist(trws_whole_results_path,'file') && ismember(best_result, parameter_tested_before, 'rows') 
    continue;
  else
    best_parameters('random') = [best_parameters('random'); best_result];
  end
end

all_parameters = [best_parameters('grid'); best_parameters('random')];
%% run the entire season 

% get the total number of best combination
num_combs = size(all_parameters,1);
results = double.empty(0,13); % 13 for TRWS

parfor idx = 1:num_combs
  results(idx,:) = run_entire_season_TRWS_3D(all_parameters(idx,:));
end

results = [parameter_tested_results_before; results]; % concatenate the results we've had
result_TRWS_3D = rank_result(results, 'TRWS', '3D');
save(trws_whole_results_path, 'result_TRWS_3D');