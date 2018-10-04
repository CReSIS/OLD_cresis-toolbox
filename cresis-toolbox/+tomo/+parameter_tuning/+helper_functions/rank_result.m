function [table_ranked_data] = rank_result( data_matrix, algo_name, data_type)
  if strcmp(data_type, '3D')
    if strcmp(algo_name, 'Viterbi')
      info = {'smooth_weight', 'smooth_variance', 'repulsion',...
        'ice_dist_thr', 'egt_weight', 'percentage_correct', 'bin_5_pass_rate', ...
        'bin_10_pass_rate', 'bin_15_pass_rate', ...
        'bin_20_pass_rate', 'bin_25_pass_rate', ...
        'rmse', 'mean_error','median_error'};

      till = min(size(data_matrix,1), 10);
      table_ranked_data = table(data_matrix(1:till,1),data_matrix(1:till,2), data_matrix(1:till,3), ...
        data_matrix(1:till,4), data_matrix(1:till,5), data_matrix(1:till,6), ...
        data_matrix(1:till,7), data_matrix(1:till,8), ...
        data_matrix(1:till,9),  data_matrix(1:till,10) ,...
        data_matrix(1:till,11),  data_matrix(1:till,12),  data_matrix(1:till,13), ...
        data_matrix(1:till,14), ...
        'VariableNames',info);
    end
    if strcmp(algo_name, 'TRWS')
       info = {'correlation', 'smooth_1', 'smooth_2',...
        'smooth_3', 'percentage_correct', 'bin_5_pass_rate', ...
        'bin_10_pass_rate', 'bin_15_pass_rate', ...
        'bin_20_pass_rate', 'bin_25_pass_rate', ...
        'rmse', 'mean_error','median_error'};
      
        till = min(size(data_matrix,1), 10);
        table_ranked_data = table(data_matrix(1:till,1),data_matrix(1:till,2), data_matrix(1:till,3), ...
                data_matrix(1:till,4), data_matrix(1:till,5), data_matrix(1:till,6), ...
                data_matrix(1:till,7), data_matrix(1:till,8), ...
                data_matrix(1:till,9),  data_matrix(1:till,10) ,...
                data_matrix(1:till,11),  data_matrix(1:till,12),  data_matrix(1:till,13), ...
                'VariableNames',info);
    end

  elseif strcmp(data_type, '2D')
    if strcmp(algo_name, 'Viterbi')
      info = {'smooth_weight', 'smooth_variance', 'repulsion',...
        'ice_dist_thr', 'percentage_correct', 'bin_5_pass_rate', ...
        'bin_10_pass_rate', 'bin_15_pass_rate', ...
        'bin_20_pass_rate', 'bin_25_pass_rate', ...
        'rmse', 'mean_error','median_error'};
      till = min(size(data_matrix,1), 10);
      table_ranked_data = table(data_matrix(1:till,1),data_matrix(1:till,2), data_matrix(1:till,3), ...
        data_matrix(1:till,4), data_matrix(1:till,5), data_matrix(1:till,6), ...
        data_matrix(1:till,7), data_matrix(1:till,8), ...
        data_matrix(1:till,9),  data_matrix(1:till,10) ,...
        data_matrix(1:till,11),  data_matrix(1:till,12),  data_matrix(1:till,13), ...
        'VariableNames',info);
    end
  end
  table_ranked_data
end