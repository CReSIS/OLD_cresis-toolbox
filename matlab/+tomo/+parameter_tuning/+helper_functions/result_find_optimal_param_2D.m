function [ result ] = result_find_optimal_param_2D( data_matrix, params, param_override, options, geotiff_fn, geotiff2_fn,...
  OPS_Surface, OPS_Bottom, OPS_data, OPS_crossover_data)
%RESULT_FIND_OPTIMAL_PARAM_2D 
%   Process the result from the test and find the optimal parameters using
%   fmincon

%   Input:
%     data_matrix: a matrix; result we got from the grid/random test
%     Viterbi_param: for getting the default parameter values
%     sources: the source (image) files we've setup according to the
%     user_setting 
%     reference: the reference (ground truth) files we've setup according to the
%     user_setting 
%   Output:
%     result: the optimal set of parameter and its error measures returned by 
%     the constraint minimization function in matlab (fmincon)

    iterations = 100;

    fprintf('Search complete. Now use fmincon to find the best set of parameters close to the region of the best set of parameters we have found with respect to the same data set.\n\n'); 
    
    [~,rank_indices] = sort(data_matrix(:,end-1), 'ascend');  % rank by mean error from low to high 
    data_matrix = data_matrix(rank_indices,:);                % make data matrix according to the ranking
    
    best_detect_param.sw = data_matrix(1,1);
    best_detect_param.sv = data_matrix(1,2);
    best_detect_param.repls = data_matrix(1,3);
    best_detect_param.ibt = data_matrix(1,4);
    
    best_stats_param.counter_f = 1;
    best_stats_param.total_diff = [];
    best_stats_param.hit_ratio_flag = true;
    best_stats_param.hit_ratio = [];
    best_stats_param.comb_vector = data_matrix(1,1:4);
        
    ObjFcn = @(x) Viterbi_handle_2D(params, param_override, options, geotiff_fn, geotiff2_fn,...
      best_detect_param, best_stats_param, OPS_Surface, OPS_Bottom, OPS_data, OPS_crossover_data, x(1), x(2), x(3));
    
    x0 = [best_detect_param.sw;...
      best_detect_param.sv;...
      best_detect_param.repls];                              % Initial point for the optimizer 
    A = [];
    B = [];
    Aeq = [];
    Beq = [];
    LB = [0; 0; 0];                                          % Lower bound for the input vector     
    UB = x0 + 500*ones(3,1);                                 % Upper bound for the input vector
    NONLCON = [];
    OPTIONS = optimset('MaxFunEval',iterations ,'DiffMinChange', 1); % step size is one
    
    x = fmincon(ObjFcn, x0, A, B, Aeq, Beq, LB, UB, NONLCON, OPTIONS);
    
    best_detect_param.sw = x(1);
    best_detect_param.sv = x(2);  
    best_detect_param.repls = x(3); 
    best_stats_param.comb_vector = [x(1) x(2) x(3) data_matrix(1,4)];
    best_stats_param = cluster_kernel_Viterbi_2D(params, param_override, options, geotiff_fn, geotiff2_fn, best_detect_param, best_stats_param,...
      OPS_Surface, OPS_Bottom, OPS_data, OPS_crossover_data);
    best_stats_param = compute_hit_ratios(best_stats_param);
    best_stats_param = compute_errors(best_stats_param);
    
    result = [best_stats_param.comb_vector , best_stats_param.hit_ratios, ...
          best_stats_param.rmse, best_stats_param.mean_difference, ...
          best_stats_param.median_difference];
    result_visualize(result, [], 'Viterbi', '2D');  
end


function [mean_error] = Viterbi_handle_2D(params, param_override, options, geotiff_fn, geotiff2_fn, ...
  best_detect_param, best_stats_obj, OPS_Surface, OPS_Bottom, OPS_data, OPS_crossover_data, sw, sv, repls)
  best_detect_param.sw = sw;
  best_detect_param.sv = sv;
  best_detect_param.repls = repls;  
  best_stats_obj.comb_vector = [sw sv repls best_detect_param.ibt];
  stats_struct = cluster_kernel_Viterbi_2D(params, param_override, options, geotiff_fn, geotiff2_fn,...
  best_detect_param, best_stats_obj, OPS_Surface, OPS_Bottom, OPS_data, OPS_crossover_data);
  mean_error = nanmean(stats_struct.total_diff(:));
  fprintf('Mean Error: %d \n', mean_error);
end
