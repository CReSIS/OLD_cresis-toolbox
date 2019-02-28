function [ result ] = result_find_optimal_param_3D(data_matrix, algo_param, sources, references, num_slices, algorithm_name, data_type)
%RESULT_FIND_OPT_PARAM 
%   Process the result from the test and find the optimal parameters using
%   fmincon

%   Input:
%     data_matrix: a matrix; result we got from the grid/random test
%     algo_param: for getting the default parameter values
%     sources: the source (image) files we've setup according to the
%     user_setting 
%     reference: the reference (ground truth) files we've setup according to the
%     user_setting 
%   Output:
%     result: the optimal set of parameter and its error measures returned by 
%     the constraint minimization function in matlab (fmincon)


    iterations = 100;

    fprintf('Search complete. Now use fmincon to find the best set of parameters close to the region of the best set of parameters we have found with respect to the same data set.\n\n'); 

    [~,rank_indices] = sort(data_matrix(:,end-1), 'ascend'); % rank by mean error from low to high 
    data_matrix = data_matrix(rank_indices,:);             % make data matrix according to the ranking
    
      
    if strcmp(algorithm_name, 'Viterbi')
      if strcmp(data_type, '3D')
        stats_struct = find_Viterbi_3D_best_param(algo_param, data_matrix, sources, references, num_slices);
        result = [stats_struct.comb_vector , stats_struct.hit_ratios, ...
          stats_struct.rmse, stats_struct.mean_difference, ...
          stats_struct.median_difference];
        result_visualize(result, [], 'Viterbi', '3D');
      elseif strcmp(data_type, '2D')
      end
    elseif strcmp(algorithm_name, 'TRWS')
      stats_struct = find_TRWS_3D_best_param(algo_param, data_matrix, sources, references, num_slices);
      result = [stats_struct.comb_vector , stats_struct.hit_ratios, ...
        stats_struct.rmse, stats_struct.mean_difference, ...
        stats_struct.median_difference];
      result_visualize(result, [], 'TRWS', '3D');
    end
           
    % Setup the parameters according to the best result from our test
    
end

function [best_stats_param] = find_TRWS_3D_best_param(algo_param, data_matrix, sources, references, num_slices)
    best_detect_param.threshold = algo_param.threshold;
    best_detect_param.slope = algo_param.slope;
    best_detect_param.cols = algo_param.num_columns;    
    best_detect_param.correlation = data_matrix(1,1); % this neesd to be an integer
    best_detect_param.smooth_1 = data_matrix(1,2);
    best_detect_param.smooth_2 = data_matrix(1,3);
    best_detect_param.smooth_3 = data_matrix(1,4);
       
    best_stats_param.rmse_f = double.empty(0, num_slices);
    best_stats_param.diff_f = double.empty(0, num_slices);
    best_stats_param.med_f = double.empty(0, num_slices);
    best_stats_param.max_df = double.empty(0, num_slices);
    best_stats_param.min_df = double.empty(0, num_slices);
    best_stats_param.counter_f = 1;
    best_stats_param.total_diff = [];
    best_stats_param.comb_vector = data_matrix(1,1:4);
    
      
    % Setup for fmincon
    ObjFcn = @(x) TRWS_handle_3D(sources, references, best_detect_param, ...
    best_stats_param, x(1), x(2), x(3));

    x0   = [best_detect_param.smooth_1;...
      best_detect_param.smooth_2;...
      best_detect_param.smooth_3];          % Initial point for the optimizer 
    A = [];
    B = [];
    Aeq = []; 
    Beq = [];
    LB = x0 - 30*ones(3,1);                 % Lower bound for the input vector     
    UB = x0 + 30*ones(3,1);                 % Upper bound for the input vector
    NONLCON = [];
    OPTIONS = optimset('MaxFunEval',iterations ,'DiffMinChange', 1); % step size is one
    
    x = fmincon(ObjFcn, x0, A, B, Aeq, Beq, LB, UB, NONLCON, OPTIONS);       
    
    best_detect_param.smooth_1 = x(1);
    best_detect_param.smooth_2 = x(2);  
    best_detect_param.smooth_3 = x(3); 
    best_stats_param.comb_vector = [data_matrix(1,1) x(1) x(2) x(3)];
    best_stats_param = cluster_kernel_TRWS_3D(sources, references, best_detect_param, best_stats_param);
    best_stats_param = compute_hit_ratios(best_stats_param);
    best_stats_param = compute_errors(best_stats_param);
end

function [mean_error] = TRWS_handle_3D(sources, refs, best_detect_param, best_stats_obj, s1, s2, s3)  
  best_detect_param.smooth_1 = s1;
  best_detect_param.smooth_2 = s2;
  best_detect_param.smooth_3 = s3;
  best_stats_obj.comb_vector = [best_detect_param.correlation s1 s2 s3];
  stats_struct = cluster_kernel_TRWS_3D(sources, refs, best_detect_param, best_stats_obj);
  mean_error = nanmean(stats_struct.total_diff(:));
  fprintf('Mean Error: %d \n', mean_error);
end

function [best_stats_param] = find_Viterbi_3D_best_param(algo_param, data_matrix, sources, references, num_slices)
    best_detect_param.previous = algo_param.previous;  
    best_detect_param.threshold = algo_param.threshold;
    best_detect_param.slope = algo_param.slope;
    best_detect_param.cols = algo_param.num_columns;
    best_detect_param.viterbi_weight = [];
    best_detect_param.smooth_weight = data_matrix(1,1);
    best_detect_param.smooth_var = data_matrix(1,2);
    best_detect_param.repulsion = data_matrix(1,3); 
    best_detect_param.ice_b_threshold = data_matrix(1,4);
    best_detect_param.egt_weight = data_matrix(1,5);  

    best_stats_param.previous = algo_param.previous; 
    best_stats_param.rmse_f = double.empty(0, num_slices);
    best_stats_param.diff_f = double.empty(0, num_slices);
    best_stats_param.med_f = double.empty(0, num_slices);
    best_stats_param.max_df = double.empty(0, num_slices);
    best_stats_param.min_df = double.empty(0, num_slices);
    best_stats_param.counter_f = 1;
    best_stats_param.total_diff = [];

    info = sprintf('smooth weight: %d; smooth var: %d; repulsion: %d, ice distance threshold: %d egt_weight: %d', ...
      best_detect_param.smooth_weight, ...
      best_detect_param.smooth_var, ...
      best_detect_param.repulsion, ...
      best_detect_param.ice_b_threshold, ...
      best_detect_param.egt_weight );

    best_stats_param.title = info;  
    best_stats_param.comb_vector = data_matrix(1,1:5);

    % Run the constraint minimization
    % We exclude the parameter 'Ice distance threshold' here because it is an integer
    
    % Setup for fmincon
    ObjFcn = @(x) Viterbi_handle_3D(sources, references, best_detect_param, ...
    best_stats_param, x(1), x(2), x(3), x(4));

    x0   = [best_detect_param.smooth_weight;...
      best_detect_param.smooth_var;...
      best_detect_param.egt_weight;...
      best_detect_param.repulsion];         % Initial point for the optimizer 
    A = [];
    B = [];
    Aeq = []; 
    Beq = [];
    NONLCON = []; 
    LB = x0 - 50*ones(4,1);                 % Lower bound for the input vector     
    UB = x0 + 50*ones(4,1);                 % Upper bound for the input vector
    OPTIONS = optimset('MaxFunEval',iterations ,'DiffMinChange', 1); % step size is one
    
    x = fmincon(ObjFcn, x0, A, B, Aeq, Beq, LB, UB, NONLCON, OPTIONS);       
    
    best_detect_param.smooth_weight = x(1);
    best_detect_param.smooth_var = x(2);
    best_detect_param.egt_weight = x(3);  
    best_detect_param.repulsion = x(4); 
    best_stats_param.comb_vector = [x(1) x(2) x(4) data_matrix(1,4) x(3)];
    best_stats_param = cluster_kernel_Viterbi_3D(sources, references, best_detect_param, best_stats_param);
    best_stats_param = compute_hit_ratios(best_stats_param);
    best_stats_param = compute_errors(best_stats_param);
end


function  [ mean_error ]  = Viterbi_handle_3D(sources, refs, best_detect_param, best_stats_obj, sw, sv, gt, repl)  
  best_detect_param.smooth_weight = sw;
  best_detect_param.smooth_var = sv;
  best_detect_param.egt_weight = gt;
  best_detect_param.repulsion = repl;
  best_stats_obj.comb_vector(1,1) = sw;
  best_stats_obj.comb_vector(1,2) = sv;
  best_stats_obj.comb_vector(1,3) = repl;
  best_stats_obj.comb_vector(1,5) = gt;
  stats_struct = cluster_kernel_Viterbi_3D(sources, refs, best_detect_param, best_stats_obj);
  mean_error = nanmean(stats_struct.total_diff(:));
  fprintf('Mean Error: %d \n', mean_error);
end


