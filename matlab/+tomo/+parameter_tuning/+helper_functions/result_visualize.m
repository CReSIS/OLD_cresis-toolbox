function [] = result_visualize( data_matrix, test, algo_name, data_type, best_point)
%DATA_STATISTICS Summary of this function goes here
%   statistical analysis and data visualization on the processed data
  
  [~,rank_indices] = sort(data_matrix(:,end-2), 'ascend'); % rank by rmse in ascending order 
  data_matrix = data_matrix(rank_indices,:);

  %% rank the result  
  rank_result(data_matrix, algo_name, data_type); 
  
  %% grid result plot (for grid tests only)
  if strcmp(test, 'grid') && strcmp(data_type, '3D') && strcmp(algo_name, 'Viterbi')
    grid_plots(data_matrix, best_point);  
  end
end



function [] = grid_plots( data_matrix, best_point )

  best_smoothness_weight = data_matrix(1,1);
  best_smoothness_variance = data_matrix(1,2);
  best_repulsion = data_matrix(1,3);
  best_ice_b_threshold = data_matrix(1,4);
  best_egt_weight = data_matrix(1,5);

  best_of_best = [best_smoothness_weight, ...
    best_smoothness_variance, ...
    best_repulsion, ...
    best_ice_b_threshold, ...
    best_egt_weight];


  Xlabel = {
    'Smoothness Weight', ...
    'Smoothness Variance', ...
    'Repulsion', ...
    'Ice Distance Threshold', ...
    'Ground Truth Weight'};
  
   Ylabel = {
     'percentage_correct', ...
     '% error less than 5', ...
     '% error less than 10',...
     '% error less than 15',...
     '% error less than 20',...
     '% error less than 25',...
     'RMSE', 'Mean Error', 'Median Error'};

  array_index = [1 2 3 4 5];
  
%   optimal_params = [];
  
  for f = 1:9 % 8 measures
    figure(f);
    for i = 1:5 % 5 parameters
      cur_indices = setdiff(array_index, i);
      indices = intersect(intersect(intersect(...
        find(data_matrix(:,cur_indices(1))== best_of_best(cur_indices(1))), ...
        find(data_matrix(:,cur_indices(2)) == best_of_best(cur_indices(2)))), ...
        find(data_matrix(:,cur_indices(3)) == best_of_best(cur_indices(3)))), ...
        find(data_matrix(:,cur_indices(4)) == best_of_best(cur_indices(4))));
      subplot(3,3,i); plot(data_matrix(indices,i), data_matrix(indices,f+5), 'bo', 'MarkerEdgeColor','k', ...
  'MarkerFaceColor','b', 'MarkerSize',7); 
      hold on;  plot(best_point(1,i), best_point(1,f+5), 'diamond','MarkerEdgeColor','k', ...
  'MarkerFaceColor','r', 'MarkerSize',7);
      
      xlabel(Xlabel{i});
      ylabel(Ylabel{f});            
      
%       if f == 7 && (i == 1 || i ==2 || i ==5) % (convex) take sw, sv, and egt to account for now       
%         xd = data_matrix(indices,i);
%         yd = data_matrix(indices,f+5);
%         p = polyfit(xd, yd, 4);
% 
%         optimal_params = [optimal_params poly_min_pt(p)];
%       end      
    end    
    hSub = subplot(3,3,6); plot(1, nan, 'bo', 'MarkerEdgeColor','k', 'MarkerFaceColor','b', 'MarkerSize',7); hold on;
    plot(1, nan, 'diamond','MarkerEdgeColor','k', 'MarkerFaceColor','r', 'MarkerSize',7);
    set(hSub, 'Visible', 'off');    
    legend(hSub, 'grid search', 'fmincon', 'Location','east');
  end
end

function [min_pt] = poly_min_pt(polynomial)

    d1p = polyder(polynomial);                  % First Derivative
    d2p = polyder(d1p);                         % Second Derivative
    ips = roots(d1p);                           % Inflection Points
    xtr = polyval(d2p, ips);                    % Evaluate 'd2p' at 'ips'
        
    minpts = ips((xtr > 0) & (imag(xtr)==0));   % Find Minima    
    min_pt = min(minpts);   
    
    
%     keyboard
%     x = linspace(0,300);
%     ep = polyval(polynomial,x);
%     figure(1)
%     plot(x, ep, '-r')
%     hold on
%     plot(minpts, polyval(polynomial,minpts), 'bp')
%     hold off
%     grid
end




