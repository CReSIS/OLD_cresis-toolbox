function [data_matrix] = result_make_data_matrix(cluster_result)
%Summary of this function goes here
%   format the data after running the tests on the cluster

  %% preprocessed the data
  reverseStr = '';
  for idx = 1:length(cluster_result)
    cluster_result{idx}.argsout{1} = compute_hit_ratios(cluster_result{idx}.argsout{1});
    cluster_result{idx}.argsout{1} = compute_errors(cluster_result{idx}.argsout{1});

    % Display the progress
    percentDone = 100 * idx / length(cluster_result);
    msg = sprintf('Preprocessing the data -- percent done: %3.1f\n', percentDone); 
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
  end
  
  data_matrix = make_big_matrix(cluster_result);
end

function [data_frame_result] = make_big_matrix( cluster_result_processed ) 
  data_frame_result = [];
  % make a numerical matrix for the ease of data visualization purposes

  for idx = 1:length(cluster_result_processed)
    row = [cluster_result_processed{idx}.argsout{1}.comb_vector , ...
      cluster_result_processed{idx}.argsout{1}.hit_ratios, ...
      cluster_result_processed{idx}.argsout{1}.rmse, ...
      cluster_result_processed{idx}.argsout{1}.mean_difference, ...
      cluster_result_processed{idx}.argsout{1}.median_difference];
    data_frame_result = [data_frame_result; row];
  end
end

