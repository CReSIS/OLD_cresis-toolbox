function RMSE = doa_rmse(param,results)
% RMSE = sim.doa_rmse(param,results)
%
% Compute root mean squared error (RMSE) for each source from sim.doa
% param and results.
%
% RMSE(method,test index,source index)
%
% Author: John Paden, Theresa Stumpf
%
% See also: doa.m, doa_example.m
method_idx = 1:length(param.method.list);
    
RMSE = zeros(length(method_idx),size(results.theta_est{method_idx(1)},2), ...
  size(results.theta_est{method_idx(1)},3));

for method_idx = 1:length(param.method.list)
  for test_idx = 1:size(results.theta_est{method_idx},2)
    for src_idx = 1:size(results.theta_est{method_idx},3)
      Err = abs(results.theta_est{method_idx}(:,test_idx,src_idx)*180/pi - param.monte.DOA(test_idx,src_idx));
      sigma = std(Err(~isnan(Err)));
      Err(Err > inf*sigma) = 0;
      RMSE(method_idx,test_idx,src_idx) = sqrt(mean(abs(Err).^2));
    end
  end
end
