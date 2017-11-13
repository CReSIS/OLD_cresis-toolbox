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

%appending zeroes after true directions so that the error direction is
%going to be accounted in RMSE completely
param.monte.DOA = [param.monte.DOA zeros(1,param.Nc-1-length(param.monte.DOA))];

RMSE = zeros(length(param.method.list),size(results.theta_est{param.method.list(1)},2), ...
  size(results.theta_est{param.method.list(1)},3));

for method_idx = 1:length(param.method.list)
  method = param.method.list(method_idx);
  for test_idx = 1:size(results.theta_est{method},2)
    
    
    %for one source case DOA estimator may give 0 or 20. it is not fair if
    %ans is 20 and we calculate RMSE w.r.t true doa at 0 deg.
    %TRY ADDRESSING THE ABOVE ISSUE
    
    for src_idx = 1:size(results.theta_est{method},3)
      RMSE(method_idx,test_idx,src_idx) ...
        = sqrt(mean(abs(results.theta_est{method}(:,test_idx,src_idx)*180/pi - ...
        param.monte.DOA(test_idx,src_idx)).^2));
    end
    
  end
end

return


