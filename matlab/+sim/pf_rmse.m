function [rmse_results] =  pf_rmse(param)
% This function is called from the main particle filter sumulator script.
% It calculates the RMSE of the estimated DoAs before and after smoothing
% is applied.
%
% Author: Mohanad Al-Ibadi
%
%%                               Define parameters
% =========================================================================
doa_error          = param.est_doa_error; % Cell array of the DoA errors for all trials
Nt                 = param.Nt; % Total number of range-bins
Nx                 = param.Nx; % total number of range-lines
rbins              = param.rbins; % Range-bins to be processed
rlines             = param.rlines; % Range-lines to be processed

%%                               RMSE calculation
% =========================================================================
%% RMSE before smoothing is applied
% ---------------------------------
% Remove error outliers using 3 sigma rule (error>3 std from the mean is
% ignored). Outliers occur in cases where the number of particles is small
% or the Gaussian posterior doesn't cover the entire support of the actual
% posterior.
tmp_doa_error = cat(3,doa_error{:});
% tmp_doa_error(isnan(tmp_doa_error)) = 0;
% Setting threshold = inf ignores its effect (no thresholding). The
% standard value is threshold = 3.
threshold = 3* std(tmp_doa_error(~isnan(tmp_doa_error(:))));
% threshold = inf* std(tmp_doa_error(~isnan(tmp_doa_error(:))));
tmp_doa_error(tmp_doa_error>threshold) = 0;

% RMSE over runs for each doa and range-bin 
for rbin = 1:size(tmp_doa_error,1)
  xx = squeeze(tmp_doa_error(rbin,:,:));
  est_error_doa(rbin,:) = sqrt(nanmean(xx.^2,2));
end
% Average RMSE over all DOAs
rmse = mean(est_error_doa,2);
rmse(rmse==0) = NaN;
% rmse = NaN(size(tmp_doa_error));
rmse_other = sqrt(nanmean(tmp_doa_error.^2,3));
rmse_other(rmse_other==0) = NaN;
avg_rline_rmse = nanmean(rmse_other,1); % Average RMSE over range-bins
avg_rbin_rmse = nanmean(rmse_other,2); % Average RMSE over range-lines

% Average RMSE over range-bins
% avg_rline_rmse = nansum(rmse,1);
% for rline = 1:size(rmse,2)
%   avg_rline_rmse(rline) = avg_rline_rmse(rline)./nnz(rmse(:,rline)) ;
% end
% avg_rline_rmse(avg_rline_rmse==0) = NaN;

% Average RMSE over range-lines
% avg_rbin_rmse = nansum(rmse,2);
% for rbin =1:size(rmse,1)
%   avg_rbin_rmse(rbin) = avg_rbin_rmse(rbin)/nnz(rmse(rbin,:)) ;
% end
% avg_rbin_rmse(avg_rbin_rmse==0) = NaN;
% 
% rmse(rmse==0) = NaN;

%% RMSE after smoothing is applied
% --------------------------------
if isfield(param,'smoothing_method') && ~isempty(param.smoothing_method)
  smoothed_doa_error = param.smoothed_doa_error; % Cell array of the DoA errors after smoothing for all trials
  
  tmp_smoothed_doa_error = cat(3,smoothed_doa_error{:});
  tmp_smoothed_doa_error(isnan(tmp_smoothed_doa_error)) = 0;
  threshold = 3 *std(tmp_smoothed_doa_error(~isnan(tmp_smoothed_doa_error(:))));
%   threshold = 3 *std(tmp_smoothed_doa_error(:));
  tmp_smoothed_doa_error(tmp_smoothed_doa_error>threshold) = 0;
  
  % MSE over all runs
%   rmse_smoothed_doa = NaN(Nt,Nx);
  rmse_smoothed_doa = sqrt(nanmean(tmp_smoothed_doa_error.^2,3));
  
  % Average RMSE over range-bins
  avg_rline_rmse_smoothed_doa = nansum(rmse_smoothed_doa,1);
  for rline = 1:size(rmse_smoothed_doa,2)
    avg_rline_rmse_smoothed_doa(rline) = avg_rline_rmse_smoothed_doa(rline)./nnz(rmse_smoothed_doa(:,rline)) ;
  end
  avg_rline_rmse_smoothed_doa(avg_rline_rmse_smoothed_doa==0) = NaN;
  
  % Average RMSE over range-lines
  avg_rbin_rmse_smoothed_doa = nansum(rmse_smoothed_doa,2);
  for rbin =1:size(rmse_smoothed_doa,1)
    avg_rbin_rmse_smoothed_doa(rbin) = avg_rbin_rmse_smoothed_doa(rbin)/nnz(rmse_smoothed_doa(rbin,:)) ;
  end
  avg_rbin_rmse_smoothed_doa(avg_rbin_rmse_smoothed_doa==0) = NaN;
  
  rmse_smoothed_doa(rmse_smoothed_doa==0) = NaN;
end

%%                                 Results
% =========================================================================
rmse_results.rmse           = rmse;
rmse_results.avg_rline_rmse = avg_rline_rmse;
rmse_results.avg_rbin_rmse  = avg_rbin_rmse;

if isfield(param,'smoothing_method') && ~isempty(param.smoothing_method)
  rmse_results.rmse_smoothed_doa           = rmse_smoothed_doa;
  rmse_results.avg_rline_rmse_smoothed_doa = avg_rline_rmse_smoothed_doa;
  rmse_results.avg_rbin_rmse_smoothed_doa  = avg_rbin_rmse_smoothed_doa;
end

return