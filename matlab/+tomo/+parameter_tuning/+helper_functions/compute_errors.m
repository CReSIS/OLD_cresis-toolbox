function [ stats_struct ] = compute_errors( stats_struct )
%   compute the average rmse, average error, and median (of the median
%   error) for each set of parameters applied to the test (the set of frames we've choosed)

%   input: 
%   stats_struct: a structure used for statistics purposes (after
%   hyperparameters tests)
    stats_struct.rmse = sqrt(nanmean(stats_struct.total_diff(:).^2));
    stats_struct.mean_difference = nanmean(stats_struct.total_diff(:));
    stats_struct.median_difference = nanmedian(stats_struct.total_diff(:));
    
end
