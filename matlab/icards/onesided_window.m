function [ window_matrix ] = onesided_window( td_window,td_window_side,data_size )
% this function is to create a one-sieded window to avoid big value at the begining or
% the end of records 
% td_window=a full window we use
% td_window_side='begin' or 'end'

window_single=td_window;
if strcmpi(td_window_side,'begin')
  window_onesided=[window_single(1:floor(length(window_single)/2));ones(data_size(1,1)-floor(length(window_single)/2),1)];
else
  window_onesided=[ones(data_size(1,1)-floor(length(window_single)/2),1);window_single(end-floor(length(window_single)/2)+1:end)];
end

window_matrix=repmat(window_onesided,1,data_size(1,2));

end

