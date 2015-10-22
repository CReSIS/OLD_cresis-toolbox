function [C,lags] = xcorr_noedge(A,B,MAXLAG)
% [C,lags] = xcorr_noedge(A,B,MAXLAG)
%
% Removes edge effects by taking a subset of the total dataset.
% A and B must be vectors.
% 
% Useful when you have two long vectors that you are just correlating 
% a short portion of.  MAY NOT WORK CORRECTLY... NEEDS WORK!
%
% Author: ?

lags = -MAXLAG:MAXLAG;
C = zeros(size(lags));
for lag_idx = 1:length(lags)
  lag = lags(lag_idx);
  C(lag_idx) = dot(A(1+MAXLAG:end-MAXLAG),B(1+MAXLAG+lag:end-MAXLAG+lag));
end

return;
