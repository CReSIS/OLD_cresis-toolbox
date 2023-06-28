function H = boxcar(N)
% H = boxcar(N)
%
% Returns N-point boxcar window. Replaces Matlab's signal processing boxcar
% function.

H = ones(N,1);
