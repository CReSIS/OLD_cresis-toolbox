function L = mle_compute_cost(A,Rxx)
% L = mle_compute_cost(A, Rxx)
% 
% This function is meant to be replace by a mex function,
% mle_compute_cost.mexa64.
%
% Inputs:
%   A = Nc x Nsig matrix of steering vectors,
%   Rxx = Nc x Nc sample covariance matrix
% 
% Outputs:
%   L = real valued cost of a Theta solution leading to the A matrix
%       based on the likelihood function
% 
% See also:  array_proc.m, mle_initialize.m, mle_cost_function.m
% =========================================================================


Pa  = A * inv(A'*A) * A';
L   = abs(trace(Pa*Rxx));

end