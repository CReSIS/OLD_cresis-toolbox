function [kx dkx] = gen_kx(x, Nx)
% [kx dkx] = gen_kx(x, Nx)
%
% Generates wavenumber, kx, axis from x. Assumes x is uniformly sampled.
% Nx is optional, positive integer, that specifies the number of points
% to generate rather than the default number of points = length(x).

if ~exist('Nx','var') || isempty(Nx)
  Nx     = length(x);
end

dx     = (x(end) - x(1)) / (length(x)-1);
X      = Nx*dx;
dkx    = 2*pi/X;
kx     = dkx * ifftshift(-floor(Nx/2) : floor((Nx-1)/2));

return