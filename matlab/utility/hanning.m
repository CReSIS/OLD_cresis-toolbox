function H = hanning(N)
% H = hanning(N)
%
% Returns N-point Hanning window. Replaces signal processing toolbox
% function.

if ~rem(N,2)
  % Even length
  H = .5*(1 - cos(2*pi*(1:N/2)'/(N+1)));
  H = [H; H(end:-1:1)];
else
  % Odd length
  H = .5*(1 - cos(2*pi*(1:(N+1)/2)'/(N+1)));
  H = [H; H(end-1:-1:1)];
end

