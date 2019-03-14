function H = chebwin(N, R)
% H = chebwin(N, R)
%
% Returns N-point chebwin with sidelobe attenuation of R. Replaces Matlab's
% signal processing toolbox chebwin function.

N1 = N - 1;

alpha = cosh(acosh(10.^(R/20))/N1);

m = (0:N1-1).';

Am = alpha*cos(pi*m/N1);

H = zeros(N1,1);
mask = abs(Am)<=1;
H(mask) = cos(N1*acos(Am(mask)));
H(~mask) = cosh(N1*acosh(Am(~mask)));
H(2:2:end) = -1*H(2:2:end);

H = real(ifft(H));
H(1) = H(1)/2;
H(end+1) = H(1);

H = H ./ max(H);
