function bubble_shape_scatter
% From: M. R. Vant, R. O. Ramseier, and V. Makios, "The complex-dielectric
% constant of sea ice at frequencies in the range 0.1-40 GHz," Journal of Applied
% Physics, vol. 49, no. 3, Mar. 1978.
% -- Specifically page 1266.  We have hard copies of this article.
%
% Vant cites: M. Kerker, The Scattering of Light and Other Electromagnetic
% Radiation, Academic, New York, 1969.

% Consider air bubble prolate spheroids of minor-axis 1 mm and major-axis 10 mm
a = 10e-3;
b = 1e-3;

physicalConstants;
e_ice = 1.78^2;

% Ratio of the complex refractive index of the particle to that of the medium
m = e_ice/1;

ecc = sqrt(1 - (b/a).^2);
if 1
  % E-field lines with major-axis
  Pa = 4*pi*(b./a).^2*(1./(2*ecc.^3)).*(-2*ecc + log((1+ecc)./(1-ecc)));
  a_prime_cubed = (4/3*a*b^2)*(m^2+2)*(4*pi + (m^2-1)*Pa).^-1;
else
  % E-field lines with one of minor-axes
  Pb = 4*pi*(b./a).^2*(1./(4*ecc.^3)).*(2*ecc/(b/a)^2 + log((1-ecc)./(1+ecc)));
  a_prime_cubed = (4/3*a*b^2)*(m^2+2)*(4*pi + (m^2-1)*Pb).^-1;
end

fprintf('Effective radius: %.3f mm\n',a_prime_cubed^(1/3)*1000);

% Effective volume per particle
V = 4/3*pi*a_prime_cubed;

% Number of paricles
N = logspace(3,7.3,101);
density_fraction = (1-4/3*pi*a*b^2*N);

% Scattering cross-section
f = [10e6 100e6 500e6 1000e6];
lambda = 3e8./f;
Csca = (24*pi^3*V^2)*lambda.^-4*abs((m^2-1)/(m^2+2))^2

if 0
hand1 = semilogx(N,10*log10(exp(-N*Csca(1)*1)),'k-');
hold on;
hand2 = semilogx(N,10*log10(exp(-N*Csca(2)*1)),'k-+');
hand3 = semilogx(N,10*log10(exp(-N*Csca(3)*1)),'k:');
hand4 = semilogx(N,10*log10(exp(-N*Csca(4)*1)),'k-o');
hold off;
xlabel('Number of Bubbles per m^3');
ylabel('Power Scattered (dB/m)');
legend([hand1 hand2 hand3 hand4],'10 Mhz','100 MHz','500 MHz','1 GHz');
title('Figure 10 in Robin et al, 1969 (using Vant''s equations)');
pause;
end

hand1 = semilogx(N,-10*log10(exp(-N*Csca(1)*1000)),'k-');
hold on;
hand2 = semilogx(N,-10*log10(exp(-N*Csca(2)*1000)),'k-+');
hand3 = semilogx(N,-10*log10(exp(-N*Csca(3)*1000)),'k:');
hand4 = semilogx(N,-10*log10(exp(-N*Csca(4)*1000)),'k-o');
hold off;
xlabel('Number of Bubbles per m^3');
ylabel('Effective Attenuation per km (dB/km)');
legend([hand1 hand2 hand3 hand4],'10 Mhz','100 MHz','500 MHz','1 GHz');
pause;

semilogx(N,density_fraction,'k-');
xlabel('Number of Bubbles per m^3');
ylabel('Fraction of Volume That is Air');
