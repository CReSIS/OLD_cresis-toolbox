function bubbles
% From: G. De Q. Robin, S. Evans, J.T. Bailey, "Interpretation of Radio Echo
% Sounding in Polar Ice Sheets," Philosophical Transactions of the Royal Society
% of Lonson. Series A, Mathematical and Physical Sciences, Vol. 265, No. 1166,
% (Dec. 18, 1969), 437-505.
% -- Specifically section 4. pages 447-454.  We have a pdf version of this
%    article.
%
% From: B.M. Smith, and S. Evans, "Radio echo sounding: absorption and
% scattering by water inclusion and ice lenses," Journal of Glaciology, Vol. 11,
% No. 61, 1972.
% -- Specifically pages 141-142.  We have hard copies of this article.

% Consider air bubble spheres of radius 1 mm (according to Smith, 1972 this
% is large)
b = 1e-3;

physicalConstants;
e_ice = 1.78^2;

% gamma is the polarizability.  The polarizability is the ratio of the induced
% dipole moment divided by applied electric field.  For a sphere of AIR, with
% radius b and placed in snow:
gamma = 4*pi*b^3*e0*e_ice*(1-e_ice)/(1+2*e_ice);

m = logspace(3,8,21);
density_fraction = (1-4/3*pi*b^3*m);

f = [10e6 100e6 500e6 1000e6];
w = 2*pi*f;
ss = 8/3*pi*10^-14*w.^4*gamma^2;

% How much energy is scattered in a 1m x 1m x 1m cube with the energy density
% directed from the top of the cube to the bottom of the cube.

hand1 = semilogx(m,10*log10(m*ss(1)),'k-');
hold on;
hand2 = semilogx(m,10*log10(m*ss(2)),'k-+');
hand3 = semilogx(m,10*log10(m*ss(3)),'k:');
hand4 = semilogx(m,10*log10(m*ss(4)),'k-o');
hold off;
xlabel('Number of Bubbles per m^3');
ylabel('Power Scattered (dB/m)');
legend([hand1 hand2 hand3 hand4],'10 Mhz','100 MHz','500 MHz','1 GHz');
title('Figure 10 in Robin et al, 1969');
pause;

hand1 = semilogx(m,-10*log10(1-m*ss(1)*1000),'k-');
hold on;
hand2 = semilogx(m,-10*log10(1-m*ss(2)*1000),'k-+');
hand3 = semilogx(m,-10*log10(1-m*ss(3)*1000),'k:');
hand4 = semilogx(m,-10*log10(1-m*ss(4)*1000),'k-o');
hold off;
xlabel('Number of Bubbles per m^3');
ylabel('Effective Attenuation per km (dB/km)');
legend([hand1 hand2 hand3 hand4],'10 Mhz','100 MHz','500 MHz','1 GHz');
pause;

semilogx(m,density_fraction);
xlabel('Number of Bubbles per m^3');
ylabel('Fraction of Volume That is Air');
