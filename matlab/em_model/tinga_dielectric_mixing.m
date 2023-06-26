function er = tinga_dielectric_mixing(er_1,er_2,V2,a,b,theta)
% er = tinga_dielectric_mixing(er_1,er_2,V2,a,b,theta)
%
% Support function for sea_ice.m.
%
% See M. R. Vant, R. O. Ramseier, and V. Makios, "The complex dielectric
% constant of sea ice at frequencies in the range 0.1–40 GHz", Journal
% of Applied Physics, vol. 49, no. 1264 (1978); doi: 10.1063/1.325018
%    - OR -
% Ulaby, Moore, Fung, Vol 3, Appendix E.
%
% er_1 = dielectric constant of the host medium
% er_2 = dielectric constant of the inclusions
% V2 = volume fraction of the inclusions
% a = major axis of ellipsoidal inclusions
% b = minor axis of ellipsoidal inclusions (relative values a/b are all that
%     matter in the dielectric mixing model), b == c is assumed
% theta = angle to vertical of ellipsoidal inclusions
%
% a = 1;
% b = a/20;
% er_1 = 3.14;
% er_2 = 90-30j;
% V1 = 0.7;
% theta = 40/180*pi;
% er = tinga_dielectric_mixing(er_1,er_2,V1,a,b,theta)
% figure(1); clf;
% h1 = semilogx(freq/1e9,real(er),'k-');
% hold on;
% h2 = semilogx(freq/1e9,-imag(er),'k:');
% hold off;
% legend([h1 h2], 'er_m''', 'er_m''''');
% xlabel('frequency (GHz)');
% ylabel('dielectric constant of ice/brine mixture, er_m');
% xlim([1 20]);
%
% Author: John Paden
%
% See also brine.m, sea_ice.m, rock.m, salt_water.m, soil.m, water.m
%   and salinity_to_normality.m, salinity_to_volume_fraction.m

% V1 = volume fractions of the host medium
V1 = 1-V2;

%% Compute eccentricity of ellipsoids
c = b;
e = sqrt(1-(b./a).^2);

%% Compute n1 and n2 = depolarization coefficients, for ellipsoidal inclusions
if e == 0
  % For spherical inclusions: n1 = 1/3
  n1_a = 1/3;
  n1_b = n1_a;
%   n1_c = n1_b;
else
  n1_b = (b./c).^2 ./(4*e.^3) .* (2*e./(b./a).^2 + log((1-e)./(1+e)));
%   n1_c = n1_b;
  n1_a = (b./a).^2 ./ (2*e.^3) .* (-2*e + log((1+e)/(1-e)));
end
n2_a = n1_a;
n2_b = n1_b;

%% This mixing equation is valid as long as quasistatic assumption is met,
% i.e. scattering effects are minimal.
er_a = er_1 + er_1 .* V2/V1 .* (er_2 - er_1) ./ (-(V2./V1) .* n1_a .* (er_2-er_1) + n2_a.*(er_2-er_1) + er_1);
er_b = er_1 + er_1 .* V2/V1 .* (er_2 - er_1) ./ (-(V2./V1) .* n1_b .* (er_2-er_1) + n2_b.*(er_2-er_1) + er_1);

%% Account for ellipsoidal tilt
er_theta = (sin(theta).^2 ./ er_a + cos(theta).^2 ./ er_b).^-1;

er = 0.5*(er_b + er_theta);

return;

