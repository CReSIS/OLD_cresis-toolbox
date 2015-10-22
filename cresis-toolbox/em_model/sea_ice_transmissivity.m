% function T = sea_ice_transmissivity(scatter_case)
fprintf('==================================================\n');
fprintf('Sea Ice Transmissivity\n');
fprintf('  Checking Vant J. Applied Physics 1978 pg 1266 results\n');

if 1
  scatter_case = 2;
end

% Ellipsoid dimensions, a = major axis, b=c are minor axes
if scatter_case == 1
  % Brine in first-year ice
  a = 5e-3;
  b = 0.025e-3;
  c = b;
  
elseif scatter_case == 2
  % Air bubbles in first-year ice
  a = 0.5e-3;
  b = 0.5e-3;
  c = b;
  
elseif scatter_case == 3
  % Air bubbles in first-year ice
  a = 5e-3;
  b = 5e-3;
  c = b;
end

% e = eccentricity
e = sqrt(1 - (b./a)^2);

% Dielectric of inclusions (brine or air)
if scatter_case == 1
  er_inc = 70-20j;
elseif scatter_case > 1
  er_inc = 1;
end

% Dielectric of medium (ice)
er_medium = 3.14;

% m = ratio of the complex refractive index of the particle to that of the
%     surrounding medium
%     n is the index of refraction
%     For real and complex dielectrics: n^2 = er (Ulaby, Moore, Fung Vol 1
%     pg 67)
m = abs(sqrt(er_inc/er_medium)); % ADDED ABSOLUTE VALUE!

if e > 0
  % Ellipsoidal Inclusions
  
  % Depolarization factor
  Pb = 4*pi*(b./a).^2 / (4*e.^3) * (2*e*(b./a).^-2 + log( (1-e)/(1+e) ));
  Pc = Pb;
  Pa = 4*pi*(b./a).^2 / (2*e.^3) * (-2*e + log( (1+e)/(1-e) ));

  % V = volume of scattering
  V = 4/3*a*b*c*pi;

  % V_eff is the effective particle volume for scattering
  V_eff = (4/3*a.*b.*c) * (m.^2 + 2) * (4*pi + (m.^2 - 1) * Pb).^-1;
else
  % Spherical Inclusions
  
  % V = volume of scattering
  V = 4/3*a*b*c*pi;
  
  % V_eff is the effective particle volume for scattering
  V_eff = V;
end

lambda = logspace(log10(0.001),log10(0.2),101);
% lambda = 10e-2;

% Csca = scattering cross section
Csca = (24*pi^3*V_eff.^2) * lambda.^-4 * abs((m.^2 - 1) / (m.^2 + 2)).^2;

% T is the magnitude of the power-transmission coefficient
% prop_len is the propagation length
% N is the number of particles per unit volume

% Length of sample
prop_len = 0.08;

% VF = volume fraction
if scatter_case == 1
  VF = 0.11; % volume fraction for brine in first-year ice
elseif scatter_case == 2
  VF = 0.015; % volume fraction for air in first-year ice
elseif scatter_case == 3
  VF = 0.244; % volume fraction for air in multiyear ice
end

N = VF/V;

T = exp(-N*Csca*prop_len);

fprintf('  er inclusions: %.2f + %.2fj\n', real(er_inc), imag(er_inc));

Transmissivity = abs(T).^2;

figure(1); clf;
plot(lambda*1000, lp(Transmissivity));
ylim([-1 0]);
grid on;
if scatter_case == 1
  xlim([0.001 0.010]*1000);
elseif scatter_case == 2
  xlim([0.001 0.020]*1000);
elseif scatter_case == 3
  xlim([0.01 0.20]*1000);
end
xlabel('wavelength (mm)');
ylabel('transmissivity (dB)');

return




