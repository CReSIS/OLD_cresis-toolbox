function er = sea_ice(freq,temp,salinity,density,method)
% er = sea_ice(freq,temp,salinity,density)
%
% Returns dielectric of sea ice.
%
% See M. R. Vant, R. O. Ramseier, and V. Makios, "The complex dielectric
% constant of sea ice at frequencies in the range 0.1–40 GHz", Journal
% of Applied Physics, vol. 49, no. 1264 (1978); doi: 10.1063/1.325018
%    - OR -
% Ulaby, Moore, Fung, Vol 3, Appendix E.
%
% freq = frequency in Hz
% temp = temperature in Celcius
% salinity = salinity in kg/L or kg/kg (often given in parts per thousand,
%   such as 35 parts per thousand which would be 0.035 for this function)
% density = density of sea ice in g/cm^3
% method = 1, 2, or 3 (default is 2)
%   See code for descriptions.
%
% freq = logspace(9,10.3);
% er = sea_ice(freq,-10,0.0008,0.95,2);
% figure(1); clf;
% h1 = semilogx(freq/1e9,real(er),'k-');
% hold on;
% h2 = semilogx(freq/1e9,-imag(er),'k:');
% hold off;
% legend([h1 h2], 'er_b''', 'er_b''''');
% xlabel('frequency (GHz)');
% ylabel('dielectric constant of sea ice, er_b');
% xlim([1 20]);
%
% Some examples referenced from Ulaby et al:
%   FY Frazil: salinity = 4.4 g/L, density = 0.836 g/cm^3
%   FY Frazil: salinity = 3.2 g/L, density = 0.836 g/cm^3
%   FY Columnar: salinity = 3.2 g/L, density = 0.878 g/cm^3
%   FY Columnar: salinity = 4.6 g/L, density = 0.896 g/cm^3
%   MY: salinity = 0.61 g/L, density = 0.771 g/cm^3
%   MY: salinity = 0.7 g/L, density = 0.770 g/cm^3
%
% Author: John Paden
%
% See also brine.m, sea_ice.m, rock.m, salt_water.m, soil.m, water.m
%   and salinity_to_normality.m, salinity_to_volume_fraction.m

%% Input Checking
if ~exist('method','var')
  method = 2;
end

%% Compute volume fraction of air
%   rho_seaice = 0.92-0.96 for first-year ice
%   rho_seaice = 0.7 for multiyear ice
rho_seaice_nonporous = 0.926;
vf_air = 1 - density/rho_seaice_nonporous;
if vf_air < 0
  vf_air = 0;
end
vf_notair = 1-vf_air;

%% Compute Sb in parts per thousand (g/L)
Sb = zeros(size(temp));
n = find(temp>-8.2);
Sb(n) = 1.725 - 18.756*temp(n) - 0.3964*temp(n).^2;
n = find(temp <= -8.2 & temp > -22.9);
Sb(n) = 57.041 - 9.929*temp(n) - 0.16204*temp(n).^2 - 0.002396*temp(n).^3;
n = find(temp <= -22.9 & temp > -36.8);
Sb(n) = 242.94 + 1.5299*temp(n) + 0.0429*temp(n).^2;
n = find(temp <= -36.8);
Sb(n) = 508.18 + 14.535*temp(n) + 0.2018*temp(n).^2;

%% Compute dielectric of brine
er_b = brine(freq,temp,Sb/1000);

%% Compute volume fraction of brine
salinity = salinity * 1000;
vf_brineinice = salinity * salinity_to_volume_fraction(temp);
vf_host = 1-vf_air-vf_notair*vf_brineinice;
vf_brine = 1-vf_host-vf_air;

%% Compute er of ice
rho_ice = 0.917; % This is solid ice for the iceCond equation
er_ice = iceCond(temp+273.15,rho_ice,freq,0,22,-20+273.15);

%% Compute er of air
er_air = 1;

%% Apply mixing formula to determine er of sea ice
if method == 1
  %% Empirical (no air bubble component)
  % er = er_ice + 3*vf_brine .* er .* (er_b - er_ice) ./ (2*er + er_b);
  er = real(er_ice) / (1 - 3*vf_brine) + j*vf_brine*imag(er_b);
  
elseif method == 2
  %% W. R. Tinga, W. A. G. Voss, and D. F. Blossey, J. Appl. Phys. 44,
  % 3897 (1973).
  
%   er_ice_brine = tinga_dielectric_mixing(er_ice,er_b,vf_brine,20,1,89/180*pi);
  er_ice_brine = tinga_dielectric_mixing(er_ice,er_b,vf_brine,20,1,45/180*pi);
  er = er_ice_brine;
%   er = tinga_dielectric_mixing(er_ice_brine,er_air,vf_air,1,1,0/180*pi);
  
elseif method == 3
  %% Gloersen and Larabee (1981) referenced from
  % Ulaby, Moore, Fung, Vol 3, Appendix E, pp. 2051.
  
  
end

return;

% This does not match Figure E.20... should go up near 0 C.
temp = linspace(-60,-1,101);
freq = 10e9;
salinity = 0.0044;
density = 0.836;

er = sea_ice(freq,temp,salinity,density,2);
figure(1); clf;
h1 = plot(-temp,real(er),'k-');
xlabel('temperature (C)');
ylabel('permittivity of sea ice, er_si''');

% Vant fig 5
temp = -5.1;
freq = logspace(8,10.3);
salinity = 0.0051;
density = 0.91;

er = sea_ice(freq,temp,salinity,density,2);
figure(1); clf;
h1 = semilogx(freq/1e9,real(er),'k-');
xlabel('frequency (GHz)');
ylabel('permittivity of sea ice, er_si''');

% Vant fig 7
temp = -5.20;
freq = logspace(8,10.3);
salinity = 0.0105;
density = 0.91;

er = sea_ice(freq,temp,salinity,density,2);
figure(1); clf;
h1 = semilogx(freq/1e9,real(er),'k-');
xlabel('frequency (GHz)');
ylabel('permittivity of sea ice, er_si''');

% Vant fig 9
temp = -3.99;
freq = logspace(8,10.3);
salinity = 0.0075;
density = 0.91;

er = sea_ice(freq,temp,salinity,density,2);
figure(1); clf;
h1 = semilogx(freq/1e9,real(er),'k-');
xlabel('frequency (GHz)');
ylabel('permittivity of sea ice, er_si''');

% Vant fig 11
temp = -5.10;
freq = logspace(8,10.3);
salinity = 0.0051;
density = 0.91;

er = sea_ice(freq,temp,salinity,density,2);
figure(1); clf;
h1 = loglog(freq/1e9,-imag(er),'k-');
xlabel('frequency (GHz)');
ylabel('permittivity of sea ice, er_si''''');
xlim([0.1 5]);
ylim([0.01 10]);

