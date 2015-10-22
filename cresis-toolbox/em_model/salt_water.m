function er = salt_water(freq,temp,salinity)
% er = salt_water(freq,temp,salinity)
%
% Returns dielectric of salt water according to Ulaby, Moore, Fung, Vol 3, Appendix E.
% Valid for 0 <= salinity <= 0.04, 0 <= T <= 40.
%
% freq = frequency in Hz
% temp = temperature in Celcius
% salinity = salinity in kg/L or kg/kg (often given in parts per thousand,
%   such as 35 parts per thousand which would be 0.035 for this function)
%
% Salinity is the total mass of solid salt in grams dissolved in one
% kilogram of solution. Usually measured in parts per thousand on a weight
% basis.
%
% Oceans have an average salinity of 32.54 parts per thousand.  So
% salinity = 0.03254.
%
% Relating salinity to the normality of a solution: salinity_to_normality.m
%
% freq = logspace(9,12);
% er = salt_water(freq,0,0.03254);
% figure(1); clf;
% loglog(freq/1e9,real(er),'k-');
% xlabel('frequency (GHz)');
% ylabel('loss factor of water, er_w''');
% figure(2); clf;
% loglog(freq/1e9,-imag(er),'k-');
% xlabel('frequency (GHz)');
% ylabel('permittivity of water, er_w''''');
%
% Author: John Paden
%
% See also brine.m, sea_ice.m, rock.m, salt_water.m, soil.m, water.m
%   and salinity_to_normality.m

% Get e0 (permittivity of free space)
physical_constants;

% Convert to parts per thousand
salinity = salinity * 1000;

%% Assymptotic dielectric of water
er_water_inf = 4.9;

%% First Term
eb0_water = 87.134 - 1.949e-1*temp - 1.276e-2*temp.^2 + 2.491e-4.^3;
a1 = 1.0 + 1.613e-5*temp.*salinity - 3.656e-3*salinity ...
  + 3.210e-5*salinity.^2 - 4.232e-7*salinity.^3;
er_b0 = eb0_water .* a1;

%% Second Term
tao_water = 1/(2*pi) * (1.1109e-10 - 3.824e-12*temp + 6.938e-14*temp.^2 - 5.096e-16*temp.^3);
b1 = 1.0 + 2.282e-5*temp.*salinity - 7.638e-4*salinity - 7.760e-6*salinity.^2 + 1.105e-8*salinity.^3;
tao_b = tao_water .* b1;

%% Third Term
sigma_b_25 = salinity*(0.18252 - 1.4619e-3*salinity + 2.093e-5*salinity.^2 ...
  - 1.282e-7*salinity.^3);
dtemp = 25-temp;
c1 = dtemp.*(2.033e-2 + 1.266e-4*dtemp + 2.464e-6*dtemp.^2 ...
  - salinity*(1.849e-5 - 2.551e-7*dtemp + 2.551e-8*dtemp.^2));
sigma_b = sigma_b_25 .* c1;

%% Put it all together
er_real = er_water_inf + (er_b0 - er_water_inf) ./ (1 + (2*pi*freq*tao_b).^2);

er_imag = (2*pi*freq*tao_b) .* (er_b0 - er_water_inf) ./ (1 + (2*pi*freq*tao_b).^2) + sigma_b ./ (2*pi*freq*e0);

er = er_real - j*er_imag;

return;
