function er = brine(freq,temp,salinity)
% er = brine(freq,temp,salinity)
%
% Returns dielectric of brine according to Ulaby, Moore, Fung, Vol 3, Appendix E.
% Valid for 0 <= normality <= 3 which corresponds to 0 <= salinity <= 170
% parts per thousand. Temperature limits depend on the method used in
% salinity_to_volume_fraction.m and the assumption that the data for
% temperature > 0 Celcius can be extrapolated down... valid to -8 C at
% least.
%
% freq = frequency in Hz
% temp = temperature in Celcius
% salinity = salinity in kg/L or kg/kg (often given in parts per thousand,
%   such as 35 parts per thousand which would be 0.035 for this function)
%
% Brine is highly saline water. See salt_water.m for more information.
%
% Relating salinity to the normality of a solution: salinity_to_normality.m
%
% freq = logspace(9,11);
% er = brine(freq,-5,0.0856);
% figure(1); clf;
% h1 = loglog(freq/1e9,real(er),'k-');
% hold on;
% h2 = loglog(freq/1e9,-imag(er),'k:');
% hold off;
% legend([h1 h2], 'er_b''', 'er_b''''');
% xlabel('frequency (GHz)');
% ylabel('dielectric constant of liquid brine, er_b');
%
% Author: John Paden
%
% See also brine.m, sea_ice.m, rock.m, salt_water.m, soil.m, water.m
%   and salinity_to_normality.m, salinity_to_volume_fraction.m

% Get e0 (permittivity of free space)
physical_constants;

% Convert to parts per thousand
salinity = salinity * 1000;

% Convert from salinity to normality
Nb = salinity_to_normality(salinity);

%% Assymptotic dielectric of water
er_water_inf = 4.9;

%% First Term
eb0_water = 88.045 - 0.4147*temp + 6.295e-4*temp.^2 + 1.075e-5*temp.^3;
a1 = 1.0 - 0.255*Nb + 5.15e-2*Nb.^2 - 6.89e-3*Nb.^3;
er_b0 = eb0_water .* a1;

%% Second Term
tao_water = 1/(2*pi) * (1.1109e-10 - 3.824e-12*temp + 6.938e-14*temp.^2 - 5.096e-16*temp.^3);
b1 = 1.0 + 0.146e-2*temp.*Nb - 4.89e-2*Nb - 2.97e-2*Nb.^2 + 5.64e-3*Nb.^3;
tao_b = tao_water .* b1;

%% Third Term
sigma_b_25 = Nb .* (10.39 - 2.378*Nb + 0.683*Nb.^2 - 0.135*Nb.^3 + 1.01e-2*Nb.^4);
dtemp = 25-temp;
c1 = 1.0 - 1.96e-2*dtemp + 8.08e-5*dtemp.^2 - Nb.*dtemp ...
  .* (3.02e-5 + 3.92e-5*dtemp + Nb.*(1.72e-5 - 6.58e-6*dtemp));
sigma_b = sigma_b_25 .* c1;

%% Put it all together
er_real = er_water_inf + (er_b0 - er_water_inf) ./ (1 + (2*pi*freq*tao_b).^2);

er_imag = (2*pi*freq*tao_b) .* (er_b0 - er_water_inf) ./ (1 + (2*pi*freq*tao_b).^2) + sigma_b ./ (2*pi*freq.*e0) ;

er = er_real - j*er_imag;

return;



