function er = water(T,freq)
% er = water(T,freq)
%
% T = temperature in Kelvin
% freq = frequency in Hz
%
% e0*er = e0*(e' - j*e'')
%
% Information from pg 2020-2022 of Ulaby, Moore, and Fung, 
% "Microwave Remote Sensing," vol 3.
%
% Author: John Paden
%
% See also brine.m, sea_ice.m, rock.m, salt_water.m, soil.m, water.m
%   and salinity_to_normality.m, salinity_to_volume_fraction.m

T = T - 273.15;

einf = 4.9;
tao = 1/(2*pi) * (1.1109e-10 - 3.824e-12*T + 6.938e-14*T.^2 - 5.096e-16*T.^3);
estatic = 88.045 - 0.4147*T + 6.295e-4*T.^2 + 1.075e-5*T.^3;
er = einf + (estatic-einf) ./ (1 + j*2*pi*freq*tao);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1 (recreates figure E.1 in Ulaby, Moore, and Fung)
clf; clear all;
freq = logspace(9,12,100);
w = 2*pi*freq;
physicalConstants;
sigma = 0;
er1 = water(273.15+0,freq);
er2 = water(273.15+20,freq);
subplot(2,1,1);
loglog(freq/1e9,real(er1),'r');
hold on;
loglog(freq/1e9,real(er2),'k');
hold off;
ylabel('Real permittivity of water \epsilon_w''');
title('Black = 20 C\circ, Red = 0 C\circ');
subplot(2,1,2);
% loglog(freq/1e9,-imag(er1)./real(er1),'r');
loglog(freq/1e9,-imag(er1),'r');
hold on;
% loglog(freq/1e9,-imag(er2)./real(er2),'k');
loglog(freq/1e9,-imag(er2),'k');
hold off;
ylabel('Imag permittivity of water \epsilon_w''''');
xlabel('Freq');

% Example 2
% reflection coefficient for 110-500 MHz from air to water
% water is -2 centigrade
clf; clear all;
freq = linspace(110e6,500e6,20);
w = 2*pi*freq;
physicalConstants;
sigma = 0;
er = water(273.15-2,freq);
atten = real(j*sqrt(-j*w.*u0.*(sigma+j*w.*er*e0)));
plot(freq,20*log10(exp(-atten*1)));
e = er*e0;
eta_air = sqrt(u0./e0);
eta_water = sqrt(u0./e);
refl = (eta_water - eta_air) ./ (eta_water + eta_air);
subplot(2,1,1);
plot(freq,20*log10(abs(refl)));
subplot(2,1,2);
plot(freq,phase(refl));

% Example 2
% reflection coefficient for 110-500 MHz from air to ice
% water is -2 centigrade
clf; clear all;
freq = linspace(110e6,500e6,20);
w = 2*pi*freq;
physicalConstants;
sigma = 0;
e_water = e0*water(273.15-2,freq);
e_ice = e0*ice(273.15-2,0.917,freq,2.5e-6);
eta_water = sqrt(u0./e_water);
eta_ice = sqrt(u0./e_ice);
refl = (eta_ice - eta_water) ./ (eta_ice + eta_water);
subplot(2,1,1);
plot(freq,20*log10(abs(refl)));
subplot(2,1,2);
plot(freq,phase(refl));




