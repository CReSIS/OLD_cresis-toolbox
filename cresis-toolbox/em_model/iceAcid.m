function er = ice(T,rho,freq,acidity,eV)
% er = ice(T,rho,freq,acidity,eV)
%
% T = temperature in Kelvin
% rho = density of ice/snow (g/cm^3)
% freq = frequency in Hz
% acidity = molarity of acid concentration (moles of solute to liters of solution)
% eV = activation energy in electron Volts (typical: 0.22 eV)
%
% T, rho, and freq are N-dim matrices of the same size.  Except for T, if the
% parameter doesn't change, then it can be input as a scalar.
% Basically, you just have to satisfy the point-wise operators (e.g .*, .^, ./).
%
% e0*er = e0*(e' - j*e'')

if (ispc)
  ice_prop = load('P:\prism\radar\radarSimulator\profiles\matsuoka_ice_prop.txt');
else
  ice_prop = load('/projects/prism/radar/radarSimulator/profiles/matsuoka_ice_prop.txt');
end

A = interp1(ice_prop(:,1),ice_prop(:,2)/1e4,T,'spline','extrap');
B = interp1(ice_prop(:,1),ice_prop(:,4)/1e5,T,'spline','extrap');
C = interp1(ice_prop(:,1),ice_prop(:,6),T,'spline','extrap');

% Real part from Matzler 1987 and ref. in Fujita 2000
% Imag part from Matzler 1987 and ref. in Fujita 2000, Matsuoka 1997
% er = (3.1884+0.00091*(T-273.15)) + (log10(freq)-7.18)*-0.012 ...
%    -j*(A./(freq/1e9) + B.*(freq/1e9).^C);
er = (3.1884+0.00091*(T-273.15)) -j*(A./(freq/1e9) + B.*(freq/1e9).^C);

if (ispc)
  acid_prop = load('P:\prism\radar\radarSimulator\profiles\matsuoka_acid_prop.txt');
else
  acid_prop = load('/projects/prism/radar/radarSimulator/profiles/matsuoka_acid_prop.txt');
end
A = interp1(acid_prop(:,1),acid_prop(:,2),T-273.15,'spline','extrap');
B = interp1(acid_prop(:,1),acid_prop(:,4),T-273.15,'spline','extrap');
er = er + acidity.*10.^A.*freq.^B;
% Assume molar conductivity of ice to be 3.3 Siemans per meter
% per molarity (data ref Fujita).  This is at -20 C, so adjust using
% Arrhenius eqn. and activation energy of 0.22 eV.  Convert
% conductivity to e''.
E = eV*1.602e-19*6.02e23;
R = 8.3143;
physicalConstants;
Tref = 273.15-20;
er = er - j*acidity.*3.3.*exp(E/R*(-1./T+1/Tref))./(2*pi*freq)./e0;

% Specific gravity of ice is 0.917 g/cm^3
% Matsuoka et al, J of Applied Physics, vol 80, no 10, 15 Nov 1996 gives
% specific gravity of ice to be 0.9169 g/cm^3 at 270 K.
er = er .* (1 + 0.851*rho).^2/3.1697;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1
clf; clear all;
freq = linspace(12e9,18e9,20).';
w = 2*pi*freq;
physicalConstants;
sigma = 0;
er = ice(268.15,0.5,freq,2e-6,0.22);
atten = real(j*sqrt(-j*w.*u0.*(sigma+j*w.*er*e0)));
plot(freq,20*log10(exp(-atten*1)))

% Example 2
clf; clear all;
freq = linspace(50e6,500e6,20).';
w = 2*pi*freq;
physicalConstants;
sigma = 0;
er = ice(253.15,0.917,freq,2e-6);
atten = real(j*sqrt(-j*w.*u0.*(sigma+j*w.*er*e0)));
plot(freq,20*log10(exp(-atten*1000)))

% Example 3
% This plot agrees with pg 205 Fujita figure 9a results
% sigma = -2*pi*freq.*e0.*imag(er);
clf; clear all;
physicalConstants;
T = 192:263;
freq = 30e6; w = 2*pi*freq;
er = ice(T,0.917,freq,2e-6);
atten = 1e3*8.686*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
semilogy(T-273.15,atten,'r')
set(get(gcf,'CurrentAxes'),'XDir','reverse')

hold on;
freq = 60e6; w = 2*pi*freq;
er = ice(T,0.917,freq,2e-6);
atten = 1e3*8.686*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
semilogy(T-273.15,atten)
freq = 179e6; w = 2*pi*freq;
er = ice(T,0.917,freq,2e-6);
atten = 1e3*8.686*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
semilogy(T-273.15,atten)
freq = 300e6; w = 2*pi*freq;
er = ice(T,0.917,freq,2e-6);
atten = 1e3*8.686*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
semilogy(T-273.15,atten)
freq = 600e6; w = 2*pi*freq;
er = ice(T,0.917,freq,2e-6);
atten = 1e3*8.686*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
semilogy(T-273.15,atten)
freq = 1270e6; w = 2*pi*freq;
er = ice(T,0.917,freq,2e-6);
atten = 1e3*8.686*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
semilogy(T-273.15,atten,'g')
hold off;
xlabel('Temperature (C)');
ylabel('Loss (dB/km)');
axis([-80 -10 10^-1 10^1.5]);
title('Matches Fujita 2000 results with a small unexplained error');

% Example 4
clf; clear all;
physicalConstants;
T = 243:273;
freq = 150e6; w = 2*pi*freq;
er = ice(T,0.917,freq,0);
plot(T,real(er));
xlabel('Temperature (C)');
ylabel('Re\{\epsilon_r\}')

% Example 5
clf; clear all;
physicalConstants;
T = 253;
freq = logspace(7,9,1001); w = 2*pi*freq;
er = ice(T,0.917,freq,0);
semilogx(freq,real(er));
xlabel('Frequency (Hz)');
ylabel('Re\{\epsilon_r\}');


