function er = ice_westphal(T,rho,freq)
% er = ice_westphal(T,rho,freq)
%
% T = temperature in Kelvin
% rho = density of ice/snow (g/cm^3)
% freq = frequency in Hz
%
% T, rho, and freq are N-dim matrices of the same size.  Except for T, if the
% parameter doesn't change, then it can be input as a scalar.
% Basically, you just have to satisfy the point-wise operators (e.g .*, .^, ./).
%
% e0*er = e0*(e' - j*e'')

% Use Tuto Tunnel data collected by Westphal
load_westphal;
real_part = interp2(repmat(tuto.freq,[size(tuto.T,1) 1]), ...
  repmat(tuto.T,[1 size(tuto.freq,2)]), ...
  tuto.e_real,freq,T-273.15,'linear');
imag_part = interp2(repmat(tuto.freq,[size(tuto.T,1) 1]), ...
  repmat(tuto.T,[1 size(tuto.freq,2)]), ...
  tuto.e_imag.*repmat(tuto.freq,[size(tuto.T,1) 1]),freq,T-273.15,'linear');
er = real_part -j*imag_part./freq.';

% Specific gravity of ice is 0.917 g/cm^3
% Matsuoka et al, J of Applied Physics, vol 80, no 10, 15 Nov 1996 gives
% specific gravity of ice to be 0.9169 g/cm^3 at 270 K.
% Real Part:
% Robin Q., S. Evans, and J. Bailey, "Interpretation of radio echo sounding
% in polar ice sheets," Philosophical Transactions of the Royal Society of
% London. Series A, Mathematical and Physical Sciences, vol. 265, no. 1166,
% pp. 437-505, Dec. 18, 1969.
%     er' = (1 + 0.851*rho).^2 ==> Only gives 0.85 in the paper
% Real Part:
% Austin Kovacs, Anthony J. Gow, Rexford M. Morey, "The in-situ dielectric
% constant of polar firn revisited," Cold Regions Science and Technology,
% vol. 23, pp. 245-256, 1995.
%     er' = (1 + 0.845*rho).^2 ==> Gives nearly identical results to Tiuri
% Real and Imaginary Parts:
% Martti E. Tiuri, Ari H. Sihvola, Ebbe F. Nyfors, Martti T. Hallikaiken, "The
% complex dielectric constant of snow at microwave frequencies," IEEE Journal of
% Oceanic Engineering, vol. OE-9, no. 5, pp. 377-382, December 1984.
%     er' = 1 + 1.7*rho + 0.7*rho.^2
%     er_d''/er_ice'' = 0.52*rho + 0.62*rho.^2
rho_ice = 0.917;
if (size(er) ~= size(rho))
  er = real(er.') .* (1 + 1.7*rho + 0.7*rho.^2)./(1 + 1.7*rho_ice + 0.7*rho_ice.^2) ...
    + j*imag(er.') .* (0.52*rho + 0.62*rho.^2) ./ (0.52*rho_ice + 0.62*rho_ice.^2);
else
  er = real(er) .* (1 + 1.7*rho + 0.7*rho.^2)./(1 + 1.7*rho_ice + 0.7*rho_ice.^2) ...
    + j*imag(er) .* (0.52*rho + 0.62*rho.^2) ./ (0.52*rho_ice + 0.62*rho_ice.^2);
end

return;

clf; clear all;
freq = linspace(150e6,500e6,20).';
w = 2*pi*freq;
physicalConstants;
sigma = 0;
er = ice_westphal(263.15,0.917,freq).';
atten = real(j*sqrt(-j*w.*u0.*(sigma+j*w.*er*e0)));
plot(freq,20*log10(exp(-atten*1)))

clear all;
freq = 500e6;
w = 2*pi*freq;
physicalConstants;
sigma = 0;
T = 243.15:5:268.15;
er = ice_westphal(T,0.917,freq).';
atten = real(j*sqrt(-j*w.*u0.*(sigma+j*w.*er*e0)));
plot(T,20*log10(exp(-atten*1)),'r')

