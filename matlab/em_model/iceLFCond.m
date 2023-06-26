function er = ice_jezek(T,rho,freq,conductivity,eV,Tref)
% er = ice(T,rho,freq,conductivity,eV,Tref)
%
% T = temperature in Kelvin
% rho = density of ice/snow (g/cm^3)
% freq = frequency in Hz
% conductivity = total high-frequency limit conductivity (e.g. measured by DEP)
%    This conductivity includes ice and impurities.
% eV = activation energy in electron Volts (typical: 0.22 eV)
% Tref = reference temperature in Kelvin that the conductivity was
%        measured at (typical: -15 C or -20 C)
%
% T, rho, and freq are N-dim matrices of the same size.  Except for T, if the
% parameter doesn't change, then it can be input as a scalar.
% Basically, you just have to satisfy the point-wise operators (e.g .*, .^, ./).
%
% e0*er = e0*(e' - j*e'')

% Real part from:
% Christian Matzler and Urs Wegmuller, "Dielectric properties of fresh-water ice
% at microwave frequencies," Journal of Physics. D: Applied Physics, vol. 20,
% pp. 1623-1630, 1987.
% Imag part is also from Matzler/Wegmuller, but we use the coefficients determined
% by Matsuoka:
% Takeshi Matsuoka, Shuji Fujita, and Shinji Mae, "Effect of temperature on
% dielectric properties of ice in the range 5-39 GHz," Journal of Applied
% Physics, vol. 80, no. 10, pp. 5884-5890, 1996.
er = (3.1884+0.00091*(T-273.15));

% Fujita et al, Physics of Ice Core Records, Hokkaido Press, 2000, pp. 185-212.
% Fujita suggests a frequency dependents of er in the microwave region which
% John Paden has tried to duplicate here... it is very small.  Also John Paden
% did not measure any dispersion from 110-500 MHz.
% er = (3.1884+0.00091*(T-273.15)) + (log10(freq)-7.18)*-0.012 ...
%    -j*(A./(freq/1e9) + B.*(freq/1e9).^C);

% Wolff, Physics of Ice Core Records, Hokkaido Press, 2000, pp. 155-171.
%
% There are two parts to the conductivity.  One part is from the ice
% and has an activation energy of 0.58 eV (although Wolff, page 165
% suggests that 0.5 eV be used for Greenland and Antarctic ice).  Wolff
% and we are assuming that at -15 C (the Tref), the ice conductivity is
% 9 uS.
%
% The other part is from impurities which has an activation energy
% given by the dataset author (the most common value is around 0.22 eV
% as suggested on page 165 of Wolff).  The Arrenhius equation used to
% scale the result is shown explicitly in (eqn 3):
% Hugh Corr, John C. Moore, and Keith E. Nicholls, "Radar absorption due to
% impurities in Antarctic ice," Geophysical Research Letters, vol. 20, no. 11,
% pp. 1071-1074, June 7, 1993.
% Note that in Corr's eqn (mieu_ssCl * [ssCl] = sigma_infinity at the ref temp)
E = eV*1.602e-19*6.02e23;
E2 = 0.57*1.602e-19*6.02e23;
R = 8.3143;
physicalConstants;
impurity_conductivity = conductivity - 9e-6;
impurity_conductivity(find(impurity_conductivity<0)) = 0;
er = er - j*impurity_conductivity.*exp(E/R*(-1./T+1/Tref))./(2*pi*freq)./e0;
er = er - j*9e-6.*exp(E2/R*(-1./T+1/Tref))./(2*pi*freq)./e0;

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
er = real(er) .* (1 + 1.7*rho + 0.7*rho.^2)./(1 + 1.7*rho_ice + 0.7*rho_ice.^2) ...
      + j*imag(er) .* (0.52*rho + 0.62*rho.^2) ./ (0.52*rho_ice + 0.62*rho_ice.^2);

return;

