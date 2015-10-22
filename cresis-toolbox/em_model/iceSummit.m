function er = iceSummit(rho)
% er = iceSummit(rho)
%
% rho = density of ice/snow (g/cm^3)
%
% rho is an N-dim matrix. If the
% parameter doesn't change, then it can be input as a scalar.
% Basically, you just have to satisfy the point-wise operators (e.g .*, .^, ./).
%
% e0*er = e0*e' (does not find the imaginary part)

% Real part from:
% Christian Matzler and Urs Wegmuller, "Dielectric properties of fresh-water ice
% at microwave frequencies," Journal of Physics. D: Applied Physics, vol. 20,
% pp. 1623-1630, 1987.
% Imag part is also from Matzler/Wegmuller, but we use the coefficients determined
% by Matsuoka:
% Takeshi Matsuoka, Shuji Fujita, and Shinji Mae, "Effect of temperature on
% dielectric properties of ice in the range 5-39 GHz," Journal of Applied
% Physics, vol. 80, no. 10, pp. 5884-5890, 1996.

er = 3.11;
er = 3.15;

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
er = real(er) .* (1 + 1.7*rho + 0.7*rho.^2)./(1 + 1.7*rho_ice + 0.7*rho_ice.^2);

return;

