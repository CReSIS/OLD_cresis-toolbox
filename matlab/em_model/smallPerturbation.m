function scatter = smallPerturbation(er1,er2,ur1,ur2,corrLen,h,pdf,freq,thetaInc,phiInc,thetaScat,phiScat,etaInc,etaScat)
% scatter = smallPerturbation(er1,er2,ur1,ur2,corrLen,h,pdf,freq,thetaInc,phiInc,thetaScat,phiScat,etaInc,etaScat)
% er1,ur1: relative permittivity and permeability of medium 1 (radar located in this medium)
% er2,ur2: relative permittivity and permeability of medium 2
% corrLen: mean correlation length of roughness (m)
% h: RMS height of surface (m)
% pdf: 'gauss' or 'exp' for correlation coefficient model
% freq: frequency (Hz)
% thetaInc: elevation incidence angle (rad)
% phiScat: azimuth incidence angle (rad)
% thetaScat: elevation scattering angle (rad)
% phiScat: azimuth scattering angle (rad)
% etaInc: incidence polarization angle (0 is v, pi/2 is h)
% etaScat: scattering polarization angle (0 is v, pi/2 is h)
%
% scatter: normalized backscatter in dB/m^2
%
% Assumptions, pp. 704
% (a) k0*h < 1.0, i.e. roughness height is small
% Let zeta = member function of surface height
% (b) abs(derivative of zeta w.r.t. x) < 1.0
%     and abs(derivative of zeta w.r.t. y) < 1.0, i.e. surface slopes are small
% (c) Expectation{abs(derivative of zeta w.r.t. x).^2}
%     = Expectation{abs(derivative of zeta w.r.t. y).^2}, i.e. roughness is isotropic
%
% Example and recreation of plots from Ulaby, Moore, and Fung vol. 2 and
% Ruck et al. vol. 2 are at the bottom of the this file.

physicalConstants;
er = er2/er1;
ur = ur2/ur1;

% Bistatic backscatter (multi-polarization)
%  Note that hh is dual of vv and hv is dual of vh
alpha_hh = -( (ur-1).*(ur.*sin(thetaInc).*sin(thetaScat) ...
                      + cos(phiScat-phiInc).*sqrt(er.*ur-sin(thetaInc).^2) ...
                                    .*sqrt(er.*ur-sin(thetaScat).^2)) ...
             - ur.^2.*(er-1).*cos(phiScat-phiInc) ) ...
           ./ (ur.*cos(thetaInc) + sqrt(er.*ur-sin(thetaInc).^2)) ...
           ./ (ur.*cos(thetaScat) + sqrt(er.*ur-sin(thetaScat).^2));

alpha_vv = ( (er-1).*(er.*sin(thetaInc).*sin(thetaScat) ...
                      + cos(phiScat-phiInc).*sqrt(er.*ur-sin(thetaInc).^2) ...
                                    .*sqrt(er.*ur-sin(thetaScat).^2)) ...
             - er.^2.*(ur-1).*cos(phiScat-phiInc) ) ...
           ./ (er.*cos(thetaInc) + sqrt(er.*ur-sin(thetaInc).^2)) ...
           ./ (er.*cos(thetaScat) + sqrt(er.*ur-sin(thetaScat).^2));

alpha_vh = -sin(phiScat-phiInc).*(er.*(ur-1).*sqrt(er.*ur-sin(thetaInc).^2) ...
                          - ur.*(er-1).*sqrt(er.*ur-sin(thetaScat).^2)) ...
           ./ (ur.*cos(thetaInc) + sqrt(er.*ur-sin(thetaInc).^2)) ...
           ./ (er.*cos(thetaScat) + sqrt(er.*ur-sin(thetaScat).^2));

alpha_hv = -sin(phiScat-phiInc).*(ur.*(er-1).*sqrt(er.*ur-sin(thetaInc).^2) ...
                          - er.*(ur-1).*sqrt(er.*ur-sin(thetaScat).^2)) ...
           ./ (er.*cos(thetaInc) + sqrt(er.*ur-sin(thetaInc).^2)) ...
           ./ (ur.*cos(thetaScat) + sqrt(er.*ur-sin(thetaScat).^2));

% Roughness spectral density
k1 = 2*pi*freq.*sqrt(er1.*e0.*ur1.*u0);
xi_x = -sin(thetaInc).*cos(phiInc) - sin(thetaScat).*cos(phiScat);
xi_y = -sin(thetaInc).*sin(phiInc) - sin(thetaScat).*sin(phiScat);
% xi_z = -cos(thetaInc) - cos(phiScat); <-- Incorrect form in Ruck et al.'s book
xi_z = -cos(thetaInc) - cos(thetaScat);
if strcmp(pdf,'gauss')
  % Gaussian
  I = pi*corrLen.^2.*exp(-k1.^2.*corrLen.^2.*(xi_x.^2+xi_y.^2)./4);
else
  % Exponential
  I = 2*pi*corrLen.^2./(1 + k1.^2.*corrLen.^2.*(xi_x.^2+xi_y.^2)).^1.5;
end

alpha = alpha_vv.*cos(etaInc).*cos(etaScat) + alpha_vh*sin(etaInc)*cos(etaScat) ...
  + alpha_hv.*cos(etaInc).*sin(etaScat) + alpha_hh.*sin(etaInc).*sin(etaScat);
  
% abs(alpha) insures that complex er/ur work
scatter = 4/pi.*k1.^4.*h.^2.*cos(thetaInc).^2.*cos(thetaScat).^2.*abs(alpha).^2.*I;

return;

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% Examples
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

function test
% Slightly Rough Surface
%   Based on small perturbation theory
%
%   Ruck, Barrick, Stuart, and Krichbaum, "Radar Cross Section Handbook Volume 2," 1970, pp. 703-718.
%     - References the generic bistatic situation, but no proofs
%   Ulaby, Moore, and Fung, "Microwave Remote Sensing, Active and Passive,
%   Volume II, Radar Remote Sensing and Surface Scattering and Emission Theory,"
%   Artech House, 1986, pp. 949-966.
%     - Proves the generic bistatic situation, but only finished the proof for monostatic (backscatter)
%     - Result is identical to Ruck for backscatter
%
%   Ruck references:
%   Barrick, D. E. and W. H. Peake, "Scattering from surfaces with different roughness scales: analysis and interpretation," Batelle Memorial Institute, Columbus Laboratories, Columbus Ohio, Research Report BAT-197-A-10-3 (November 1967), AD 662 751.

er1 = 1;
er2 = [1+1/12 1+1/4 2 5 55*(1+0.55i)];
er2 = 20;
ur1 = 1;
ur2 = 1;
physicalConstants;
% Set freq so that k1 = 1;
freq = 1./real(2*pi*sqrt(er1.*e0.*ur1.*u0));
k = 2*pi*freq.*sqrt(e0.*er1.*u0.*ur1);
phiScat = pi;
etaInc = 0;
etaScat = 0;
corrLen = 2.5;
h = 0.025;

pdf = 'gauss';

thetaInc = [1:1:89]/180*pi;
for ind = 1:length(thetaInc)
  thetaScat = thetaInc(ind);
  freq = 120e6;
  scatterRuck1(ind,:) = smallPerturbation(er1,er2,ur1,ur2,corrLen,h,pdf,freq,thetaInc(ind),pi,thetaScat,phiScat,etaInc,etaScat);
  freq = 240e6;
  scatterRuck2(ind,:) = smallPerturbation(er1,er2,ur1,ur2,corrLen,h,pdf,freq,thetaInc(ind),pi,thetaScat,phiScat,etaInc,etaScat);
end

plot(10*log10(scatterRuck1));
hold on;
plot(10*log10(scatterRuck2));
hold off;
title('Frequency Dependence');
xlabel('Angle from vertical (deg)');
ylabel('(dB)');
axis([0 80 -35 5]);
return;


thetaInc = [1:1:89]/180*pi;
for ind = 1:length(thetaInc)
  thetaScat = thetaInc(ind);
  corrLen = 0.2;
  h = 1;
  scatterRuck(ind,:) = smallPerturbation(er1,er2,ur1,ur2,corrLen,h,pdf,freq,thetaInc(ind),pi,thetaScat,phiScat,etaInc,etaScat);
  corrLen = 1/k;
  h = 0.1/k;
  scatterUlaby1(ind,:) = smallPerturbation(er1,er2,ur1,ur2,corrLen,h,pdf,freq,thetaInc(ind),pi,thetaScat,phiScat,etaInc,etaScat);
  corrLen = 2/k;
  h = 0.2/k;
  scatterUlaby2(ind,:) = smallPerturbation(er1,er2,ur1,ur2,corrLen,h,pdf,freq,thetaInc(ind),pi,thetaScat,phiScat,etaInc,etaScat);
  corrLen = 3/k;
  h = 0.3/k;
  scatterUlaby3(ind,:) = smallPerturbation(er1,er2,ur1,ur2,corrLen,h,pdf,freq,thetaInc(ind),pi,thetaScat,phiScat,etaInc,etaScat);
  corrLen = 1/k;
  h = 0.1/k;
  scatterUlaby1_hh(ind,:) = smallPerturbation(er1,er2,ur1,ur2,corrLen,h,pdf,freq,thetaInc(ind),pi,thetaScat,phiScat,pi/2,pi/2);
  corrLen = 2/k;
  h = 0.2/k;
  scatterUlaby2_hh(ind,:) = smallPerturbation(er1,er2,ur1,ur2,corrLen,h,pdf,freq,thetaInc(ind),pi,thetaScat,phiScat,pi/2,pi/2);
  corrLen = 3/k;
  h = 0.3/k;
  scatterUlaby3_hh(ind,:) = smallPerturbation(er1,er2,ur1,ur2,corrLen,h,pdf,freq,thetaInc(ind),pi,thetaScat,phiScat,pi/2,pi/2);
end

plot(10*log10(scatterRuck));
title('Fig. 9-12 of Ruck');
xlabel('Angle from vertical (deg)');
ylabel('(dB)');
pause;

plot(10*log10(scatterUlaby1));
hold on;
plot(10*log10(scatterUlaby2));
plot(10*log10(scatterUlaby3));
hold off;
title('Fig. 12.6 of Ulaby, Moore, and Fung');
xlabel('Angle from vertical (deg)');
ylabel('(dB)');
axis([0 80 -35 5]);
pause;

plot(10*log10(scatterUlaby1_hh));
hold on;
plot(10*log10(scatterUlaby2_hh));
plot(10*log10(scatterUlaby3_hh));
hold off;
title('Fig. 12.6 of Ulaby, Moore, and Fung');
xlabel('Angle from vertical (deg)');
ylabel('(dB)');
axis([0 80 -50 5]);

return;
