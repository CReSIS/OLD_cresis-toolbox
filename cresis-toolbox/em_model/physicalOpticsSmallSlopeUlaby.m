function scatter = physicalOpticsSmallSlopeUlaby(er1,er2,ur1,ur2,corrLen,h,pdf,freq, ...
    thetaInc,phiInc,thetaScat,phiScat,etaInc,etaScat)
% scatter = physicalOpticsSmallSlopeUlaby(er1,er2,ur1,ur2,s,pdf,freq, ...
%   thetaInc,phiInc,thetaScat,phiScat,etaInc,etaScat)
%
% er1,ur1: relative permittivity and permeability of medium 1 (radar located in this medium)
% er2,ur2: relative permittivity and permeability of medium 2
% corrLen: mean correlation length of roughness (m)
% h: RMS height of surface (m)
% pdf: 'gauss' or 'exp' for correlation coefficient model
% freq: frequency (Hz) <-- there is no actual frequency dependence as it all cancels out
% thetaInc: elevation incidence angle (rad)
% phiInc:azimuth incidence angle (rad)
% thetaScat: elevation scattering angle (rad)
% phiScat: azimuth scattering angle (rad)
% etaInc: incidence polarization angle (0 is v, pi/2 is h)
% etaScat: scattering polarization angle (0 is v, pi/2 is h)
%
% scatter: normalized to unit area backscatter in dB/m^2
%   --> DOES NOT INCLUDE COHERENT TERM, even for specular.  Coherent term
%   is included on page 1009&938 AND IS NOT the usual coherent term found from
%   small perturbation theory... it is the square of that term.
%
% Assumptions
% (a) s < 0.25, i.e. slopes are small
% Let zeta = member function of surface height
% (b) Expectation{abs(derivative of zeta w.r.t. x).^2}
%     = Expectation{abs(derivative of zeta w.r.t. y).^2}, i.e. roughness is isotropic
% (c) corrLen << sqrt(A) where A is illuminated area (requires many
%     surface features in the area of illumination to get "averaged" effect
% (d) Multiple scattering and shadowing effects are neglected
% (e) Radius of curvature is large with respect to wavelength
% (f) k1*corrLen > 6 and corrLen^2 > 2.76*sigma*wavelength
%
% SOLUTION INVOLVES AN INFINITE SERIES... IT MUST BE TRUNCATED TO A FINITE SERIES
%   You should set N (termination point) to a value such that N+M produces the
%   same result (to the precision you are interested in) for all positive M.
%   i.e. all terms that are being truncated are insignificant.
%
% phiInc is poorly defined in Ulaby (it is defined as the normal "phiInc"
% definition + pi). This is corrected here... i.e. use the normal definition.
%
% Example and recreation of plots from Ulaby, Moore, and Fung vol. 2 are at
% the bottom of the this file

physicalConstants;
k1 = 2*pi*freq.*sqrt(e0.*er1.*u0.*ur1);
k2 = 2*pi*freq.*sqrt(e0.*er2.*u0.*ur2);
eta1 = sqrt(u0*ur1./(e0.*er1));
eta2 = sqrt(u0*ur2./(e0.*er2));

qx = k1*(sin(thetaScat)*cos(phiScat) + sin(thetaInc)*cos(phiInc));
qy = k1*(sin(thetaScat)*sin(phiScat) + sin(thetaInc)*sin(phiInc));
qz = k1*(cos(thetaScat) + cos(thetaInc));
q = sqrt(qx.^2 + qy.^2 + qz.^2);

thetaTran = asin(real(k1)./real(k2).*sin(thetaInc));

Rperp0 = (eta2.*cos(thetaInc) - eta1.*cos(thetaTran)) ...
  ./ (eta2.*cos(thetaInc) + eta1.*cos(thetaTran));
Rperp1 = -Rperp0.*(eta2.*sin(thetaInc) + eta1.*sin(thetaTran)) ...
  ./ (eta2.*cos(thetaInc) + eta1.*cos(thetaTran));

Rpara0 = (eta1.*cos(thetaInc) - eta2.*cos(thetaTran)) ...
  ./ (eta1.*cos(thetaInc) + eta2.*cos(thetaTran));
Rpara1 = -(eta1.*sin(thetaInc) - eta2.*sin(thetaTran) ...
    - Rpara0.*(eta1.*sin(thetaInc) + eta2.*sin(thetaTran))) ...
  ./ (eta1.*cos(thetaInc) + eta2.*cos(thetaTran));

a0_hh = Rperp0.*(cos(thetaInc) + cos(thetaScat)).*cos(phiScat - phiInc);
a_hh = Rperp0.*(sin(thetaScat) + sin(thetaInc).*cos(phiScat-phiInc)) ...
  + Rperp1.*(cos(thetaScat) + cos(thetaInc)).*cos(phiScat-phiInc);
a1_hh = a_hh*cos(phiInc);
a2_hh = a_hh*sin(phiInc);

a0_vh = +Rperp0.*(1+cos(thetaInc).*cos(thetaScat)).*sin(phiScat-phiInc);
a_vh = -(-Rperp0.*sin(thetaInc).*cos(thetaScat) ...
    - Rperp1.*(1+cos(thetaInc).*cos(thetaScat))).*sin(phiScat-phiInc);
a1_vh = a_vh*cos(phiInc);
a2_vh = a_vh*sin(phiInc);

a0_hv = +Rpara0.*(1+cos(thetaInc).*cos(thetaScat)).*sin(phiScat-phiInc);
a_hv = -(-Rpara0.*sin(thetaInc).*cos(thetaScat) ...
    - Rpara1.*(1+cos(thetaInc).*cos(thetaScat))).*sin(phiScat-phiInc);
a1_hv = a_hv*cos(phiInc);
a2_hv = a_hv*sin(phiInc);

% Ulaby, Moore, and Fung have conflicting equations for the a0 and a variables
%   This particular combination seems to reproduce the backscattering results
%   that they plotted, but does not exactly match either of the sets of equations
%   that they show:
%     pg 941 and pg 1004
a0_vv = Rpara0.*(cos(thetaInc) + cos(thetaScat)).*cos(phiScat - phiInc);
a_vv = -Rpara0.*(sin(thetaScat) + sin(thetaInc).*cos(phiScat-phiInc)) ...
  - Rpara1.*(cos(thetaScat) + cos(thetaInc)).*cos(phiScat-phiInc);
a1_vv = a_vv*cos(phiInc);
a2_vv = a_vv*sin(phiInc);

if strcmp(pdf,'gauss')
  % Gaussian
  N = 140;
  n = 1:N;
  sigmaNonCoh_hh = (abs(a0_hh)*k1*corrLen/2).^2 .* exp(-qz.^2.*h.^2) ...
   .* sum([ (qz.^2.*h.^2).^n./(gamma(1+n).*n) .* exp(-(qx.^2+qy.^2).*corrLen.^2./(4.*n)) ]);

  sigmaSlope_hh = -(k1.*h.*corrLen).^2.*(qz/2).*exp(-qz.^2.*h.^2) ...
    .* real(a0_hh.*(qx.*conj(a1_hh)+qy.*conj(a2_hh))) ...
    .* sum([ (qz.^2.*h.^2).^(n-1)./(gamma(1+n).*n).*exp(-(qx.^2+qy.^2).*corrLen.^2./(4*n)) ]);

  sigmaNonCoh_vh = (abs(a0_vh)*k1*corrLen/2).^2 .* exp(-qz.^2.*h.^2) ...
   .* sum([ (qz.^2.*h.^2).^n./(gamma(1+n).*n) .* exp(-(qx.^2+qy.^2).*corrLen.^2./(4.*n)) ]);

  sigmaSlope_vh = -(k1.*h.*corrLen).^2.*(qz/2).*exp(-qz.^2.*h.^2) ...
    .* real(a0_vh.*(qx.*conj(a1_vh)+qy.*conj(a2_vh))) ...
    .* sum([ (qz.^2.*h.^2).^(n-1)./(gamma(1+n).*n).*exp(-(qx.^2+qy.^2).*corrLen.^2./(4*n)) ]);

  sigmaNonCoh_hv = (abs(a0_hv)*k1*corrLen/2).^2 .* exp(-qz.^2.*h.^2) ...
   .* sum([ (qz.^2.*h.^2).^n./(gamma(1+n).*n) .* exp(-(qx.^2+qy.^2).*corrLen.^2./(4.*n)) ]);

  sigmaSlope_hv = -(k1.*h.*corrLen).^2.*(qz/2).*exp(-qz.^2.*h.^2) ...
    .* real(a0_hv.*(qx.*conj(a1_hv)+qy.*conj(a2_hv))) ...
    .* sum([ (qz.^2.*h.^2).^(n-1)./(gamma(1+n).*n).*exp(-(qx.^2+qy.^2).*corrLen.^2./(4*n)) ]);

  sigmaNonCoh_vv = (abs(a0_vv)*k1*corrLen/2).^2 .* exp(-qz.^2.*h.^2) ...
   .* sum([ (qz.^2.*h.^2).^n./(gamma(1+n).*n) .* exp(-(qx.^2+qy.^2).*corrLen.^2./(4.*n)) ]);

  sigmaSlope_vv = -(k1.*h.*corrLen).^2.*(qz/2).*exp(-qz.^2.*h.^2) ...
    .* real(a0_vv.*(qx.*conj(a1_vv)+qy.*conj(a2_vv))) ...
    .* sum([ (qz.^2.*h.^2).^(n-1)./(gamma(1+n).*n).*exp(-(qx.^2+qy.^2).*corrLen.^2./(4*n)) ]);
end

sigma_hh = sigmaNonCoh_hh + abs(sigmaSlope_hh);
sigma_hv = sigmaNonCoh_hv + sigmaSlope_hv;
sigma_vh = sigmaNonCoh_vh + sigmaSlope_vh;
sigma_vv = sigmaNonCoh_vv + abs(sigmaSlope_vv);

scatter = sigma_vv.*cos(etaInc).*cos(etaScat) + sigma_vh*sin(etaInc)*cos(etaScat) ...
  + sigma_hv.*cos(etaInc).*sin(etaScat) + sigma_hh.*sin(etaInc).*sin(etaScat);

return;

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% Examples
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

function test
% Small slopes, large radius of curvature surface, "undulating surface"
%   Based on geometric optics (ray optics, physical optics, "kirchoff method")
%   and scalar approximation assumptions
%   Ulaby, Moore, and Fung, "Microwave Remote Sensing, Active and Passive,
%   Volume II, Radar Remote Sensing and Surface Scattering and Emission Theory,"
%   Artech House, 1986, pp. 925-929 and 936-943.
%
%   Ulaby, Moore, and Fung reference Beckmann and Spizzichino... but derived in the text.
% 
clear all; format compact;

er1 = 1;
er2 = 1.6;
ur1 = 1;
ur2 = 1;
physicalConstants;
% Set freq so that k1 = 1;
freq = 1./real(2*pi*sqrt(er1.*e0.*ur1.*u0));
phiScat = pi;
phiInc = pi;
s = 0.1;
h = 1.5;
corrLen = 1.25*h/s;

  freq = 120e6;
  freq = 480e6;
  h = 0.25;
  corrLen = 2.0;

pdf = 'gauss';

k1 = 2*pi*freq.*sqrt(e0.*er1.*u0.*ur1);
lambda = 2*pi./real(k1);
if (k1*corrLen > 6)
  fprintf('k1*corrLen = %.2f > 6\n',k1*corrLen);
else
  fprintf('k1*corrLen = %.2f ~> 6 (VIOLATION!)\n',k1*corrLen);
end
if (corrLen^2 > 2.76*h*lambda)
  fprintf('corrLen^2 = %.2f > 2.76*h*lambda = %.2f\n',corrLen^2,2.76*h*lambda);
else
  fprintf('corrLen^2 = %.2f ~> 2.76*h*lambda = %.2f (VIOLATION!)\n',corrLen^2,2.76*h*lambda);
end
fprintf('Surface slope %.2f < 0.25\n',sqrt(4*h.^2./corrLen.^2));

thetaInc = [0:1:80]/180*pi;
for ind = 1:length(thetaInc)
  thetaScat = thetaInc(ind);
  freq = 120e6;
  scatter1(ind,:) = physicalOpticsSmallSlopeUlaby(er1,er2,ur1,ur2,corrLen,h, ...
      pdf,freq,thetaInc(ind),phiInc,thetaScat,phiScat,pi/2,pi/2);
  freq = 240e6;
  scatter2(ind,:) = physicalOpticsSmallSlopeUlaby(er1,er2,ur1,ur2,corrLen,h, ...
      pdf,freq,thetaInc(ind),phiInc,thetaScat,phiScat,pi/2,pi/2);
  freq = 480e6;
  scatter3(ind,:) = physicalOpticsSmallSlopeUlaby(er1,er2,ur1,ur2,corrLen,h, ...
      pdf,freq,thetaInc(ind),phiInc,thetaScat,phiScat,pi/2,pi/2);
  freq = 960e6;
  scatter4(ind,:) = physicalOpticsSmallSlopeUlaby(er1,er2,ur1,ur2,corrLen,h, ...
      pdf,freq,thetaInc(ind),phiInc,thetaScat,phiScat,pi/2,pi/2);
end

plot(thetaInc*180/pi,10*log10(scatter1),'k-');
hold on;
plot(thetaInc*180/pi,10*log10(scatter2),'r-');
plot(thetaInc*180/pi,10*log10(scatter3),'g-');
plot(thetaInc*180/pi,10*log10(scatter4),'b-');
hold off;
title('Frequency Dependence');
xlabel('Angle of Incidence (deg)');
ylabel('Backscattering coefficient (dB)');
axis([0 35 -50 0]);
return;

thetaInc = [0:1:80]/180*pi;
for ind = 1:length(thetaInc)
  thetaScat = thetaInc(ind);
  s = 0.1;
  h = 0.5/k1;
  corrLen = 1.35*h/s;
  scatter0_5hh(ind,:) = physicalOpticsSmallSlopeUlaby(er1,er2,ur1,ur2,corrLen,h, ...
      pdf,freq,thetaInc(ind),phiInc,thetaScat,phiScat,pi/2,pi/2);
  scatter0_5vv(ind,:) = physicalOpticsSmallSlopeUlaby(er1,er2,ur1,ur2,corrLen,h, ...
      pdf,freq,thetaInc(ind),phiInc,thetaScat,phiScat,0,0);
  s = 0.1;
  h = 0.6/k1;
  corrLen = 1.35*h/s;
  scatter0_6hh(ind,:) = physicalOpticsSmallSlopeUlaby(er1,er2,ur1,ur2,corrLen,h, ...
      pdf,freq,thetaInc(ind),phiInc,thetaScat,phiScat,pi/2,pi/2);
  scatter0_6vv(ind,:) = physicalOpticsSmallSlopeUlaby(er1,er2,ur1,ur2,corrLen,h, ...
      pdf,freq,thetaInc(ind),phiInc,thetaScat,phiScat,0,0);
  s = 0.1;
  h = 1.5/k1;
  corrLen = 1.35*h/s;
  scatter1_5hh(ind,:) = physicalOpticsSmallSlopeUlaby(er1,er2,ur1,ur2,corrLen,h, ...
      pdf,freq,thetaInc(ind),phiInc,thetaScat,phiScat,pi/2,pi/2);
  scatter1_5vv(ind,:) = physicalOpticsSmallSlopeUlaby(er1,er2,ur1,ur2,corrLen,h, ...
      pdf,freq,thetaInc(ind),phiInc,thetaScat,phiScat,0,0);
  scatter(ind,:) = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,s, ...
      pdf,freq,thetaInc(ind),phiInc,thetaScat,phiScat,0,0);
end

plot(thetaInc*180/pi,10*log10(scatter0_5hh),'ko');
hold on;
plot(thetaInc*180/pi,10*log10(scatter0_5vv),'k-x');
plot(thetaInc*180/pi,10*log10(scatter0_6hh),'ro');
plot(thetaInc*180/pi,10*log10(scatter0_6vv),'r-x');
plot(thetaInc*180/pi,10*log10(scatter1_5hh),'bo');
plot(thetaInc*180/pi,10*log10(scatter1_5vv),'b-x');
plot(thetaInc*180/pi,10*log10(scatter),'g-');
hold off;
title('Fig. 12.4 Ulaby, Moore, and Fung');
xlabel('Angle of Incidence (deg)');
ylabel('Backscattering coefficient (dB)');
axis([0 35 -50 0]);
return;
