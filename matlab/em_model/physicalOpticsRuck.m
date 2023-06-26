function scatter = physicalOpticsRuck(er1,er2,ur1,ur2,h,s,pdf,thetaInc,phiInc,thetaScat,phiScat,etaInc,etaScat)
% scatter = physicalOpticsRuck(er1,er2,ur1,ur2,h,s,pdf,thetaInc,phiInc,thetaScat,phiScat,etaInc,etaScat)
%
% er1,ur1: relative permittivity and permeability of medium 1 (radar located in this medium)
% er2,ur2: relative permittivity and permeability of medium 2
% h: RMS height of interface (including only wavelength components of the surface
%    that are near the radar's operating wavelengths) <-- use for 2-scale model
% s: RMS slope of interface
% pdf: 'gauss' or 'exp' for correlation coefficient model
% thetaInc: elevation incidence angle (rad)
% phiInc:azimuth incidence angle (rad)
% thetaScat: elevation scattering angle (rad)
% phiScat: azimuth scattering angle (rad)
% etaInc: incidence polarization angle (0 is v, pi/2 is h)
% etaScat: scattering polarization angle (0 is v, pi/2 is h)
%
% scatter: normalized to unit area backscatter in dB/m^2
%
% Assumptions
% (a) k0*h > 5.0, i.e. roughness height is large
% Let zeta = member function of surface height
% (b) Expectation{abs(derivative of zeta w.r.t. x).^2}
%     = Expectation{abs(derivative of zeta w.r.t. y).^2}, i.e. roughness is isotropic
% (c) corrLen << sqrt(A) where A is illuminated area (requires many
%     surface features in the area of illumination to get "averaged" effect
% (d) Multiple scattering, shadowing and local diffraction effects are neglected
% (e) Radius of curvature is large with respect to wavelength
%
% Example and recreation of plots from Ruck et al. vol. 2 are at
% the bottom of the this file
physicalConstants;
er = er2/er1;
ur = ur2/ur1;

% Handle limiting cases separately (only beta terms need to be handled)
if (thetaScat == thetaInc & phiInc == phiScat)
  inc = 0;
  Rpara = (er.*cos(inc) - sqrt(er.*ur - sin(inc).^2)) ...
          ./ (er.*cos(inc) + sqrt(er.*ur - sin(inc).^2)) .* exp(-(qz.*h).^2);
  Rperp = (ur.*cos(inc) - sqrt(er.*ur - sin(inc).^2)) ...
          ./ (ur.*cos(inc) + sqrt(er.*ur - sin(inc).^2)) .* exp(-(qz.*h).^2);
  % beta_hh = sec(thetaInc)*Rpara; <-- incorrect from in Ruck et al.'s book
  beta_hh = -sec(thetaInc)*Rpara;
  beta_vv = -sec(thetaInc)*Rpara;
  beta_vh = 0;
  beta_hv = 0;
else
  % Bistatic backscatter (multi-polarization)
  %  Note that hh is dual of vv and hv is dual of vh
  inc = acos(1/sqrt(2)*sqrt(1 + sin(thetaInc).*sin(thetaScat).*cos(phiScat-phiInc) + cos(thetaInc).*cos(thetaScat)));
  a1 = 1 - sin(thetaInc).*sin(thetaScat).*cos(phiScat-phiInc) - cos(thetaInc).*cos(thetaScat);
  a2 = cos(thetaInc).*sin(thetaScat) - sin(thetaInc).*cos(thetaScat).*cos(phiScat-phiInc);
  a3 = sin(thetaInc).*cos(thetaScat) - cos(thetaInc).*sin(thetaScat).*cos(phiScat-phiInc);
  a4 = cos(thetaInc) + cos(thetaScat);

  Rpara = (er.*cos(inc) - sqrt(er.*ur - sin(inc).^2)) ...
          ./ (er.*cos(inc) + sqrt(er.*ur - sin(inc).^2)) .* exp(-(qz.*h).^2);
  Rperp = (ur.*cos(inc) - sqrt(er.*ur - sin(inc).^2)) ...
          ./ (ur.*cos(inc) + sqrt(er.*ur - sin(inc).^2)) .* exp(-(qz.*h).^2);

  beta_vv = (a2.*a3.*Rpara + sin(thetaInc).*sin(thetaScat).*sin(phiScat-phiInc).^2.*Rperp) ./ (a1.*a4);
  beta_vh = -sin(phiScat-phiInc).*(sin(thetaScat).*a2.*Rperp - sin(thetaInc).*a3.*Rpara) ./ (a1.*a4);
  beta_hv = -sin(phiScat-phiInc).*(sin(thetaScat).*a2.*Rpara - sin(thetaInc).*a3.*Rperp) ./ (a1.*a4);
  beta_hh = (-a2.*a3.*Rperp - sin(thetaInc).*sin(thetaScat).*sin(phiScat-phiInc).^2.*Rpara) ./ (a1.*a4);

end

beta = beta_vv.*cos(etaInc).*cos(etaScat) + beta_vh*sin(etaInc)*cos(etaScat) ...
  + beta_hv.*cos(etaInc).*sin(etaScat) + beta_hh.*sin(etaInc).*sin(etaScat);

% Roughness spectral density
xi_x = +sin(thetaInc).*cos(phiInc) + sin(thetaScat).*cos(phiScat);
xi_y = +sin(thetaInc).*sin(phiInc) + sin(thetaScat).*sin(phiScat);
% xi_z = -cos(thetaInc) - cos(phiScat); <-- Incorrect form in Ruck et al.'s book
xi_z = +cos(thetaInc) + cos(thetaScat);
if strcmp(pdf,'gauss')
  % Gaussian
  J = 4./(s.^2.*xi_z.^2).*exp(-1./s.^2.*(xi_x.^2+xi_y.^2)./xi_z.^2);
else
  % Exponential
  J = 12./(s.^2.*xi_z.^2).*exp(-sqrt(6)./s.*sqrt((xi_x.^2+xi_y.^2)./xi_z.^2));
end
  
% abs(beta) insures that complex er/ur work
scatter = abs(beta).^2.*J;

return;

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% Examples
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

function test
% Very rough surface, "undulating surface"
%   Based on geometric optics (ray optics, physical optics, "kirchoff method")
%   and stationary phase assumptions
%   Ruck, Barrick, Stuart, and Krichbaum, "Radar Cross Section Handbook Volume 2," 1970, pp. 719-727.
%
%   Ruck references:
%  --> In the process of finding the proofs for these equations

er1 = 1;
er2 = 2;
ur1 = 1;
ur2 = 1;
physicalConstants;
% Set freq so that k1 = 1;
freq = 1./real(2*pi*sqrt(er1.*e0.*ur1.*u0));
s = tan([2.5 5 10 15 20 25 30]/180*pi)
etaInc = 0;
etaScat = 0;
phiInc = pi;
phiScat = pi;

thetaInc = [0:1:89]/180*pi;
for ind = 1:length(thetaInc)
  thetaScat = thetaInc(ind);
  scatterGauss(ind,:) = physicalOpticsRuck(er1,er2,ur1,ur2,s,'gauss', ...
      thetaInc(ind),phiInc,thetaScat,phiScat,etaInc,etaScat);
  scatterExp(ind,:) = physicalOpticsRuck(er1,er2,ur1,ur2,s,'exp', ...
      thetaInc(ind),phiInc,thetaScat,phiScat,etaInc,etaScat);
end

eta1 = sqrt(ur1./er1);
eta2 = sqrt(ur2./er2);
R0 = (eta2 - eta1)./(eta2 + eta1);

plot(thetaInc*180/pi,10*log10(scatterGauss./abs(R0).^2));
title('Fig. 9-16 of Ruck');
xlabel('Angle from vertical (deg)');
ylabel('(dB)');
axis([0 90 -50 30]);
pause;

plot(thetaInc*180/pi,10*log10(scatterExp./abs(R0).^2));
title('Fig. 9-16 of Ruck');
xlabel('Angle from vertical (deg)');
ylabel('(dB)');
axis([0 90 -50 30]);
pause;

% ------------------------------------------------------------------------------
% Ulaby Figure 12.3
% ------------------------------------------------------------------------------
er2 = 81;
s = 0.4*sqrt(2);
pdf = 'gauss';

thetaInc = 4.5/180*pi;
thetaScat = thetaInc;

phiScat = [0:1:180]/180*pi;
for ind = 1:length(phiScat)
  scatter_vv(ind,:) = physicalOpticsRuck(er1,er2,ur1,ur2,s,pdf,...
    thetaInc,pi,thetaScat,phiScat(ind),0,0);
  scatter_hh(ind,:) = physicalOpticsRuck(er1,er2,ur1,ur2,s,pdf,...
    thetaInc,pi,thetaScat,phiScat(ind),pi/2,pi/2);
  scatter_hv(ind,:) = physicalOpticsRuck(er1,er2,ur1,ur2,s,pdf,...
    thetaInc,pi,thetaScat,phiScat(ind),0,pi/2);
end

plot(phiScat*180/pi,10*log10(scatter_vv),'k-');
hold on;
plot(phiScat*180/pi,10*log10(scatter_hh),'r-.');
plot(phiScat*180/pi,10*log10(scatter_hv),'g:');
hold off;
title('Fig. 12.3a Ulaby, Moore, and Fung');
xlabel('Azimuth Angle From Specular (deg)');
ylabel('Scattering coefficient (dB)');
axis([0 180 -40 10]);
pause;

thetaInc = 22.5/180*pi;
thetaScat = thetaInc;

phiScat = [0:1:180]/180*pi;
for ind = 1:length(phiScat)
  scatter_vv(ind,:) = physicalOpticsRuck(er1,er2,ur1,ur2,s,pdf,...
    thetaInc,pi,thetaScat,phiScat(ind),0,0);
  scatter_hh(ind,:) = physicalOpticsRuck(er1,er2,ur1,ur2,s,pdf,...
    thetaInc,pi,thetaScat,phiScat(ind),pi/2,pi/2);
  scatter_hv(ind,:) = physicalOpticsRuck(er1,er2,ur1,ur2,s,pdf,...
    thetaInc,pi,thetaScat,phiScat(ind),0,pi/2);
end

plot(phiScat*180/pi,10*log10(scatter_vv),'k-');
hold on;
plot(phiScat*180/pi,10*log10(scatter_hh),'k-.');
plot(phiScat*180/pi,10*log10(scatter_hv),'k:');
hold off;
title('Fig. 12.3b Ulaby, Moore, and Fung');
xlabel('Azimuth Angle From Specular (deg)');
ylabel('Scattering coefficient (dB)');
axis([0 180 -40 10]);
pause;

thetaInc = 49.5/180*pi;
thetaScat = thetaInc;

phiScat = [0:1:180]/180*pi;
for ind = 1:length(phiScat)
  scatter_vv(ind,:) = physicalOpticsRuck(er1,er2,ur1,ur2,s,pdf,...
    thetaInc,pi,thetaScat,phiScat(ind),0,0);
  scatter_hh(ind,:) = physicalOpticsRuck(er1,er2,ur1,ur2,s,pdf,...
    thetaInc,pi,thetaScat,phiScat(ind),pi/2,pi/2);
  scatter_hv(ind,:) = physicalOpticsRuck(er1,er2,ur1,ur2,s,pdf,...
    thetaInc,pi,thetaScat,phiScat(ind),0,pi/2);
end

semilogy(phiScat*180/pi,scatter_vv,'k-');
hold on;
semilogy(phiScat*180/pi,scatter_hh,'k-.');
semilogy(phiScat*180/pi,scatter_hv,'k:');
hold off;
title('Fig. 12.3c Ulaby, Moore, and Fung');
xlabel('Azimuth Angle From Specular (deg)');
ylabel('Scattering coefficient (dB)');
axis([0 180 1e-4 10]);

return;

