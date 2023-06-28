function scatter = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,h,s,pdf,thetaInc,phiInc,thetaScat,phiScat,etaInc,etaScat)
% scatter = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,h,s,pdf,thetaInc,phiInc,thetaScat,phiScat,etaInc,etaScat)
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
% (a) (qz*sigma).^2 > 10.0, i.e. roughness height is large
% Let zeta = member function of surface height
% (b) Expectation{abs(derivative of zeta w.r.t. x).^2}
%     = Expectation{abs(derivative of zeta w.r.t. y).^2}, i.e. roughness is isotropic
% (c) corrLen << sqrt(A) where A is illuminated area (requires many
%     surface features in the area of illumination to get "averaged" effect
% (d) Multiple scattering, shadowing, and local diffraction effects are neglected
% (e) Radius of curvature is large with respect to wavelength
% (f) k1*corrLen > 6 and corrLen^2 > 2.76*sigma*wavelength
%
% Example and recreation of plots from Ulaby, Moore, and Fung vol. 2 are at
% the bottom of the this file

freq = 1;
physical_constants;
k1 = 2*pi*freq.*sqrt(e0.*er1.*u0.*ur1);
k2 = 2*pi*freq.*sqrt(e0.*er2.*u0.*ur2);
eta1 = sqrt(u0*ur1./(e0.*er1));
eta2 = sqrt(u0*ur2./(e0.*er2));

qx = k1*(sin(thetaScat)*cos(phiScat) + sin(thetaInc)*cos(phiInc));
qy = k1*(sin(thetaScat)*sin(phiScat) + sin(thetaInc)*sin(phiInc));
qz = k1*(cos(thetaScat) + cos(thetaInc));
q = sqrt(qx.^2 + qy.^2 + qz.^2);

% Handle limiting cases separately (only beta terms need to be handled)
if (thetaScat == thetaInc & mod(phiInc,2*pi) == mod(phiScat,2*pi))
  Rpara = (k2-k1)./(k2+k1) .* exp(-(qz.*h).^2);
  % Deal with polarization!
  Uhh = Rpara;
  Uvv = Rpara;
  Uhv = 0;
  Uvh = 0;
  U = Uvv.*cos(etaInc).*cos(etaScat) + Uvh*sin(etaInc)*cos(etaScat) ...
    + Uhv.*cos(etaInc).*sin(etaScat) + Uhh.*sin(etaInc).*sin(etaScat);
  scatter = abs(U).^2.*exp(-(tan(thetaInc).^2./(2*s.^2)))./(2*s.^2.*cos(thetaInc).^4);
else
  % Bistatic backscatter (multi-polarization)
  %  Note that hh is dual of vv and hv is dual of vh

  inc = acos(q*abs(qz)/(2*k1*qz));
  Rpara = (k2.*cos(inc) - k1.*cos(asin(k1/k2*sin(inc)))) ...
    ./ (k2.*cos(inc) + k1*cos(asin(k1/k2*sin(inc)))) .* exp(-(qz.*h).^2);
  Rperp = (k1.*cos(inc) - k2.*cos(asin(k1/k2*sin(inc)))) ...
    ./ (k1.*cos(inc) + k2*cos(asin(k1/k2*sin(inc)))) .* exp(-(qz.*h).^2);
  % Rpara = (eta2.*cos(inc) - eta1.*cos(asin(eta1/eta2*sin(inc)))) ./ (eta2.*cos(inc) + eta1*cos(asin(eta1/eta2*sin(inc))))
  % Rperp = -(eta1.*cos(inc) - eta2.*cos(asin(eta1/eta2*sin(inc)))) ./ (eta1.*cos(inc) + eta2*cos(asin(eta1/eta2*sin(inc))))

  vs_ni = +sin(thetaInc).*cos(thetaScat).*cos(phiScat-phiInc) - cos(thetaInc).*sin(thetaScat);
  v_ns = -cos(thetaInc).*sin(thetaScat).*cos(phiScat-phiInc) + sin(thetaInc).*cos(thetaScat);
  hs_ni = +sin(thetaInc).*sin(phiScat-phiInc);
  h_ns = -sin(thetaScat).*sin(phiScat-phiInc);

  Ml = q.*abs(qz)./((hs_ni.^2 + vs_ni.^2).*k1.*qz);
  Uhh = Ml .* (Rpara.*hs_ni.*h_ns +  Rperp.*vs_ni.*v_ns);
  Uvh = Ml .* (Rpara.*vs_ni.*h_ns -  Rperp.*hs_ni.*v_ns);
  Uhv = Ml .* (Rpara.*hs_ni.*v_ns -  Rperp.*vs_ni.*h_ns);
  Uvv = Ml .* (Rpara.*vs_ni.*v_ns +  Rperp.*hs_ni.*h_ns);
  U = Uvv.*cos(etaInc).*cos(etaScat) + Uvh*sin(etaInc)*cos(etaScat) ...
    + Uhv.*cos(etaInc).*sin(etaScat) + Uhh.*sin(etaInc).*sin(etaScat);

  if strcmp(pdf,'gauss')
    % Gaussian
    scatter = (k1.*q.*abs(U)).^2./(2*qz.^4.*s.^2) .* exp(-(qx.^2+qy.^2)./(2.*qz.^2.*s.^2));
  end

end
  
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
%   Ulaby, Moore, and Fung, "Microwave Remote Sensing, Active and Passive,
%   Volume II, Radar Remote Sensing and Surface Scattering and Emission Theory,"
%   Artech House, 1986, pp. 925-936.
%
%   Ulaby, Moore, and Fung reference nobody... derived in the text.
% 
clear; close all; format compact;

% ------------------------------------------------------------------------------
% Ulaby Figure 12.2 and Ruck Figure 9-16
% ------------------------------------------------------------------------------
er1 = 1;
% er2 = 1.4; % Test Ruck
er2 = 1.6;
ur1 = 1;
ur2 = 1;
physical_constants;
% Set freq so that k1 = 1;
freq = 1./real(2*pi*sqrt(er1.*e0.*ur1.*u0));
% s = tan([2.5 5 10 15 20 25 30]/180*pi)/sqrt(2) % Test Ruck Fig. 9-16 case <-- requires 1/sqrt(2) to match
s = [0.1 0.2 0.3 0.4]; % Ulaby Fig. 12.2
etaInc = 0;
etaScat = 0;

pdf = 'gauss';

thetaInc = [0:1:89]/180*pi;
for ind = 1:length(thetaInc)
  thetaScat = thetaInc(ind);
  scatter(ind,:) = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,0,s,pdf,thetaInc(ind),0,thetaScat,0,etaInc,etaScat);
end

eta1 = sqrt(ur1./er1);
eta2 = sqrt(ur2./er2);
R0 = (eta2 - eta1)./(eta2 + eta1);

% plot(thetaInc*180/pi,10*log10(scatter./abs(R0).^2)); % RUCK
plot(thetaInc*180/pi,10*log10(scatter)); % Ulaby
title('Fig. 12.2 Ulaby, Moore, and Fung');
xlabel('Angle of Incidence (deg)');
ylabel('Backscattering coefficient (dB)');
% axis([0 90 -50 30]); % RUCK
axis([0 60 -50 0]); % Ulaby
pause;

% ------------------------------------------------------------------------------
% Ulaby Figure 12.3
% ------------------------------------------------------------------------------
er2 = 81;
s = 0.4;

thetaInc = 4.5/180*pi;
thetaScat = thetaInc;

phiScat = [0:1:180]/180*pi;
for ind = 1:length(phiScat)
  scatter_vv(ind,:) = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,0,s,pdf, ...
    thetaInc,pi,thetaScat,phiScat(ind),0,0);
  scatter_hh(ind,:) = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,0,s,pdf, ...
    thetaInc,pi,thetaScat,phiScat(ind),pi/2,pi/2);
  scatter_hv(ind,:) = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,0,s,pdf, ...
    thetaInc,pi,thetaScat,phiScat(ind),0,pi/2);
end

plot(phiScat*180/pi,10*log10(scatter_vv),'k-');
hold on;
plot(phiScat*180/pi,10*log10(scatter_hh),'k-.');
plot(phiScat*180/pi,10*log10(scatter_hv),'k:');
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
  scatter_vv(ind,:) = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,0,s,pdf, ...
    thetaInc,pi,thetaScat,phiScat(ind),0,0);
  scatter_hh(ind,:) = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,0,s,pdf, ...
    thetaInc,pi,thetaScat,phiScat(ind),pi/2,pi/2);
  scatter_hv(ind,:) = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,0,s,pdf, ...
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
  scatter_vv(ind,:) = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,0,s,pdf, ...
    thetaInc,pi,thetaScat,phiScat(ind),0,0);
  scatter_hh(ind,:) = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,0,s,pdf, ...
    thetaInc,pi,thetaScat,phiScat(ind),pi/2,pi/2);
  scatter_hv(ind,:) = physicalOpticsLargeUlaby(er1,er2,ur1,ur2,0,s,pdf, ...
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
