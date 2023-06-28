function modRefl = coherentRefl(origRefl,k0,h,thetaInc)
% modRefl = coherentRefl(origRefl,k0,h,thetaInc)
%
% origRefl: original power reflection coefficient
% k0: wavenumber in incident medium
% h: RMS height of surface
% thetaInc: elevation incidence angle (rad)
%
% modRefl: modified reflection coeff for rough surface
%
% Coherent reflection coefficients for rough planar surface of infinite extent
%
% Reference:
% Ruck, Barrick, Stuart, and Krichbaum, "Radar Cross Section Handbook Volume 2," 1970, pp. 700-701.
% Ulaby, Moore, and Fung, "Microwave Remote Sensing," pp. 951.
%
% Original References
% Ament, W. W., "Toward a theory of reflection by a rough surface," Proceedings IRE, vol. 41:142 (1953)
%
% Davies, H., "The reflection of electromagnetic waves from a rough surface," Proceedings of IEE (GB), V. 101, Part IV:209 1954.

modRefl = origRefl.*exp(-2*(k0.*h.*cos(thetaInc)).^2);

return;

