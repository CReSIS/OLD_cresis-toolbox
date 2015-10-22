function [refl,tran,tx_angle] = genReflTranFromPerm(thick,ur,er,freq,inc)
% [refl,tran,tx_angle] = genReflTranFromPerm(thick,ur,er,freq,inc)
%
% Generates the reflection coefficient, transmission coefficient and
% transmitted angle for arbitrary incidence angle plane waves incident
% on layered media.
% WARNING: CURRENT IMPLEMENTATION ONLY SUPPORTS NORMAL INCIDENCE AND REFL!
%
% M = length of freq
% N = length of er
%
% refl: 1 by M vector containing complex reflection coefficients
%   from the layered media at each frequency
%
% thick: N-2 by 1 vector with thickness of each layer
% ur: N by M vector containing relative permeabiliity of free space
%   (relative to u0)
% er: N by M length matrix containing relative permittivity of free space
%   (relative to e0)
%   e.g. e0*er = e0*(e' - j*e'')
% freq: 1 by M length vector frequency (Hz)
%   You have to provide an er profile for each frequency such that:
%   er(:,n) corresponds to freq(n)
% inc = vector of incidence angles (rad)
% 
% For each freq,inc pair, you get the following outputs:
% refl = reflection coefficient (Er/Ei)
% tran = transmission coefficient (Et/Ei)
% tx_angle= transmission angle in er(n) media (rad)
%
% IMPLEMENTATION NOT DONE YET:
% pg. 473-484 in Radar Cross Section Handbook, Vol. 2 (blue book)
% Depth Index
%                           er(1)
%     1       ----------------------------------
%                           er(2)
%     2       ----------------------------------
%                           er(3)
%     3       ----------------------------------
%     .
%     .
%     .
%    n-1      ----------------------------------
%                           er(n-1)
%     n       ----------------------------------
%                           er(n)
%
% The reflection coefficient and incident angle is defined on interface 1 in
% medium er(1).
% The transmission coefficient and transmission angle is defined on interface n
% in medium er(n).
%
% CURRENT IMPLEMENTATION:
% Demarest, Engineering Electromagnetics, Prentice Hall, New Jersey,
% 1998, pp. 480-490.
%
% Depth Index
%                           er(1,:),ur(1,:)
%             ----------------------------------
%     thick(1)              er(2,:),ur(2,:)
%             ----------------------------------
%     thick(2)              er(3,:),ur(3,:)
%             ----------------------------------
%     .
%     .
%     .
%             ----------------------------------
%     thick(n-2)            er(n-1,:),ur(n-1,:)
%             ----------------------------------
%                           er(n,:),ur(n,:)
%
% Example:
% physical_constants;
% freq = 100e6;
% lambda = c/freq;
% er = [1 6.25 2.25 1].';
% refl_jl2([1 lambda/8 lambda/5 1]', ones(size(er)), er, freq)
%
% Author: Joe Lilek, John Paden
%
% See also: ?

M = length(freq); % columns
N = length(er); % rows

% Load u0 and e0
physical_constants

% omega: 1 by M vector of angular frequencies (radians/sec)
omega = 2*pi*freq;

% gamma: propagation constant (1/m)
gamma = 1j.*repmat(omega, [N 1]).*repmat(sqrt(u0*ur.*er*e0),[1 M]);

% eta: intrinsic impedance (Ohms)
eta = sqrt(u0*ur./(er*e0));
eta = repmat(eta, [1 M]);

% Temporary variables
TAN = tanh(gamma(1:N,:) .* repmat(thick, [1 M]));

% Determine the effective intrinsic impedance for all the layers by
% combining layers two at a time starting from the bottom.
for n = size(eta,1) : -1 : 3
  % Effective impedance from equation 12.159, pg 487
  eta(n-1,:) = eta(n-1,:).*((eta(n,:)+eta(n-1,:).*TAN(n-2,:)) ./ ((eta(n-1,:)+eta(n,:).*TAN(n-2,:))));
  if 0
    % Debug print
    fprintf('%.0f: effective impedance: %.2f %+.2fi\n', n-1, real(eta(n-1,1)), imag(eta(n-1,1)));
  end
end

refl = (eta(2,:) - eta(1,:)) ./ (eta(2,:) + eta(1,:));

end

% See example 12-10 in Engineering Electromagnetics (Demarest), pg. 489-90




