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
% genReflTranFromPerm([1 lambda/8 lambda/5 1]', ones(size(er)), er, freq)
%
% Author: Joe Lilek, John Paden
%
% See also: ?

%function [refl,tran,tx_angle] = genReflTranFromPerm(thick,ur,er,freq,inc)
physical_constants
f = freq.';
l = thick/c;
n = sqrt(er).';
theta = inc*180/pi;

[Gamma,Z] = multidiel2(n,l,f,theta,'tm');
refl = Gamma;

return;

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



% multidiel2.m - reflection response of lossy isotropic multilayer dielectric structures
%
%          na | n1 | n2 | ... | nM | nb
% left medium | l1 | l2 | ... | lM | right medium 
%   interface 1    2    3     M   M+1
%
% Usage: [Gamma,Z] = multidiel2(n,l,f,theta,pol)
%        [Gamma,Z] = multidiel2(n,l,f,theta)       (equivalent to pol='te')
%        [Gamma,Z] = multidiel2(n,l,f)             (equivalent to theta=0)
%
% f     = N-dimensional vector of frequencies in units of f0, f = [f(1), f(2),..., f(N)]
% n     = Nx(M+2) matrix of complex refractive indices,  
%         i-th row represents the refractive indices at the i-th frequency f(i), 
%         that is, n(i,:) =  [na(i),n1(i),n2(i),...,nM(i),nb(i)]
% l     = M-dimensional vector of physical lengths of layers in units of la0, l = [l(1),...,l(M)]
% theta = incidence angle from left medium (in degrees)
% pol   = 'tm' or 'te', for parallel/perpendicular polarizations
%
% Gamma = reflection response (at interface 1) evaluated at the N frequencies
% Z     = input impedance at interface-1 in units of eta_a (left medium)
%
% notes: M is the number of layers (must be >=0)
%        it assumes isotropic layers
%
%        f is in units of some f0, i.e. f/f0
%        physical (not optical) layer thicknesses are in units of la0=c0/f0, i.e., l/la0

% S. J. Orfanidis - 2003 - www.ece.rutgers.edu/~orfanidi/ewa

% function [refl,tran,tx_angle] = genReflTranFromPerm(thick,ur,er,freq,inc)
% 
% n = sqrt(er);


function [Gamma,Z] = multidiel2(n,l,f,theta,pol)

if nargin==0, help multidiel2; return; end
if nargin<=4, pol='te'; end
if nargin==3, theta=0; end

M = size(n,2)-2;                                    % number of layers
theta = theta * pi/180;

for i=1:length(f),                                  % frequency loop
   nf = n(i,:);                                     % complex refractive indices at i-th frequency   

   costh = sqrt(1 - (nf(1)*sin(theta)./nf).^2);     % nf(1) is na

   if pol=='te' | pol=='TE',
      nT = nf .* costh;                             % transverse refractive indices
   else
      nT = nf ./ costh;                             % TM case, bombs at 90 deg for left medium
   end

   r = -diff(nT) ./ (diff(nT) + 2*nT(1:M+1));       % transverse reflection coefficients

   Gamma1 = r(M+1);                                 % initialize Gamma at right-most interface

   for m = M:-1:1,
      delta = 2*pi*f(i)*l(m) * sqrt(nf(m)^2 - nf(1)^2*sin(theta)^2);    % phase thickness in m-th layer
      z = exp(-2*j*delta);                                        
      Gamma1 = (r(m) + Gamma1.*z) ./ (1 + r(m)*Gamma1.*z);           % backward recursion
   end

   Gamma(i) = Gamma1;                               % reflection response at i-th frequency    
   Z(i) = (1 + Gamma1) ./ (1 - Gamma1);
end  

end

% 
% 
% % multidiel.m - reflection response of isotropic or birefringent multilayer structure
% %
% %          na | n1 | n2 | ... | nM | nb
% % left medium | L1 | L2 | ... | LM | right medium 
% %   interface 1    2    3     M   M+1
% %
% % Usage: [Gamma,Z] = multidiel(n,L,lambda,theta,pol)
% %        [Gamma,Z] = multidiel(n,L,lambda,theta)       (equivalent to pol='te')
% %        [Gamma,Z] = multidiel(n,L,lambda)             (equivalent to theta=0)
% %
% % n      = isotropic 1x(M+2), uniaxial 2x(M+2), or biaxial 3x(M+2), matrix of refractive indices
% % L      = vector of optical lengths of layers, in units of lambda_0
% % lambda = vector of free-space wavelengths at which to evaluate the reflection response
% % theta  = incidence angle from left medium (in degrees)
% % pol    = for 'tm' or 'te', parallel or perpendicular, p or s, polarizations
% %
% % Gamma = reflection response at interface-1 into left medium evaluated at lambda 
% % Z     = transverse wave impedance at interface-1 in units of eta_a (left medium)
% %
% % notes: M is the number of layers (M >= 0)
% %
% %        n = [na, n1, n2, ..., nM, nb]        = 1x(M+2) row vector of isotropic indices
% %
% %            [ na1  n11  n12  ...  n1M  nb1 ]   3x(M+2) matrix of birefringent indices, 
% %        n = [ na2  n21  n22  ...  n2M  nb2 ] = if 2x(M+2), it is extended to 3x(M+2)
% %            [ na3  n31  n32  ...  n3M  nb3 ]   by repeating the top row
% %
% %        optical lengths are in units of a reference free-space wavelength lambda_0:
% %        for i=1,2,...,M,  L(i) = n(1,i) * l(i), for TM, 
% %                          L(i) = n(2,i) * l(i), for TE,
% %        TM and TE L(i) are the same in isotropic case. If M=0, use L=[].
% %
% %        lambda is also in units of lambda_0, that is, lambda/lambda_0 = f_0/f
% %
% %        reflectance = |Gamma|^2, input impedance = Z = (1+Gamma)./(1-Gamma)
% %
% %        delta(i) = 2*pi*[n(1,i) * l(i) * sqrt(1 - (Na*sin(theta))^2 ./ n(3,i).^2))]/lambda, for TM
% %        delta(i) = 2*pi*[n(2,i) * l(i) * sqrt(1 - (Na*sin(theta))^2 ./ n(2,i).^2))]/lambda, for TE
% %
% %        if n(3,i)=n(3,i+1)=Na, then will get NaN's at theta=90 because of 0/0, (see also FRESNEL)
% 
% % S. J. Orfanidis - 2000 - www.ece.rutgers.edu/~orfanidi/ewa
% 
% function [Gamma,Z] = multidiel(n,L,lambda,theta,pol)
% 
% if nargin==0, help multidiel; return; end
% if nargin<=4, pol='te'; end
% if nargin==3, theta=0; end
% 
% if size(n,2)==1, n = n'; end                            % in case n is entered as column 
% 
% K = size(n,1);                                          % birefringence dimension
% M = size(n,2)-2;                                        % number of layers
% 
% if K==1, n = [n; n; n]; end                             % isotropic case
% if K==2, n = [n(1,:); n]; end                           % uniaxial case
% 
% if M==0, L = []; end                                    % single interface, no slabs
% 
% theta = theta * pi/180;
% 
% if pol=='te',
%     Nsin2 = (n(2,1)*sin(theta))^2;                      % (Na*sin(tha))^2              
%     c = sqrt(1 - Nsin2 ./ n(2,:).^2);                   % coefficient ci, or cos(th(i)) in isotropic case
%     nT = n(2,:) .* c;                                   % transverse refractive indices
%     r = n2r(nT);                                        % r(i) = (nT(i-1)-nT(i)) / (nT(i-1)+nT(i))
% else
%     Nsin2 = (n(1,1)*n(3,1)*sin(theta))^2 / (n(3,1)^2*cos(theta)^2 + n(1,1)^2*sin(theta)^2);
%     c = sqrt(1 - Nsin2 ./ n(3,:).^2);
%     nTinv = c ./ n(1,:);                                % nTinv(i) = 1/nT(i) to avoid NaNs
%     r = -n2r(nTinv);                                    % minus sign because n2r(n) = -n2r(1./n)
% end
% 
% if M>0,
%     L = L .* c(2:M+1);                                  % polarization-dependent optical lengths
% end
% 
% Gamma = r(M+1) * ones(1,length(lambda));                % initialize Gamma at right-most interface
% 
% for i = M:-1:1,                                         % forward layer recursion 
%     delta = 2*pi*L(i)./lambda;                          % phase thickness in i-th layer
%     z = exp(-2*j*delta);                          
%     Gamma = (r(i) + Gamma.*z) ./ (1 + r(i)*Gamma.*z);
% end
% 
% Z = (1 + Gamma) ./ (1 - Gamma);
% 
