function Gf = refractionGain(depth,er,er0,inc,z)
% Gf = refractionGain(depth,er,er0,inc,z)
% 
% depth = depth (m)
% er = dielectric profile corresponding to depth
% er0 = scalar which represent the medium the antenna is in
%   Primaryly two cases:
%     Antenna is air (airborne or antenna mounted above ground), er0 = 1
%       - depth,er should NOT include the air layer... start at firn
%     Antenna on snow surface er0 = er(1)
% inc = incidence angles (rad)
% z = depths to find refraction Gain for (m)
%
% depth and er are vectors of the same length and are used
% to determine the exponential density profile
%
% inc and z are N-dim matrices of the same size
% size(Gf) = size(inc) = size(z)
% Gf = refractionGain(depth,er,inc,z)
%
% Reference:
% West, James C. and Kenneth R. Demarest, "The radiation characteristics
% of an arbitrary antenna positioned on a polar ice sheet," Geophysics,
% vol. 52, no. 12, Dec 1987, pp. 1689-1696.
%   --> Typos out the wazoo.  See proofs/corrections from John, plus
%   derivation of incidence angle = 0 case.

format compact;
if (ispc)
  path(path,'P:\prism\radar\radarSimulator\');
else
  path(path,'/projects/prism/radar/radarSimulator/');
end
physicalConstants;

% Reformat depth and er to be row-vectors
dim = min(find(size(depth) ~= 1));
if isempty(dim), dim = 1; end;
depth = permute(depth,[dim 1:dim-1 dim+1:length(size(dim))]).';
dim = min(find(size(er) ~= 1));
if isempty(dim), dim = 1; end;
er = permute(er,[dim 1:dim-1 dim+1:length(size(dim))]).';

% Use density - real(er) relation used in paper (eqn 4)
den = (sqrt(real(er)) - 1)/0.854;

% Only use lower density data to produce curve fit
lastGoodIndex = min(find(den>0.9));
denLow = den(1:lastGoodIndex);
depthLow = depth(1:lastGoodIndex);

% Fit an exponential profile to the density curve.  This means finding
% V and R in the following equation:
%  denLow(z) = P - V*exp(R*z)
% where P = 0.917 is the density of pure ice.
%
% To do this, we turn the equation into a linear equation:
%
% ln(P - denLow(z)) = ln(V) + R*z
% 
% where we estimate a1 = ln(V) and R.  
% [a1 R] * A ~= data
% [a1 R ] ~= data*pinv(A)
% 
% where
% A = [ 1    1    1 . . .    1;
%      z(1) z(2) z(3) . . . z(n)];

P = 0.917;
data = log(P - denLow);
res = data*pinv([ones(1,length(depthLow)); depthLow]);
a1 = res(1);
R = res(2);
V = exp(a1);

% DEBUG: Uncomment to compare original density profile and fit
% plot(depth(1:end-1),den);
% hold on;
% plot(depth,P - V*exp(R*depth),'r');
% hold off;
% legend('Original Density','Exponential Fit');

% Example from West paper (ref: Byrd Ice Core Density Measurements, Clough PhD)
% P = 0.92;
% R = -0.033;
% V = 0.520;

gamma0 = inc;
n = 1.0 + 0.854*(P - V*exp(R*z));
% n = sqrt(1.0 + 1.7*(P - V*exp(R*z)) + 0.7*(P - V*exp(R*z)).^2);

n0 = sqrt(real(er0));
zeta = n0*sin(gamma0);

% a = (1 + 0.854*P).^2 - c.^2; % Eqn is wrong in paper
a = (1 + 0.854*P).^2 - zeta.^2;
b = -2*0.854*(0.854*P+1);
d = 0.854^2;
theta = V.*exp(R.*z);
X = a + b.*theta + d.*theta.^2;
Xv = a + b.*V + d.*V.^2;
gammaPrime = asin(n0./n.*sin(gamma0));

if zeta == 0 | gamma0 == 0
  % Handle normal incidence separately (requires taking limits since eqn is 0/0 form)
  rPrime = 0;
  D = z;
  rPrimeDIVsinGamma0 = - n0./(R.*sqrt(a)) .* log(V.*(2*sqrt(a.*X) + b.*theta + 2*a)./theta./(2*sqrt(a.*Xv) + b.*V + 2*a));
  Gf = D.^2 ./ (rPrimeDIVsinGamma0.^2.*cos(gamma0).*cos(gammaPrime));
else
  rPrime = - zeta./(R.*sqrt(a)) .* log(V.*(2*sqrt(a.*X) + b.*theta + 2*a)./theta./(2*sqrt(a.*Xv) + b.*V + 2*a));

  q = 4*a.*d - b.^2;
  drPrime_dGamma0 = rPrime.*cot(gamma0) + zeta.^2.*n0.*cos(gamma0)./R ...
   .*(rPrime.*R./(a.*zeta) - 2./(a.*q).*((b.*d.*theta - 2*a.*d+b.^2)./sqrt(X) - (b.*d.*V-2*a.*d+b.^2)./sqrt(Xv)));
   
  D = sqrt(rPrime.^2 + z.^2);

  Gf = D.^2 .* sin(gamma0) ./ (rPrime.*drPrime_dGamma0.*cos(gammaPrime));
end

% Handle limiting case when z == 0
Gf(find(z == 0)) = (n(1)/n0).^2;

% Asymptotic solutions (two forms)
% Gf = n.^2.*cos(gammaPrime)./(n0.*sqrt(n0.^2 - n.^2.*sin(gammaPrime).^2));
% Gf = n.^2.*cos(gammaPrime)./(n0.^2.*cos(gamma0));

return;

% ------------------------------------------------------------------------------
% Examples
% ------------------------------------------------------------------------------

clear all; format compact; format long;
[depth,er] = gisp2Perm(150e6,10001);
Gf = refractionGain(depth,er,er(1),0,depth(1:end-1));
plot(depth(1:end-1),10*log10(Gf),'k-');
xlabel('Depth (m)');
ylabel('Refraction Gain (dB)');
axis([depth(1) depth(end-1) 0 3]);

inc = linspace(0,pi/6,101);
Gf = refractionGain(depth,er,er(1),inc,3000);
plot(inc*180/pi,10*log10(Gf),'k-');
xlabel('Incidence Angle (deg)');
ylabel('Refraction Gain (dB)');

z = [50 100 200 300 400 600 1000];
Gf = refractionGain(depth,er,er(1),0,z);
plot(z,10*log10(Gf));


[depth,er] = gisp2Perm(150e6,1001);
Gfair = sqrt(real(er(1))).^2.*cos(0)./(1.^2.*cos(0));
Gfice = refractionGain(depth,er,er(1),0,depth(1:end-1));
Gf1 = Gfair*Gfice;

Gf2 = refractionGain(depth,er,1,0,depth(1:end-1));

plot(Gf1);
hold on;
plot(Gf2,'r:');
hold off;

