function out = iemb(freq, sig, corr, theta, er, sp, xx)
%IEMB   Evaluates the Integral Equation Surface Backscattering Model
%
% -Inputs-
%       freq    Frequency [Hz]
%       sig     Surface RMS height [m]
%       corr    Surface correlation length [m]
%       theta   Nx1 incident angle vector [rad]
%       er      Relative dielectric constant of the surface
%       sp      Shape of surface correlation function: 1=Exp, 2=Gaussian
%       xx      Parameter used to specify the x-exponential corr. func.
% -Outputs
%       Nx6 matrix, first column Svv, second column Shh
%         cross-polarization, Svh, not supported
%         multi-scattering effects not included
% 
% Ported from Mathematica 5 implementation bundled with:
%
% Microwave Scattering and Emission Models for Users
% Adrian K. Fung, Kun-Shan Chen
% Artech House, 2010
% Chapter 3
%
% Ulrik Nielsen, May 2013
% ulrik@space.dtu.dk

% Recursive call to support vector input for 'theta'
if numel(theta)>1
    if ~exist('xx','var'), xx = [];end
    out = arrayfun(@(v)iemb(freq, sig, corr, v, er, sp, xx),theta(:), 'UniformOutput', false);
    out = cell2mat(out);
    return
end

fr = freq*1e-9;     
sig = sig*1e2;      
corr = corr*1e2;    

k = 2*pi*fr/30;
phi = 0;
sigma = sig;
L = corr;
mu_r = 1;
ks = k*sig;
kL = k*corr;
ks2 = ks*ks;

cs = cos(theta);
sf = 0.0;
cf = 1.0;
s = sin(theta);
cfs = -1.0;
csfs = -1;
sfs = 0.0;
css = cs;
ss = s;
s2 = s^2;
sq = sqrt(er - s2);

kx = k*s;
ky = 0.0;
ksx = -k*s;
ksy = 0.0;
kz = k*cs;
ksz = k*cs;
rt = sqrt(er - s^2);
Rvi = (er*cs - rt)/(er*cs + rt);
Rhi = (cs - rt)/(cs + rt);

wvnb = k*sqrt((ss*cfs - s)^2 + (ss*sfs)^2);
TS = 1;
error = 1e8;

while error > 1e-8
    TS = TS + 1;
    error = (ks^2*(cs + css)^2)^TS/factorial(TS);
end

n = (1:TS)';

if(sp == 1)
    w = L^2./n.^2.*(1 + ((wvnb*L)./n).^2).^(-1.5);
elseif(sp == 2)
    w = L^2./(2*n).*exp(-((wvnb*L)^2./(4*n)));
elseif(sp == 3)
elseif(sp == 4)
elseif(sp == 5)
elseif(sp == 7)
end

Rvo = (sqrt(er) - 1)/(sqrt(er) + 1); Rho = -Rvo;
Ft = 8*Rvo^2*s^2* ( (cs + sqrt(er - s^2))/(cs*sqrt(er - s^2)));

sumvar = (ks*cs).^(2*n)./factorial(n).*w;
St = (0.25*abs(Ft)^2 * sum(sumvar)) ...
    / sum(sumvar.*abs(Ft/2+2.^(n+1)*Rvo/cs*exp(-(ks*cs)^2)).^2);
  
St0 = 1/abs(1 + 8*Rvo/(cs*Ft))^2;
Tf = 1 - St/St0;

Rvt = Rvi + (Rvo - Rvi) * Tf;
Rht = Rhi + (Rho - Rhi) * Tf;

fvv = 2*Rvt/cs;
fvvo = 2*Rvo/cs;
fvvi = 2*Rvi/cs;
fhh = -2*Rht/cs;
fhho = -2*Rho/cs;
fhhi = -2*Rhi/cs;

Tv = 1 + Rvt; Tvm = 1 - Rvt;
Th = 1 + Rht; Thm = 1 - Rht;

Fvv = (s^2/cs - sq/er ) * Tv^2 ...
    - 2*s^2*Tv*Tvm *(1/cs + 1/sq) ...
    + (s^2/cs + (er*(1 + s^2))/sq)*Tvm^2;

Fhh = -(s^2/cs - sq )*Th^2 ...
    + 2*s^2*Th*Thm*( 1/cs + 1/sq) ...
    - ((s^2) /cs + (1 + s^2)/sq)*Thm^2;

Tvi = 1 + Rvi;
Tvmi = 1 - Rvi;
Thi = 1 + Rhi;
Thmi = 1 - Rhi;

Fvvi = (s^2/cs - sq/er )*Tvi^2 ...
    - 2*s^2*Tvi*Tvmi*(1/cs + 1/sq) ...
    + ((s^2) /cs + (er*(1 + s^2))/sq)*Tvmi^2;

Fhhi = -(s^2/cs - sq )*Thi^2 ...
    + 2*s^2*Thi*Thmi*(1/cs + 1/sq) ...
    - (s^2/cs + (1 + s^2)/sq)*Thmi^2;

Tvo = 1 + Rvo;
Tvmo = 1 - Rvo;
Tho = 1 + Rho;
Thmo = 1 - Rho;

Fvvo = (s^2/cs - sq/er )*Tvo^2 ...
    - 2*s^2*Tvo*Tvmo*( 1/cs + 1/sq) ...
    + (s^2/cs + (er*(1 + s^2))/sq)*Tvmo^2;

Fhho = -(s^2/cs - sq )*Tho^2 ...
    + 2*s^2*Tho*Thmo*(1/cs + 1/sq) ...
    - ((s^2) /cs + (1 + s^2)/sq)*Thmo^2;

ex = exp(-(ks*cs)^2);

Ivv = ks.^n.*((2*cs).^n*fvv*ex + (cs.^n*Fvv));
Ihh = ks.^n.*((2*cs).^n*fhh*ex + (cs.^n*Fhh));

Ivvo = ks.^n.*((2*cs).^n*fvvo*ex + (cs.^n*Fvvo));
Ihho = ks.^n.*((2*cs).^n*fhho*ex + (cs.^n*Fhho));

Ivvi = ks.^n.*((2*cs).^n*fvvi*ex + (cs.^n*Fvvi));
Ihhi = ks.^n.*((2*cs).^n*fhhi*ex + (cs.^n*Fhhi));

ct = cot(theta);
cts = cot(theta);
  
if sp == 1
    rss = sqrt(1)*sig/L;
elseif sp == 2
    rss = sqrt(2)*sig/L;
elseif sp == 3
    rss = sqrt(2*xx)*sig/L;
elseif sp == 4
    rss = sqrt(1)*sig/L;
elseif sp == 5
    rss = sqrt(2/(L*xx))*sig;
elseif sp == 6
    rss = sqrt(2/(L*xx))*sig;
elseif sp == 7
    rss = sqrt(2)*sig/xx;
end

rslp = rss;

ctorslp = ct/(sqrt(2)*rslp);
ctsorslp = cts/(sqrt(2)*rslp);
shadf = 1/2*(1/(sqrt(pi)*ctorslp) * exp(-ctorslp^2) - erfc(ctorslp));
shadfs = 1/2*(1/(sqrt(pi)*ctsorslp) * exp(-ctsorslp^2) - erfc(ctsorslp));
ShdwS = 1/(1 + shadf + shadfs);

ampterm = ShdwS * k^2/2 * exp(-2*sigma^2*kz^2);
sumvar = w./factorial(n);
sigmavv = ampterm * sum(abs(Ivv).^2.*sumvar);
sigmahh = ampterm * sum(abs(Ihh).^2.*sumvar);
osigmavv = ampterm * sum(abs(Ivvo).^2.*sumvar);
osigmahh = ampterm * sum(abs(Ihho).^2.*sumvar);
isigmavv = ampterm * sum(abs(Ivvi).^2.*sumvar);
isigmahh = ampterm * sum(abs(Ihhi).^2.*sumvar);

ssv = 10*log10(sigmavv);
ossv = 10*log10(osigmavv);
issv = 10*log10(isigmavv);
ssh = 10*log10(sigmahh);
ossh = 10*log10(osigmahh);
issh = 10*log10(isigmahh);

out = [ssv, ssh, ossv, ossh, issv, issh];
