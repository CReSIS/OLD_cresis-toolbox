function out = iiemb(freq, sig, corr, theta, er, sp, xx)
%IIEMB   Evaluates the Integral Equation Surface Backscattering Model
%
% -Inputs-
%       freq    Frequency [Hz]
%       sig     Surface RMS height [m]
%       corr    Surface correlation length [m]
%       theta   Incident angle [rad]
%       er      Relative dielectric constant of the surface
%       sp      Shape of surface correlation function: 1=Exp, 2=Gaussian
%       xx      Parameter used to specify the x-exponential corr. func.
%
% Ported from Mathematica 5 implementation bundled with:
%
% Microwave Scattering and Emission Models for Users
% Adrian K. Fung, Kun-Shan Chen
% Artech House, 2010
%
% Ulrik Nielsen, May 2013
% ulrik@space.dtu.dk

% Recursive call to support vector input for 'theta'
if numel(theta)>1
    if ~exist('xx','var'), xx = [];end    
    out = arrayfun(@(v)iiemb(freq, sig, corr, v, er, sp, xx),theta(:), 'UniformOutput', false);
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

cs = cos(theta);
sf = 0.0;
cf = 1.0;
s = sin(theta);
cfs = -1.0;
sfs = 0.0;
css = cs;
ss = s;
s2 = s^2;

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
elseif(sp == 6)
elseif(sp == 7)
end

Rv0 = (er-sqrt(er))/(er + sqrt(er)); 
Rh0 = -Rv0;
Ft = 8*Rv0^2*s^2* ( (cs + sqrt(er - s^2))/(cs*sqrt(er - s^2)));

sumvar = (ks*cs).^(2*n)./factorial(n).*w;
St = (0.25*abs(Ft)^2 * sum(sumvar)) ...
    / sum(sumvar.*abs(Ft/2+2.^(n+1)*Rv0/cs*exp(-(ks*cs)^2)).^2);
  
St0 = 1/abs(1 + 8*Rv0/(cs*Ft))^2;
Tf = 1 - St/St0;

Rtv = Rvi + (Rv0 - Rvi) * Tf;
Rth = Rhi + (Rh0 - Rhi) * Tf;

fvv = 2*Rtv/cs;
fhh = -2*Rth/cs;

Fvvdni = (2*k)/sqrt(er - s^2) * (-(1 + Rtv)^2 * (cs + s^2/(2*er)*(sqrt(er - s^2) - cs)) + (1 + Rtv)*(1 - Rtv)*s^2 * (sqrt(er - s^2) - cs) + er*(1 - Rtv)^2*cs);
Fvvups = Fvvdni;
Fhhdni = (2*k)/sqrt(er - s^2) * ( (1 + Rth)^2 * ( er*cs + s^2/2*(sqrt(er - s^2) - cs)) - (1 + Rth)*(1 - Rth)*s^2*(sqrt(er - s^2) - cs) - (1 - Rth)^2*cs);
Fhhups = Fhhdni;

Fvvupi = 2*k*s^2*((1 + Rtv)^2*(1 + 1/er) - (1 - Rtv)*(1 + Rtv)*(3 + cs/sqrt(er - s^2)) + (1 - Rtv)^2*(1 + ( cs*er)/sqrt(er - s^2)));
Fvvdns = 2*k*s^2*((1 + Rtv)^2*(1 + cs/(er*sqrt(er - s^2))) - (1 - Rtv)*(1 + Rtv)*(3 +  cs/sqrt(er - s^2)) + (1 - Rtv)^2*(1 + ( cs*er)/sqrt(er - s^2)));
Fhhupi = 2*k*s^2*(-(1 + Rth)^2*(2) + (1 - Rth)*(1 + Rth)*(3 + cs/sqrt(er - s^2)) - (1 - Rth)^2*(1 +  cs/sqrt(er - s^2)));
Fhhdns = 2*k*s^2*(-(1 + Rth)^2*(1 + cs/sqrt(er - s^2)) + (1 - Rth)*(1 + Rth)*(3 + cs/sqrt(er - s^2)) - (1 - Rth)^2*(1 + cs/sqrt(er - s^2)));

Ivv = (2*kz).^n*fvv + 1/4*(Fvvdni + Fvvups)*(2*kz).^(n - 1);
Ihh = (2*kz).^n*fhh + 1/4*(Fhhdni + Fhhups)*(2*kz).^(n - 1);

ct = cot(theta);
cts = cot(theta);

if sp == 1
    rss = sqrt(4)*sig/L;
elseif sp == 2
    rss = sqrt(2)*sig/L;
elseif sp == 3
    rss = sqrt(2*xx)*sig/L;
elseif sp == 4
    rss = sqrt(4)*sig/L;
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


ampterm = ShdwS * k^2/2 * exp(-4*sigma^2*kz^2);
sumvar = sigma.^(2*n).*w./factorial(n);
sumvv = abs(Ivv).^2.*sumvar;
sumhh = abs(Ihh).^2.*sumvar;

sigmavv = ampterm * (abs((2*kz*sigma)*fvv+sigma/4*(Fvvupi + Fvvdns + Fvvdni + Fvvups))^2*w(1) + sum(sumvv(2:end)));
sigmahh = ampterm * (abs((2*kz*sigma)*fhh+sigma/4*(Fhhupi + Fhhdns + Fhhdni + Fhhups))^2*w(1) + sum(sumhh(2:end)));

ssv = 10*log10(sigmavv);
ssh = 10*log10(sigmahh);

out = [ssv, ssh];