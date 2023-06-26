function out = i2emb(freq, sig1, corr1, sig2, corr2, theta, er)
%I2EMB   Evaluates the Integral Equation Two-Scale Surface Backscattering Model
%
% Gaussian shaped surface correlation function is assumed.
%
% -Inputs-
%       freq    Frequency [Hz]
%       sig1    Large scale surface RMS height [m]
%       corr1   Large scale surface correlation length [m]
%       sig2    Small scale surface RMS height [m]
%       corr2   Small scale surface correlation length [m]
%       theta   Incident angle [rad]
%       er      Relative dielectric constant of the surface
%
% Ported from Mathematica 5 implementation bundled with:
%
% Microwave Scattering and Emission Models for Users
% Adrian K. Fung, Kun-Shan Chen
% Artech House, 2010
%
% Ulrik Nielsen, April 2013
% ulrik@space.dtu.dk

% Recursive call to support vector input for 'theta'
if numel(theta)>1
    out = arrayfun(@(v)i2emb(freq, sig1, corr1, sig2, corr2, v, er),theta(:), 'UniformOutput', false);
    out = cell2mat(out);
    return
end

fr = freq*1e-9;     
sig1 = sig1*1e2;    
sig2 = sig2*1e2;    
corr1 = corr1*1e2;  
corr2 = corr2*1e2;  

k = 2*pi*fr/30;
phi = 0;
sigma12 = sig1^2;
sigma22 = sig2^2;
sigma = sqrt(sigma12+sigma22);
L1 = corr1;
L2 = corr2;
mu_r = 1;
ks = k*sigma;
kL1 = k*corr1;
kL2 = k*corr2;

cs = cos(theta);
sf = 0.0;
cf = 1.0;
s = sin(theta);
cfs = -1.0;
sfs = 0.0;
css = cs;
ss = s;
s2 = s^2;
sq = sqrt(er-s2);

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

F = @(n1)quad(@(x)(((sigma12*exp(-x.^2)+sigma22.*exp(-(x*L1/L2).^2))/sigma^2).^n1.*cos(2*k*s*L1*x)),0,10);
w = 2*L1*arrayfun(F,n);

Rv0 = (sqrt(er)-1)/(sqrt(er)+1); 
Rh0 = -Rv0;
Ft = 8*Rv0^2*s^2* ( (cs + sqrt(er - s^2))/(cs*sqrt(er - s^2)));

sumvar = (ks*cs).^(2*n)./factorial(n).*w;
St = (0.25*abs(Ft)^2 * sum(sumvar)) ...
    / sum(sumvar.*abs(Ft/2+2.^(n+1)*Rv0/cs*exp(-(ks*cs)^2)).^2);

 St0 = 1/abs(1 + 8*Rv0/(cs*Ft))^2;
Tf = 1 - St/St0;

Rvt = Rvi + (Rv0 - Rvi) * Tf;
Rht = Rhi + (Rh0 - Rhi) * Tf;

fvv = 2*Rvt/cs;
fvvo = 2*Rv0/cs;
fvvi = 2*Rvi/cs;
fhh = -2*Rht/cs;
fhho = -2*Rh0/cs;
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

Tvo = 1 + Rv0;
Tvmo = 1 - Rv0;
Tho = 1 + Rh0;
Thmo = 1 - Rh0;

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

rss = sqrt( 2*sigma12/L1^2 + 2*sigma22/L2^2 );
rslp = rss;

ctorslp = ct/(sqrt(2)*rslp);
ctsorslp = cts/(sqrt(2)*rslp);
shadf = 1/2*(1/(sqrt(pi)*ctorslp) * exp(-ctorslp^2) - erfc(ctorslp));
shadfs = 1/2*(1/(sqrt(pi)*ctsorslp) * exp(-ctsorslp^2) - erfc(ctsorslp));
ShdwS = 1/(1 + shadf + shadfs);

ampterm = ShdwS * k/4 * exp(-2*sigma^2*kz^2);
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

out = [ssv, ssh, ossv, ossh, issv, issh, Rvt, Rht];