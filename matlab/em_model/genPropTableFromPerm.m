function [time,gain,txAngle,xpos] = genPropTableFromPerm(depth,er,freq,inc);
% [time,gain,txAngle,xpos] = genPropTableFromPerm(depth,er,freq,inc);
% depth = depth axis (m)
% e0*er = e0*(e' - j*e'') (F/m)
% freq = frequency (Hz)
% inc = incidence angle in the first layer (rad)
% 
% depth must be a column vector
% er must be 2-D matrix
%   size(er,1) = length(depth)-1 or length(depth) (last value not used)
%   size(er,2) = size(freq,2)
% freq must be a N-D matrix
%   frequency should only change along column/2nd dimension because the
%   er's are matched to freq along this dimension.
% inc must be a N-D matrix
%   size(inc) = size(freq)
% 
% time = two-way time (sec)
% gain = two-way extinction through the ice (ratio of fields)
% txAngle = transmission angle in the last layer (rad)
% xpos = x-offset from nadir assuming ray-tracing (m)
%
% time, gain, txAngle, and xpos will all be the same size as freq and inc
%
% Depth Index
%     1       ----------------------------------
%                           er(1)
%     2       ----------------------------------
%                           er(2)
%     3       ----------------------------------
%                           er(3)
%     4       ----------------------------------
%     .
%     .
%     .
%    n-1      ----------------------------------
%                           er(n-1)
%     n       ----------------------------------
%                           er(n) <-- optional (not used)
%
% Antenna is assumed to be a point source on the depth index 1 interface.
% The inc argument gives the direction of propagation in er(1).
% The time and gain are calculated after going through er(1) through er(n).
% The tran is the direction of propagation in er(n).
% The xpos is the x-offset of the ray on the depth index n interface.
%
% NEED TO ADD EFFECT OF TRANSMISSIVITY LOSS!!!!!

physical_constants;
w = 2*pi*freq;

% Preallocate memory
time = zeros(size(inc));
gain = ones(size(inc));
txAngle = zeros(size(inc));
xpos = zeros(size(inc));

% gamma is the propagation constant
gamma = j*sqrt(-j*w.*u0.*(j*w.*er(1)*e0));
% beta is the phase constant
betaLayer1 = imag(gamma);

for n = 1:length(depth)-1
  % gamma is the propagation constant
  gamma = j*sqrt(-j*w.*u0.*(j*w.*er(n)*e0));
  % alpha is the attenuation constant
  alpha = real(gamma);
  % beta is the phase constant
  beta = imag(gamma);
  % velocity of propagation
  lambda = 2*pi./beta;
  vel = freq.*lambda;
  betaZ = sqrt(beta.^2 - (betaLayer1.*sin(inc)).^2);

  % Thickness of current layer
  thick = depth(n+1)-depth(n);

  txAngle = acos(betaZ./beta);
  alphaPrime = alpha./cos(txAngle);
  time = time + 2*(thick./cos(txAngle))./vel;
  gain = gain.*exp(-4*alphaPrime.*thick);
  xpos = xpos + thick.*tan(txAngle);
end

return;

% ------------------------------------------------------------------------------
% Example
% ------------------------------------------------------------------------------
depth = [0 1 2 3];
er = [1 1 1];
freq = [150e6 300e6; 150e6 300e6];
inc = [0 0; 30 30]/180*pi;
[time,gain,txAngle,xpos] = genPropTableFromPerm(depth,er,freq,inc);

freq = linspace(140e6,160e6,1024);
[depth,er] = gisp2Perm(freq,10001);
inc = linspace(0,pi/4,11).';
incMat = repmat(inc,[1 length(freq)]);
freqMat = repmat(freq,[length(inc) 1]);
[time,gain,txAngle,xpos] = genPropTableFromPerm(depth,er,freq,inc);

freq = 210e6;
[depth,er] = gisp2Perm(freq,10001);
inc = linspace(0,pi/4,401).';
[time,gain,txAngle,xpos] = genPropTableFromPerm(depth,er,freq,inc);

figure(1)
plot(xpos,time*1e6,'k-');
xlabel('Off Nadir (m)');
ylabel('Two-way Time Delay (us)');
title('a)');
figure(2)
plot(xpos,inc*180/pi,'k-');
xlabel('Off Nadir (m)');
ylabel('Initial Transmission Angle \theta_t (deg)');
title('b)');
figure(3)
plot(xpos,txAngle*180/pi,'k-');
xlabel('Off Nadir (m)');
ylabel('Bed Incidence Angle \theta_i (deg)');
title('c)');
figure(4)
plot(xpos,10*log10(abs(gain)),'k-');
xlabel('Off Nadir (m)');
ylabel('Two-way Attenuation L_t (dB)');
title('d)');




