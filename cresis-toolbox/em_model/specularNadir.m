function H = specularNadir(depth,er,freq)
% H = specularNadir(depth,er,freq)
%
% freq = frequency, must be a row vector (Hz)
%
% depth = depth profile (m) (column vector)
% er = relative permittivity profile er = e' - j*e''
%   2-D matrix with len(depth)-1 rows and len(freq) columns
%
% Depth Index
%     1       ---------------------------------- air/ice = 0
%                           er(1)
%     2       ----------------------------------
%                           er(2)
%     3       ----------------------------------
%     .
%     .
%     .
%    n-1      ----------------------------------
%                           er(n-1)
%     n       ---------------------------------- ice/bedrock
%                           er(n) (bedrock)
%
% Depth Index
%     1       ---------------------------------- antenna height
%                           er(1) = 1 + i0
%     2       ---------------------------------- air/ice = 0
%
% Because of refractionGain, you are restricted to ice-sheet
% type profiles.  If the antenna is above the air interface,
% then er(1) must be 1. If er(1) is not 1, then it is assumed
% that the antenna is sitting directly on the ice sheet.
%
% The air/ice interface must be set to a depth of zero.

format compact; format long;
if (ispc)
  path(path,'P:\prism\radar\radarSimulator\');
else
  path(path,'/projects/prism/radar/radarSimulator/');
end

physicalConstants;

% w = omega = angular frequency
w = 2*pi*freq;

% Determine path length through each layer
len = diff(depth);

% Set up initial conditions for for-loop
% g = gamma = propagation constant
gOld = j*w.*sqrt(u0*er(1,:)*e0);

% Determine refraction gain at each interface
if er(1) == 1
  % Ignore the air layer
  Gf = refractionGain(depth(2:end),real(er(2:end,1)),real(er(1)),0,depth(2:end));
  Gf(1) = 1;
else
  Gf = refractionGain(depth,real(er(:,1)),real(er(1)),0,depth(2:end));
end

H = complex(zeros([1 length(freq)]));
oldH = complex(ones([1 length(freq)]));
for n = 1:length(len)
  g = j*w.*sqrt(u0*er(n+1,:)*e0);

  A = gOld./g;
  refl = (A-1)./(A+1);
  tran12 = 2*A./(A+1);
  tran21 = 2./(A+1);

  % H = H + oldH.*refl.*exp(-2*j*imag(gOld)*len(n));
  % H = H + oldH.*refl.*exp(-2*gOld*len(n));
  lambda = 2*pi./imag(gOld);
  if 1
    % Includes effective area to gain conversion
    H = H + oldH.*Gf(n).*lambda./(4*pi*2*(depth(n+1)-depth(1))).*refl.*exp(-2*gOld*len(n));
  else
    % Does not include effective area to gain conversion
    H = H + oldH.*Gf(n)./sqrt(4*pi*(2*(depth(n+1)-depth(1))).^2).*refl.*exp(-2*gOld*len(n));
  end

  % oldH = oldH.*exp(-2*j*imag(gOld)*len(n));
  % oldH = oldH.*exp(-2*gOld*len(n));
  oldH = oldH.*tran12.*tran21.*exp(-2*gOld*len(n));

  % Recurse step
  gOld = g;
end

return;

% ------------------------------------------------------------------------------
% Examples
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% 1. Example of PRISM SAR RDS mode
format compact; format long; clear H;
if (ispc)
  path(path,'P:\prism\radar\radarSimulator\');
else
  path(path,'/projects/prism/radar/radarSimulator/');
end
physicalConstants;

freq = linspace(120e6,300e6,7200);
incSize = 500;
for ind = 1:ceil(length(freq)/incSize)
  start = (ind-1)*incSize+1;
  stop = min(ind*incSize,length(freq));
  fprintf('Indices %.0f - %.0f\n',start,stop);
  [depth,er] = gisp2Perm(freq(start:stop),4000);
  % [depth,er] = gisp2Perm(freq(start:stop),100);
  er(end+1,:) = real(rock(273.15,freq,9));
  H(start:stop) = specularNadir(depth,er,freq(start:stop));
  plot(20*log10(abs(ifft(H,10*length(H)))*10));
  % fprintf('Press Return\n'); pause;
end
fs = 180e6*10;
dt = 1/fs;
N = length(H)*10;
t = 0:dt:(N-1)*dt;
r = t.*c/1.78/2;
% +3 approximately accounts for hanning weighting
plot(r,20*log10(abs(ifft(H.'.*hanning(length(H)),10*length(H)))*10)+3,'k');
xlabel('Range (m)');
ylabel('Impulse Response (dB)');
axis([150 3050 -240 -100]);

load files60_79;
fs = 240e6;
dt = 1/fs;
N = size(iOut,1);
t = 0:dt:(N-1)*dt;
r = [t.*3e8/1.78/2].';
filtData = medfilt1(max(iOut(:,1:10)-100,[],2),20);
indices = find(r > 300 & r < 3000);
p = polyfit(r(indices),filtData(indices),4);
fitData = polyval(p,r);
h1 = plot(r,filtData,'k-');
hold on;
h2 = plot(r,fitData,'r-');
h3 = plot(r,fitData+40,'r:');
h4 = plot(r,fitData-40,'r-.');
hold off;
legend([h1 h2 h3 h4],'Data','Poly Fit','Upper Bound','Lower Bound');
axis([150 3050 -240 -100]);
xlabel('Range (m)');
ylabel('Impulse Response (dB)');
save('/users/paden/tmp/deleteThis2.mat','r','filtData','-V6');

plot(r,fitData-40,'k-');
hold on;
plot(r,fitData+40,'r-');
plot(r,filtData,'c--');
plot(r,fitData,'c--');
hold off;
legend('Loop Sensitivity','Max Power');
axis([150 3050 -240 -100]);
xlabel('Range (m)');
ylabel('P_t / P_r (dB)');

plot(r,fitData+40 - (fitData-40),'k-');
hold on;
plot(r,fitData+40 - max((fitData-40),-220),'r:');
hold off;
legend('Dynamic Range','Dynamic Range w/ 220 dB Limit');
xlabel('Range (m)');
ylabel('Dynamic Range (dB)');


% ------------------------------------------------------------------------------
% 2. Example of 2-interface (three media)
format compact; format long;
if (ispc)
  path(path,'P:\prism\radar\radarSimulator\');
else
  path(path,'/projects/prism/radar/radarSimulator/');
end

physicalConstants;

freq = linspace(0.1e6,300e6,201);
depth = [0 1 1+0.99930819333448/2].';
er = [1+j*0e-5./(2*pi*freq*e0)
      4+j*0e-5./(2*pi*freq*e0)
      1+j*0e-5./(2*pi*freq*e0)];
H = specularNadir(depth,er,freq);
plot(freq/1e6,20*log10(abs(H)));
xlabel('Frequency (MHz)');
ylabel('Reflection Coefficient (dB)');

% ------------------------------------------------------------------------------
% 3. Example of 5-interface (six media)
format compact; format long;
if (ispc)
  path(path,'P:\prism\radar\radarSimulator\');
else
  path(path,'/projects/prism/radar/radarSimulator/');
end

physicalConstants;

freq = linspace(120e6,300e6,201);
depth = [-2 0 2 4 6 8].';
er = [1
      4
      16
      64
      256
      1024];
er = repmat(er,[1 length(freq)]);
H = specularNadir(depth,er,freq);
plot(20*log10(abs(ifft(H,100*length(H)))*100));
xlabel('1/100 Range Bins (MHz)');
ylabel('Reflection Coefficient (dB)');

% ------------------------------------------------------------------------------
% 4. Example of 2-interface (three media)
format compact; format long;
if (ispc)
  path(path,'P:\prism\radar\radarSimulator\');
else
  path(path,'/projects/prism/radar/radarSimulator/');
end

physicalConstants;

freq = linspace(120e6,300e6,201);
depth = [0 500 500+1].';
er = [1.78^2+j*1e-5./(2*pi*freq*e0)
      1.78^2+j*2e-5./(2*pi*freq*e0)
      1.78^2+j*1e-5./(2*pi*freq*e0)];
H = specularNadir(depth,er,freq);
H = specular(depth,er,0,freq);
plot(20*log10(abs(ifft(H,10*length(H)))*10));
% plot(freq/1e6,20*log10(abs(H)));
xlabel('Frequency (MHz)');
ylabel('Reflection Coefficient (dB)');

freq = linspace(120e6,300e6,201);
er1 = 1.78^2+j*1e-5./(2*pi*freq*e0);
eta1 = sqrt(u0./(e0*er1));
er2 = 1.78^2+j*2e-5./(2*pi*freq*e0);
eta2 = sqrt(u0./(e0*er2));
plot(freq/1e6,20*log10(abs((eta2-eta1)./(eta2+eta1))))

% ------------------------------------------------------------------------------
% 4. Example of gisp2Perm with 640-660 MHz radar
format compact; format long;
if (ispc)
  path(path,'P:\prism\radar\radarSimulator\');
else
  path(path,'/projects/prism/radar/radarSimulator/');
end

physicalConstants;

% Radar parameters
nT = 800;           % Number of Samples
BW = 20e6;          % Bandwidth
tS = 1/BW;          % Time Step, sampling interval
T = nT*tS;          % Sample period
fStep = 1/T;        % Frequency step
fc = 650e6;         % Center frequency
f1 = fc - BW/2;     % Initial frequency

freq = f1 + [0:fStep:(nT-1)*fStep];

tAxis = [0:tS/10:(10*nT-1)*tS/10].';

[depth,er] = gisp2Perm(freq,250);
% Add air layer
depth = [-200; depth];
er = [ones(1,size(er,2)); er];
% Add rock dielectric
er(end+1,:) = 8*ones(1,size(er,2));
H = specularNadir(depth,er,freq);
figure(4);
plot(tAxis*1e6,20*log10(abs(ifft(H.*hanning(nT).',10*length(H)))*10)+3);
% plot(freq/1e6,20*log10(abs(H)));
xlabel('Time (us)');
ylabel('Reflection Coefficient (dB)');



