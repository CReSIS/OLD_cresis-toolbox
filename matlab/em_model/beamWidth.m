function [theta,xwidth] = beamWidth(freq,L,depth,er)
% Determines null-null half-beam-width of antenna array on top
% of planarly stratified media.  Also returns footprint on bedrock.
%
% freq = frequency of operation (Hz) (column vector)
% L = length of array (m)
% depth = depth profile (m)
% er = complex permittivity profile of stratified media
%     e0*er = e0 * (e' - j*e'')
%     (format described below)
%
% theta = 1/2 null-null beamwidth -- at the antenna (rad)
% xwidth = 1/2 foorprint on bedrock (m)
%
% depth must be a column vector
% er must be 2-D matrix
%   size(er,1) = length(depth)-1 or length(depth) (last value not used)
%   size(er,2) = size(freq,2)
% freq must be a row vector
%
% ------------------------------------------------------------------------------
% depth, er format
% ------------------------------------------------------------------------------
% depth Index
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
%     n       ---------------------------------- ice/bed
%   
%  length(depth) = length(er)+1

% Generate propagation table
inc = linspace(0, +pi/2 * 0.95, 501).';
freqMat = repmat(freq,[length(inc) 1]);
incMat = repmat(inc,[1 length(freq)]);
[time,gain,txAngle,xpos] = genPropTableFromPerm(depth,er,freqMat,incMat);

% Convert two-way time to one-way time (comment next line if considering SAR techniques!)
% time = time/2;

for ind = 1:length(freq)
  % For eachinc frequency, find the lambda/2 phase (i.e. halfWaveTime)
  timeOffset = interp1(xpos(:,ind),time(:,ind),xpos(:,ind)+L,'spline','extrap') - time(:,ind);
  plot(timeOffset);
  halfWaveTime = 0.5/freq(ind)*sqrt(real(er(1,ind)));
  xwidth(ind) = interp1(timeOffset,xpos(:,ind),halfWaveTime);
  theta(ind) = interp1(timeOffset,inc,halfWaveTime);
end

return;

% ------------------------------------------------------------------------------
% Examples
% ------------------------------------------------------------------------------

% PRISM SAR Example
if (ispc)
  path(path,'P:\prism\radar\radarSimulator\');
else
  path(path,'/projects/prism/radar/radarSimulator/');
end

freq = linspace(50e6,250e6,21);
[depth,er] = gisp2Perm(freq,1001);
[beamWidthTheta,beamWidthPos] = beamWidth(freq,4,depth,er);

plot(freq/1e6,beamWidthPos*2,'k-');
xlabel('Frequency (MHz)');
ylabel('Footprint (m)');
axis([50 250 300 1000]);

plot(freq/1e6,2 * beamWidthTheta*180/pi);
xlabel('Frequency (MHz)');
ylabel('Beamwidth (deg)');

% Another Test Example
freq = linspace(30e6,300e6,21);
depth = [0 3000].';
er = [1];
erMat = repmat(er,[1 length(freq)]);
[beamWidthTheta,beamWidthPos] = beamWidth(freq,4,depth,erMat);

plot(freq/1e6,2 * beamWidthTheta*180/pi);
xlabel('Frequency (MHz)');
ylabel('Beamwidth (deg)');

% plot(freq/1e6,beamWidthPos);
% xlabel('Frequency (MHz)');
% ylabel('Footprint (m)');

% Analytic solution:
c = 3e8;
L = 4;
hold on;
plot(freq/1e6,2 * asin(c./(2*L*freq*sqrt(er)))*180/pi,'r');
hold off;



