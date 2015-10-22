function [depth,er] = gripPerm(freq,depth,plotFlag)
% [depth,er] = gripPerm(freq,depth,plotFlag);
%
% freq = frequency in Hz (must be a row vector)
% depth = depth vector [optional]
%   defaults to 10001 points, zero to 3400 meters
% plotFlag = binary [optional]
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
%                           er(3)
%     4       ----------------------------------
%     .
%     .
%     .
%    n-1      ----------------------------------
%                           er(n-1)
%     n       ---------------------------------- ice/bedrock
%
% Note that er is interpolated midpoint value of each layer.

format compact;
physical_constants;
global gRadar;

if ~exist('depth','var') || isempty(depth)
  numPoints = 10001;
  depth_extent = 3400;
  depth = linspace(0,depth_extent,numPoints).';  
end

if ~exist('plotFlag','var')
  plotFlag = 0;
elseif plotFlag
  clf;
end

depthInt = (depth(1:end-1) + depth(2:end))/2; 

%-------------------------------------------------------------------------------
% Get GRIP's acidity profile (activation energy of 0.22)
%   Valid from 138 m to 3048 m
%-------------------------------------------------------------------------------
acid = load(fullfile(gRadar.data_support_path,'em_model_data','grip_dep_20030914.txt'));

% Extend the initial value of acid/conductivity to 0 m
%  - Uses the mean of the first 100 m
initCond = mean(acid(find(acid(:,1)-acid(1,1) < 100),2));
initAcid = mean(acid(find(acid(:,1)-acid(1,1) < 100),3));
acid = [[acid(1,1)-1e-5 initCond initAcid acid(1,4)-1e-5]; acid];
acid = [[0 initCond initAcid 0]; acid];
% Extend the final value of acid/conductivity to depth(end)
%  - Uses the mean of the last 100 m
finalCond = mean(acid(find(acid(end,1)-acid(:,1) < 100),2));
finalAcid = mean(acid(find(acid(end,1)-acid(:,1) < 100),3));
acid = [acid; [acid(end,1)+1e-5 finalCond finalAcid acid(end,1)+1e-5]];
acid = [acid; [depth(end) finalCond finalAcid 0]];

eV = 0.22;
acidInterp = 1e-6 * interp1(acid(:,1),acid(:,3),depthInt,'linear','extrap');
condInterp = 1e-6 * interp1(acid(:,1),acid(:,2),depthInt,'linear','extrap');
if (plotFlag)
  plot(depthInt,1e6*acidInterp,'k-');
  old = axis; axis([0 3027.6 old(3:4)]);
  xlabel('Depth (meters)');
  ylabel('Acidity (micromolarity)');
  fprintf('Mean acidity = %f micromolarity\n', mean(acidInterp)/1e-6);
  pause;
  plot(depthInt,1e6*condInterp,'k-');
  old = axis; axis([0 3027.6 old(3:4)]);
  xlabel('Depth (meters)');
  ylabel('Conductivity (uS/m)');
  fprintf('Mean conductivity = %f uS/m\n', mean(condInterp)/1e-6);
  pause;
end

%-------------------------------------------------------------------------------
% Get GRIP's temperature profile
%   Valid from 1.1 m to 3027.6 m
%-------------------------------------------------------------------------------
temp= load(fullfile(gRadar.data_support_path,'em_model_data','grip_temp_20030914.txt'));

tempInterp = interp1(temp(:,1),temp(:,2),depthInt,'linear','extrap');
if (plotFlag)
  plot(depthInt,tempInterp,'k-');
  old = axis; axis([0 3027.6 old(3:4)]);
  xlabel('Depth (meters)');
  ylabel('Temperature (Celcius)');
  fprintf('Mean temp = %f C\n', mean(tempInterp));
  pause;
end

%-------------------------------------------------------------------------------
% Get density profile estimate from GISP2
%   Valid from 0 m to 1511 m and 94 m to 333 m and 0 to 15 m
%-------------------------------------------------------------------------------
% Load B-core data
file_den = load(fullfile(gRadar.data_support_path,'em_model_data','gisp2_density_0_200_B.txt'));

denDepth = [file_den(1:95,1); file_den(112:124,1)];
den = [file_den(1:95,4); file_den(112:124,2)];

% Extend the final value of density to depth(end)
%  - Use the last value (0.917)
% finalDen = mean(den(find(denDepth(end)-denDepth < 100)));
finalDen = 0.917;
denDepth = [denDepth; depth(end)];
den = [den; finalDen];

% Load snow-pit data for initial part of density profile
[highResDepth,highResDen] = loadGispSnowpitDensity;
indexes = find(denDepth <= highResDepth(end));
denDepth = [highResDepth.'; denDepth(indexes(end)+1:end)];
den = [highResDen.'; den(indexes(end)+1:end)];

% Interpolate
denInterp = interp1(denDepth,den,depthInt,'linear','extrap');
if (plotFlag)
  plot(depthInt,denInterp,'k-');
  old = axis; axis([0 3027.6 old(3:4)]);
  xlabel('Depth (meters)');
  ylabel('Density (g/cm^3)');
  fprintf('Mean density = %f g/cm^3\n', mean(denInterp));
  pause;
  axis([0 300 0 1]);
  pause;
end

%-------------------------------------------------------------------------------
% Determine ice permittivity profile from geophysical profiles
%-------------------------------------------------------------------------------
tempInterp = repmat(tempInterp,[1 length(freq)]);
denInterp = repmat(denInterp,[1 length(freq)]);
condInterp = repmat(condInterp,[1 length(freq)]);
freqRep = repmat(freq,[length(depthInt) 1]);

er = iceCond(273.15+tempInterp,denInterp,freqRep,condInterp,eV,273.15-15);
% er = iceAcid(273.15+tempInterp,denInterp,freqRep,acidInterp,eV);
% er = iceLFCond(273.15+tempInterp,denInterp,freqRep,condInterp,eV,273.15-15);
% er = iceWestphal(273.15+tempInterp,denInterp,freqRep);
w = 2*pi*freqRep;
sigma = 0;
if (plotFlag)
  atten = real(j*sqrt(-j*w.*u0.*(sigma+j*w.*er*e0)));
  range = 3027.6;
  gain = 10.^(mean(20*log10(exp(-atten .* 2*range)))/10);
  plot(depthInt,-20*log10(exp(-2*atten*1000)),'k-')
  old = axis; axis([0 3027.6 old(3:4)]);
  xlabel('Depth (m)');
  ylabel('Two-way Loss (dB/km)');
  fprintf('Mean two-way loss = %f dB/km\n', mean(-20*log10(exp(-2*atten*1000))));
  fprintf('Total loss = %f dB\n', mean(-20*log10(exp(-2*atten*1000)))*range/1000);
end

if (plotFlag)
  % Specially Formatted Graph for GRSL publication
  clf;
  ax1 = axes;
  hl1 = line(depthInt,tempInterp,'Color','k','Parent',ax1);
  axis([0 3027.6 -140 0]);
  xlabel('Depth (m)');
  ylabel('Temperature (C)');
  set(ax1,'XColor','k','YColor','k')
  ax2 = axes('Position',get(ax1,'Position'),...
             'XAxisLocation','bottom',...
             'YAxisLocation','right',...
             'Color','none',...
             'XColor','k','YColor','k');
  hl2 = line(depthInt,denInterp,'Color','k','Parent',ax2);
  ylabel('Density (g/cm^3)');
  axis([0 3027.6 -0.1 1.4]);
  ax3 = axes('Position',get(ax1,'Position'),...
             'XAxisLocation','bottom',...
             'YAxisLocation','left',...
             'Color','none',...
             'XColor','k','YColor','k');
  hl3 = line(depthInt,1e6*condInterp,'Color','k','Parent',ax3);
  ylabel('Conductivity (uS/m)');
  axis([0 3027.6 0 140]);
  set(ax1,'YTick',[-40:10:0]);
  set(ax2,'YTick',[0.3:0.2:0.9]);
  set(ax3,'YTick',[0:15:45]);
  set(ax2,'Position',get(ax1,'Position'))
  set(ax3,'Position',get(ax1,'Position'))
  line([0 3027.6],[140 140],'Color','k','Parent',ax3);
end

return;

