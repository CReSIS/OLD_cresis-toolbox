function [depth,er] = waisPerm(freq,depth,plotFlag)
% [depth,er] = waisPerm(freq,depth,plotFlag);
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
% Get GISP2's dep profile (activation energy of 0.22)
%   Valid from 138 m to 3048 m
%-------------------------------------------------------------------------------
acid = load(fullfile(gRadar.data_support_path,'em_model_data','domeC_dep_20111230.txt'));
% Just get the DEP profile
acid = acid(:,[1 3]);

% Extend the DEP to zero depth
acid = [[0 acid(1,2)]; acid];

% Extend the final value of acid/conductivity to depth(end)
%  - Uses the mean of the last 100 m
finalCond = mean(acid(find(acid(end,1)-acid(:,1) < 100),2));
acid = [acid; [depth(end) finalCond]];

eV = 0.22;
condInterp = interp1(acid(:,1),acid(:,2),depthInt,'linear','extrap');
if (plotFlag)
  plot(depthInt,1e6*condInterp,'k-');
  old = axis; axis([0 depth(end) old(3:4)]);
  xlabel('Depth (meters)');
  ylabel('Conductivity (uS/m)');
  fprintf('Mean conductivity = %f uS/m\n', mean(condInterp)/1e-6);
  pause;
end

%-------------------------------------------------------------------------------
% Get WAIS's temperature profile
%   Valid from 100 m to 3048 m
%-------------------------------------------------------------------------------
temp = load(fullfile(gRadar.data_support_path,'em_model_data','wais_density_temp_20111229.txt'));
% Just get the temperature profile
temp = temp(:,[1 2]);

% Extend the final value of temperature to depth(end)
temp = [temp; [depth(end) temp(end,2)]];

tempInterp = interp1(temp(:,1),temp(:,2),depthInt,'linear','extrap');
if (plotFlag)
  plot(depthInt,tempInterp,'k-');
  old = axis; axis([0 3047.9 old(3:4)]);
  xlabel('Depth (meters)');
  ylabel('Temperature (Celcius)');
  fprintf('Mean temp = %f C\n', mean(tempInterp));
  pause;
end

%-------------------------------------------------------------------------------
% Get density profile estimate from GISP2
%   Valid from 0 m to 1511 m and 94 m to 333 m and 0 to 15 m
%-------------------------------------------------------------------------------
den = load(fullfile(gRadar.data_support_path,'em_model_data','wais_density_temp_20111229.txt'));
% Just get the density profile
den = den(:,[1 3]);
den(:,2) = den(:,2)/1000;

% Extend the final value of temperature to depth(end)
den = [den; [depth(end) den(end,2)]];

% Interpolate
denInterp = interp1(den(:,1),den(:,2),depthInt,'linear','extrap');
if (plotFlag)
  plot(depthInt,denInterp,'k-');
  old = axis; axis([0 3047.9 old(3:4)]);
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
w = 2*pi*freq;
sigma = 0;

tempInterp = repmat(tempInterp,[1 length(freq)]);
denInterp = repmat(denInterp,[1 length(freq)]);
condInterp = repmat(condInterp,[1 length(freq)]);
freqRep = repmat(freq,[length(depthInt) 1]);

er = iceCond(273.15+tempInterp,denInterp,freqRep,condInterp,eV,273.15-15);
% er = iceLFCond(273.15+tempInterp,denInterp,freqRep,condInterp,eV,273.15-15);
% er = iceWestphal(273.15+tempInterp,denInterp,freqRep);
% er = iceAcid(273.15+tempInterp,denInterp,freqRep,acid_interp,eV);
if (plotFlag)
  atten = real(j*sqrt(-j*w.*u0.*(sigma+j*w.*er*e0)));
  range = 3047.9;
  gain = 10^(mean(20*log10(exp(-atten * 2*range)))/10);
  plot(depthInt,-20*log10(exp(-2*atten*1000)),'k-');
  old = axis; axis([0 3047.9 old(3:4)]);
  xlabel('Depth (m)');
  ylabel('Two-way Loss (dB/km)');
  fprintf('Mean two-way loss = %f dB/km\n', mean(-20*log10(exp(-2*atten*1000))));
  fprintf('Total loss = %f dB\n', mean(-20*log10(exp(-2*atten*1000)))*range/1000);
  pause;
  plot(depthInt,real(er));
  xlabel('Depth (m)');
  ylabel('Real Permittivity e''');
  pause;
  plot(depthInt,-imag(er));
  xlabel('Depth (m)');
  ylabel('Imag Permittivity e''''');
end

if plotFlag
  % Specially Formatted Graph
  clf;
  ax1 = axes;
  hl1 = line(depthInt,tempInterp,'Color','k','Parent',ax1);
  axis([0 3047.9 -140 0]);
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
  axis([0 3047.9 -0.1 1.4]);
  ax3 = axes('Position',get(ax1,'Position'),...
             'XAxisLocation','bottom',...
             'YAxisLocation','left',...
             'Color','none',...
             'XColor','k','YColor','k');
  hl3 = line(depthInt,1e6*condInterp,'Color','k','Parent',ax3);
  ylabel('Conductivity (uS/m)');
  axis([0 3047.9 0 140]);
  set(ax1,'YTick',[-40:10:0]);
  set(ax2,'YTick',[0.3:0.2:0.9]);
  set(ax3,'YTick',[0:15:45]);
  set(ax2,'Position',get(ax1,'Position'))
  set(ax3,'Position',get(ax1,'Position'))
  line([0 3047.9],[140 140],'Color','k','Parent',ax3);
end

return;

%-------------------------------------------------------------------------------
% Example
%-------------------------------------------------------------------------------

freq = linspace(120e6,300e6,11);
[depth,er] = gisp2Perm(freq);
depthInt = (depth(1:end-1) + depth(2:end))/2;
plot(depthInt,real(er));
plot(depthInt,imag(er));

