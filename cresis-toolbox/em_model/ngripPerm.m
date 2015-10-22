function [depth,er] = ngripPerm(freq,numPoints,plotFlag)
% [depth,er] = ngripPerm(freq,numPoints,plotFlag);
%
% freq = frequency in Hz
% numPoints = number of depth points [optional]
% plotFlag = binary [optional]
%
% depth = depth profile (m)
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

if ~exist('numPoints','var')
  numPoints = 10001;
end

if ~exist('plotFlag','var')
  plotFlag = 0;
elseif plotFlag
  clf;
end

% Create discrete depth profile
depth = linspace(0,3085,numPoints).';
depthInt = (depth(1:end-1) + depth(2:end))/2;

% Get acidity profile estimate from GRIP (activation energy of 0.22)
acid = load(fullfile(gRadar.data_support_path,'em_model_data','grip_dep_20030914.txt'));
eV = 0.22;
acidInterp = 1e-6 * interp1(acid(:,1),acid(:,3),depthInt,'linear','extrap');
condInterp = 1e-6 * interp1(acid(:,1),acid(:,2),depthInt,'linear','extrap');
if (plotFlag)
  plot(depthInt,acidInterp,'k-');
  old = axis; axis([0 3085 old(3:4)]);
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

% Get NGRIP's temperature profile
temp = load(fullfile(gRadar.data_support_path,'em_model_data','ngrip_temp_20031009_Lars.txt'));

index = 1;
meas_index = 1;
while index <= size(temp,1)
   num = 1;
   tmp = temp(index,2);
   while index+1 <= size(temp,1) & temp(index,1) == temp(index+1,1)
      index = index + 1;
      num = num + 1;
      tmp = tmp + temp(index,2);
   end
   temp_cl(meas_index,1) = temp(index,1);
   temp_cl(meas_index,2) = tmp/num;
   meas_index = meas_index + 1;
   index = index + 1;
end

[temp_cl(:,1) indexes] = sort(temp_cl(:,1));
temp_cl(:,2) = temp_cl(indexes,2);
temp_clean = linspace(0,temp_cl(end,1),201).';
temp_clean(:,2) = interp1(temp_cl(:,1),temp_cl(:,2),temp_clean(:,1),'linear','extrap');

tempInterp = interp1(temp_clean(:,1),temp_clean(:,2),depthInt,'linear','extrap');
if (plotFlag)
  plot(depthInt,tempInterp,'k-');
  old = axis; axis([0 3085 old(3:4)]);
  xlabel('Depth (meters)');
  ylabel('Temperature (Celcius)');
  fprintf('Mean temp = %f C\n', mean(tempInterp));
  pause;
end

% Get NGRIP's density profile
file_den = load(fullfile(gRadar.data_support_path,'em_model_data','ngrip_density_20030914.txt'));

new_index = 1;
den(new_index,1:2) = file_den(1,:);
for index = 2:size(file_den,1)
   if (file_den(index-1,1) == file_den(index,1))
      den(new_index,:) = file_den(index,:);
   else
      new_index = new_index + 1;
      den(new_index,:) = file_den(index,:);
   end
end
den = [den(2:end,:); 3085 917];
denInterp = interp1(den(:,1),den(:,2)/1000,depthInt,'linear','extrap');
if (plotFlag)
  plot(depthInt,denInterp,'k-');
  old = axis; axis([0 3085 old(3:4)]);
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

% From DEP conductivity profile
er = iceCond(273.15+tempInterp,denInterp,freqRep,condInterp,eV,273.15-15);
% er = iceAcid(273.15+tempInterp,denInterp,freqRep,acidInterp,eV);
% er = iceLFCond(273.15+tempInterp,denInterp,freqRep,condInterp,eV,273.15-15);
% er = iceWestphal(273.15+tempInterp,denInterp,freqRep);
w = 2*pi*freqRep;
sigma = 0;
if (plotFlag)
  atten = real(j*sqrt(-j*w.*u0.*(sigma+j*w.*er*e0)));
  range = sqrt(470^2+3085^2);
  gain = 10.^(mean(20*log10(exp(-atten .* 2*range)))/10);
  plot(depthInt,-20*log10(exp(-2*atten*1000)),'k-')
  old = axis; axis([0 3085 old(3:4)]);
  xlabel('Depth (m)');
  ylabel('Two-way Loss (dB/km)');
  fprintf('Mean two-way loss = %f dB/km\n', mean(-20*log10(exp(-2*atten*1000))));
  fprintf('Total loss = %f dB\n', mean(-20*log10(exp(-2*atten*1000)))*range/1000);
  pause;
  subplot(3,1,1);
  plot(depthInt,tempInterp,'k-');
  old = axis; axis([0 3085 old(3:4)]);
  title('(a)');
  ylabel('Temperature (Celcius)');
  subplot(3,1,2);
  plot(depthInt,denInterp,'k-');
  old = axis; axis([0 3085 old(3:4)]);
  title('(b)');
  ylabel('Density (g/cm^3)');
  fprintf('Mean density = %f g/cm^3\n', mean(denInterp));
  subplot(3,1,3);
  plot(depthInt,1e6*condInterp,'k-');
  old = axis; axis([0 3085 old(3:4)]);
  title('(c)');
  xlabel('Depth (meters)');
  ylabel('Conductivity (uS/m)');
  fprintf('Mean conductivity = %f uS/m\n', mean(condInterp)/1e-6);

  % Specially Formatted Graph for GRSL publication
  clf;
  ax1 = axes;
  hl1 = line(depthInt,tempInterp,'Color','k','Parent',ax1);
  axis([0 3085 -140 0]);
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
  axis([0 3085 -0.1 1.4]);
  ax3 = axes('Position',get(ax1,'Position'),...
             'XAxisLocation','bottom',...
             'YAxisLocation','left',...
             'Color','none',...
             'XColor','k','YColor','k');
  hl3 = line(depthInt,1e6*condInterp,'Color','k','Parent',ax3);
  ylabel('Conductivity (uS/m)');
  axis([0 3085 0 140]);
  set(ax1,'YTick',[-40:10:0]);
  set(ax2,'YTick',[0.3:0.2:0.9]);
  set(ax3,'YTick',[0:15:45]);
  set(ax2,'Position',get(ax1,'Position'))
  set(ax3,'Position',get(ax1,'Position'))
  line([0 3085],[140 140],'Color','k','Parent',ax3);
end

return;

