function [depth,er] = summitPermSimp(numPoints,plotFlag,basePath)
% [depth,er] = summitPermSimp(numPoints,plotFlag,basePath);
%
% numPoints = number of depth points [optional]
% plotFlag = binary [optional]
% basePath = base path to geophysical profiles used by the function
%
% depth = depth profile (m) (column vector)
% er = relative permittivity profile er = e' (real-only)
%   vector with len(depth)-1 elements
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
%
% This SIMPLE permittivity profile is used for SAR processing
% and f-k migration. It is purely real. Assumes that the real
% part is not a function of frequency and doesn't change beyond
% a certain depth (340 m?).

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
depth_extent = 4000;
depth = linspace(0,depth_extent,numPoints).';
depthInt = (depth(1:end-1) + depth(2:end))/2;

%-------------------------------------------------------------------------------
% Get density profile estimate from GISP2
%   Valid from 0 m to 1511 m and 94 m to 333 m and 0 to 15 m
%-------------------------------------------------------------------------------
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
  old = axis; axis([0 depth_extent old(3:4)]);
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
er = iceSummit(denInterp);

return;

%-------------------------------------------------------------------------------
% Example
%-------------------------------------------------------------------------------

[depth,er] = summitPermSimp(4001);
depthInt = (depth(1:end-1) + depth(2:end))/2;
plot(depthInt,er);

