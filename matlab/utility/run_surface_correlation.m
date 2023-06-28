% Script Information\Details
% =============================================
% Finds the statistical correlation between the accumulation values and
%      their corresponding bed values.
% Author: Jilu Li, Steven Foga
% Version: 04\05\2012
% Contact Information: cresis_data@cresis.ku.edu
% Contact Subject: "surface correlation tool"
% Toolboxes Used: 
% Known Bugs: None
% Planned Updates: None
% Additional Information: 
% See also: surface_correlation.m
% =============================================

%% User input

% mcords layer data file(s)
% Example: depthData{end+1} = 'X:\mdce\mcords2\2011_Antarctica_TO\CSARP_layerData\20111201_02\Data_20111201_02_003.mat';
depthData = {};
depthData{end+1} = 'X:\mdce\mcords2\2011_Antarctica_TO\CSARP_layerData\20111201_02\Data_20111201_02_003.mat';

% accum layer data file(s)
% Example: accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_017.mat'
accumData = {};
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_017.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_018.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_019.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_020.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_021.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_022.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_023.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_024.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_025.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_026.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_027.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_028.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_029.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_028.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_030.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_031.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_032.mat';
accumData{end+1} = 'X:\mdce\accum\2011_Antarctica_TO\CSARP_layerData\20111202_01\Data_20111202_01_033.mat';

surface_correlation(depthData,accumData)