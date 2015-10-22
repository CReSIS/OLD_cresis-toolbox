% Script Information/Details
% =============================================
% Sync ATM Data and a GeoSearch CSV
% Author: Kyle Purdon
% Version: V1 01/09/2012
% Contact Information: cresis_data@cresis.ku.edu
% Contact Subject: "atm sync to geoscv"
% Toolboxes Used: 
% Known Bugs: None
% Planned Updates: None
% Additional Information:
% =============================================

%% User Input

% CSV File Input (With filename and extension)
% (Created in geographic_Search_gui.m)
param.csvinput = 'C:\Users\SomeUser\SomeInputFile.csv';

% Path to NASA ATM Data (Local to CReSIS on projects/metadata/ATM_smooth_nadir/)
param.atmpath = '\\titan\projects\metadata\ATM_smooth_nadir\';
% param.atmpath = '/cresis/projects/metadata/ATM_smooth_nadir/';

% Max Gap To Interpolate in ATM Data (m)
param.maxgap = 500;

% Output File Path (With filename and extension)
param.outpath = 'C:\Users\SomeUser\SomeOutputFile.csv';

% Debugging Level (Optional) (Default = 0 "Off")
% Level-II Creates a plot for each result
% param.debug_level = 2;

%% Automated Section

atm_sync_to_geocsv(param);
clear param;