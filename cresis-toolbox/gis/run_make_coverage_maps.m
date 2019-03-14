
% script run_make_coverage_maps

% Calls the make_coverage_maps function
%
% Author: John Paden, edits by Rohan Choudhari
%
% See also: make_coverage_maps

global gRadar

%% User Settings
% =========================================================================

% --------- STEP 1: SET THE LOCATION ---------
location = 'Arctic';
% location = 'Antarctica';




% --------- STEP 2: SELECTING THE SEASON(S) ---------
params = {};
if strcmpi(location,'Arctic')
  
  %params{end+1} = read_param_xls(ct_filename_param('rds_param_2008_Greenland_TO.xls'));
  %params{end+1} = read_param_xls(ct_filename_param('rds_param_2009_Greenland_TO.xls'));
%   params{end+1} = read_param_xls(ct_filename_param('rds_param_2010_Greenland_DC8.xls'));
%   params{end+1} = read_param_xls(ct_filename_param('rds_param_2010_Greenland_P3.xls'));
%   params{end+1} = read_param_xls(ct_filename_param('rds_param_2011_Greenland_P3.xls'));
  %params{end+1} = read_param_xls(ct_filename_param('rds_param_2011_Greenland_TO.xls'));
  %params{end+1} = read_param_xls(ct_filename_param('rds_param_2012_Greenland_P3.xls'));
  %params{end+1} = read_param_xls(ct_filename_param('rds_param_2013_Greenland_P3.xls'));
  params{end+1} = read_param_xls(ct_filename_param('rds_param_2015_Greenland_C130.xls'));
  %params{end+1} = read_param_xls(ct_filename_param('rds_param_2015_Greenland_Polar6.xls'));
%   params{end+1} = read_param_xls(ct_filename_param('rds_param_2016_Greenland_P3.xls'));
  %params{end+1} = read_param_xls(ct_filename_param('rds_param_2016_Greenland_Polar6.xls'));
  %params{end+1} = read_param_xls(ct_filename_param('rds_param_2017_Greenland_P3.xls'));
  %params{end+1} = read_param_xls(ct_filename_param('rds_param_2018_Greenland_Polar6.xls'));
  
elseif strcmpi(location,'Antarctica')
else
  return
end




% --------- STEP 3: SELECT THE DATA SOURCE ---------
% string 'gps' or 'records' (which files in csarp_support to use)
%user_variables.data_source = 'records';
user_variables.data_source = 'gps';




% --------- STEP 4 : SELECT THE GIS FILE  ---------
geotiff_fn = ct_filename_gis(gRadar,fullfile('arctic','NaturalEarth_Data','Arctic_NaturalEarth.tif'));
% geotiff_fn = ct_filename_gis(gRadar,fullfile('arctic','Landsat-7','arctic_natural_90m.tif'));
% geotiff_fn = ct_filename_gis(gRadar,fullfile('greenland','Landsat-7','Greenland_natural.tif'));
% geotiff_fn = ct_filename_gis(gRadar,fullfile('canada','Landsat-7','Canada_90m.tif'));




% --------- STEP 5 : SELECT WHETHER TO SAVE THE FIGURE ---------
% set save image flag
user_variables.save_img = 1; %to save
% save_img = 0; % to not save




% --------- STEP 6 : SELECT THE FILE FORMAT ---------
user_variables.ext = '.m'; %Default
%user_variables.ext = '.jpg';
%user_variables.ext = '.fig';






% --------- STEP 7: SELECT THE DESTINATION WHERE THE FIGURE IS TO BE SAVED AT ---------
if user_variables.save_img == 1
  user_variables.out_dir = 'H:\rohan\coverage_maps\';
end




% --------- STEP 8: SET THE DESIRED AXIS LIMIT ---------
user_variables.axis_limits = [-2300 1000 -3500 1500]; % All of Arctic (use for OIB coverage maps)
% axis_limits = [-3000 1000 -1500 2500]; % All of Antarctica (use for OIB coverage maps)




% --------- STEP 9: ALTER OTHER OPTIONAL VALUES ---------

% dt = time spacing (sec) between each plotted point
user_variables.dt = 1;

% along_track_sampling: desired along track sampling of plot (m)
user_variables.along_track_sampling = 100;

% plot_args: cell array of argument to plot function (e.g. 'b.' for blue dots)
user_variables.plot_args = {'LineWidth',1};

% label: Add text labels to the plots for specific frames
%  .day_seg: cell list of day segments to match
%  .frm: equal length cell vector to .day_seg with specific frames to label
%  .text: equal length cell vector to .day_seg with string for text label
label.day_seg = {};
label.frm = {};
label.text = {};
% label.day_seg = {'20091102_02','20091016_01','20091016_01','20091102_02','20091102_02'};
% label.frm = {8, 21, 26, 23, 32};
% label.text = {' 1',' 2',' 3',' 4',' 5'};
user_variables.label = label;

% figure_position: Set figure Position property (leave empty to use default/current figure size)
user_variables.figure_position = [];
% figure_position = [221   -60   654   704];






%% Automated section
% =========================================================================

if ~exist(out_dir,'dir')
  mkdir(out_dir);
end
  

%Looping through all the selected parameter spreadsheets
for season_idx = 1:length(params)
  
  % Call coverage maps for each tmp_params
  
  %Set the params according to what you want
  params{season_idx} = ct_set_params(params{season_idx},'cmd.generic',0,'cmd.notes','do not process');
  fig_h_idx = season_idx;
  user_variables.season_idx = season_idx;
  user_variables.proj = plot_geotiff(geotiff_fn,[],[], fig_h_idx);
  test_make_coverage_maps(params{season_idx}, user_variables);
  
  
 
end
