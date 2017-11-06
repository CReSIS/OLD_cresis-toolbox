% script run_check_data_products
%
% Script for running check_data_products.m.
%
% Author: John Paden, Logan Smith

%% User Settings
clear; clc
global gRadar
% params = read_param_xls(ct_filename_param('kuband_param_2009_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('kuband_param_2010_Greenland_DC8.xls'),[],'post');
% params = read_param_xls(ct_filename_param('kuband_param_2011_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('kuband_param_2012_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('rds_param_1993_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('rds_param_1995_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('rds_param_1996_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('rds_param_1997_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('rds_param_1998_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('rds_param_1999_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('rds_param_2001_Greenland_P3.xls'),[],'post');
params = read_param_xls(ct_filename_param('rds_param_2002_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2009_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_DC8.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2009_Antarctica_DC8.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2010_Antarctica_DC8.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2011_Antarctica_DC8.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2012_Antarctica_DC8.xls'),[],'post');

source_dir = gRadar.out_path;
backup_dirs = {};
dirs_list = [source_dir backup_dirs];
support_dir = gRadar.support_path;
support_backup_dirs = {'',''};
support_dirs_list = [support_dir support_backup_dirs];

if any(strcmp(params(1).radar_name,{'icards','mcrds','mcords','mcords2','mcords3','mcords4','mcords5'}))
  supports = {'gps','frames','records'};
  outputs = {'CSARP_qlook','CSARP_standard','CSARP_mvdr','CSARP_layerData','CSARP_out'};
%   outputs = {'CSARP_qlook','CSARP_csarp-combined'};
  outputs_post_dir = '';
  images = {'maps','echo'};
  pdf_en = 1;
  csv_outputs = {'csv','csv_good','kml','kml_good'};
  csv_en = 1;
elseif strcmp(params(1).radar_name,'accum2')
  supports = {'gps','frames','records'};
  outputs = {'CSARP_qlook','CSARP_layerData'};
  outputs_post_dir = 'CSARP_post';
  images = {'maps','echo'};
  pdf_en = 1;
  csv_outputs = {'csv','csv_good','kml','kml_good'};
  csv_en = 0;
elseif any(strcmp(params(1).radar_name,{'kaband3','kuband3','snow3','kuband2','snow2','kuband','snow'}))
  supports = {'gps','frames','records'};
  outputs = {'CSARP_qlook'};
  outputs_post_dir = 'CSARP_post';
  images = {'maps','echo'};
  pdf_en = 0;
  csv_en = 0;
end
% gps_sources = {'ATM-final_20120701'}; % Leave empty/undefined to not check gps_sources
% processing_date_check = datenum(2012,09,01); % Leave empty/undefined to not check porcessing date
gps_sources = {}; % Leave empty/undefined to not check gps_sources
processing_date_check = []; % Leave empty/undefined to not check processing date
frm_types = {-1,0,-1,-1,-1};
delete_bad_files = true;
check_for_bad_files = true;
enable_all_without_do_not_process = true;

%% Automated Section

check_data_products;

return;
