% script run_check_data_products
%
% Script for running check_data_products.m.
%
% Author: John Paden, Logan Smith

%% User Settings
clear; clc
global gRadar;

params = read_param_xls(ct_filename_param('snow_param_2017_Greenland_P3.xls'));
params = ct_set_params(params,'cmd.generic',1);
params = ct_set_params(params,'cmd.generic',0,'cmd.notes','Do not process');
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20170510_01');

source_dir = gRadar.out_path;
backup_dirs = {};
dirs_list = [source_dir backup_dirs];
support_dir = gRadar.support_path;
support_backup_dirs = {'',''};
support_dirs_list = [support_dir support_backup_dirs];

[~,~,radar_name] = ct_output_dir(params(1).radar_name);
if any(strcmp(radar_name,{'hfrds','hfrds2','icards','mcrds','mcords','mcords2','mcords3','mcords4','mcords5'}))
  supports = {'gps','frames','records'};
%   outputs = {'CSARP_qlook','CSARP_standard','CSARP_mvdr','CSARP_layerData','CSARP_out'};
  outputs = {'CSARP_qlook','CSARP_standard'};
%   outputs = {'CSARP_qlook','CSARP_standard','CSARP_mvdr'};
%   outputs = {'CSARP_qlook','CSARP_csarp-combined'};
  outputs_post_dir = 'CSARP_post';
  images = {'maps','echo'};
  pdf_en = 1;
  csv_outputs = {'csv','csv_good','kml','kml_good'};
  csv_en = 1;
elseif strcmp(radar_name,'accum2')
  supports = {'gps','frames','records'};
  outputs = {'CSARP_qlook','CSARP_layerData'};
  outputs_post_dir = 'CSARP_post';
  images = {'maps','echo'};
  pdf_en = 1;
  csv_outputs = {'kml'};
  csv_en = 1;
elseif any(strcmp(radar_name,{'kaband3','kuband3','snow3','kuband2','snow2','kuband','snow'}))
  supports = {'gps','frames','records'};
  outputs = {'CSARP_qlook'};
  outputs_post_dir = 'CSARP_post';
  images = {'maps','echo'};
  pdf_en = 0;
  csv_outputs = {'kml'};
  csv_en = 1;
elseif any(strcmp(radar_name,{'snow8'}))
  % Run 1
  supports = {'gps','frames','records'};
  outputs = {'CSARP_qlook'};
  outputs_post_dir = 'CSARP_post';
  images = {'maps','echo'};
  pdf_en = 0;
  csv_outputs = {'kml'};
  csv_en = 1;
%   % Run 2
%   supports = {};
%   outputs = {};
%   outputs_post_dir = 'CSARP_post_uwb';
%   params = ct_set_params(params,'cmd.generic',0,'cmd.notes','2-6');
%   params = ct_set_params(params,'cmd.generic',0,'cmd.notes','2-8');
%   images = {'echo'};
%   pdf_en = 0;
%   csv_outputs = {'kml'};
%   csv_en = 0;
%   % Run 3
%   supports = {};
%   outputs = {};
%   outputs_post_dir = 'CSARP_post_kuband';
%   params = ct_set_params(params,'cmd.generic',0,'cmd.notes','2-6');
%   params = ct_set_params(params,'cmd.generic',0,'cmd.notes','2-8');
%   images = {'echo'};
%   pdf_en = 0;
%   csv_outputs = {'kml'};
%   csv_en = 0;
%   % Run 4
%   supports = {};
%   outputs = {};
%   outputs_post_dir = 'CSARP_post_deconv';
%   params = ct_set_params(params,'cmd.generic',0,'cmd.notes','2-6');
%   images = {'echo'};
%   pdf_en = 0;
%   csv_outputs = {'kml'};
%   csv_en = 0;
%   % Run 5
%   supports = {};
%   outputs = {'CSARP_deconv'};
%   params = ct_set_params(params,'cmd.generic',0,'cmd.notes','2-6');
%   outputs_post_dir = 'CSARP_post';
%   images = {};
%   pdf_en = 0;
%   csv_outputs = {'kml'};
%   csv_en = 0;
%   % Run 6
%   supports = {};
%   outputs = {'CSARP_qlook_kuband','CSARP_qlook_uwb'};
%   params = ct_set_params(params,'cmd.generic',0,'cmd.notes','2-6');
%   params = ct_set_params(params,'cmd.generic',0,'cmd.notes','2-8');
%   outputs_post_dir = 'CSARP_post';
%   images = {};
%   pdf_en = 0;
%   csv_outputs = {'kml'};
%   csv_en = 0;
end
gps_sources = {'atm-final_20170620'}; % Should be checked at least once before data posting
% gps_sources = {}; % Leave empty/undefined to not check gps_sources

% processing_date_check = datenum(2012,09,01); % Leave empty/undefined to not check processing date
processing_date_check = []; % Leave empty/undefined to not check processing date

frm_types = {-1,0,-1,-1,-1};
delete_bad_files = false;
check_echogram_type = true; % Should be checked at least once before data posting
check_for_bad_files = true;
enable_all_without_do_not_process = false;

%% Automated Section

check_data_products;

return;
