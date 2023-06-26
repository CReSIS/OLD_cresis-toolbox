% Scripts run_nsidc_delivery_script
%
% Runs nsidc_delivery_script
%
% Author: Yi Zhu, John Paden

%% User Settings
% User defined directory  
USER_SPECIFIED_DIRECTORY_BASE = '/cresis/snfs1/scratch/paden/nsidc/';

% Read rds spreadsheet
% params_fn = ct_filename_param('snow_param_2009_Greenland_P3.xls');
% params_fn = ct_filename_param('snow_param_2010_Greenland_DC8.xls');
% params_fn = ct_filename_param('snow_param_2010_Greenland_P3.xls');
% params_fn = ct_filename_param('snow_param_2011_Greenland_P3.xls');
% params_fn = ct_filename_param('snow_param_2012_Greenland_P3.xls');
params_fn = ct_filename_param('snow_param_2017_Antarctica_P3.xls');

% Post L1B data
L1B_cmd = true;
data_dir_L1 = fullfile('CSARP_post','qlook'); % Snow radar, Ku-band radar
% data_dir_L1 = fullfile('CSARP_post','standard'); % Accumulation radar
% data_dir_L1 = fullfile('CSARP_post','mvdr'); % RDS

% Post extra L1B data with the primary L1B data
% data_dir_L1_extra = {fullfile('CSARP_post','deconv'), 'deconv'}; % Snow radar 2-8 GHz with deconv file
data_dir_L1_extra = {fullfile('CSARP_post','deconv'), 'deconv';fullfile('CSARP_post','qlook_uwb'), 'uwb';fullfile('CSARP_post','qlook_kuband'), 'kuband'}; % Snow radar 2-18 GHz
% data_dir_L1_extra = {}; % All others

image_extra = {'deconv','uwb','kuband'};  % Snow radar 2-18 GHz
% image_extra = {};  % All others

% Enable creation of classification mask files
L1B_supplement_cmd = true; % Snow radar with deconv file
% L1B_supplement_cmd = false; % All others
L1B_supplement_name = 'supplement';
% L1B_supplement_name_extra = {'uwb'}; % Snow: If separate supplement files for different products, then e.g. {'uwb';'kuband'}
L1B_supplement_name_extra = {''}; % Snow: If no separate supplement files

% Post L2 data
% L2_cmd = true; % RDS layer data
L2_cmd = false; % Otherwise
data_dir_L2 = fullfile('CSARP_post','csv');

% Hardcoded local version ID (our local version number, check with NSIDC
% before changing if re-sending data)
premet_param_L2.version_id = '001';
premet_param_L1B.version_id = '002';  

% MCF version ID (only change if NSIDC requests change)
mcf_version_id_L1B = '002';
mcf_version_id_L2 = '001';

% Frame types to post (frames.proc_mode), this should not be changed
frm_types = {-1,0,-1,-1,-1};

clear('param_override');
param_override = [];
param_override.sched.type = 'no scheduler';
param_override.sched.cluster_size = inf;
param_override.sched.rerun_only = false;
param_override.sched.conforming = true;
param_override.sched.submit_arguments = '-l nodes=1:ppn=1,pmem=3450mb,walltime=480:00';
param_override.sched.stop_on_fail = false;
param_override.sched.test_mode = false;
param_override.sched.max_in_queue = 4*62;
global gRadar;
param_override = merge_structs(gRadar,param_override);

%% Automated Section

nsidc_delivery_script;

return;
