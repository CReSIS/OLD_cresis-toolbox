% script run_update_mult_factor
%
% Runs update_mult_factor.m
%
% Author: John Paden

%% User Settings
% ----------------------------------------------------------------------
% params = read_param_xls(ct_filename_param('kuband_param_2014_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('rds_param_2015_Greenland_Polar6.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2009_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_DC8.xls'));
params = read_param_xls(ct_filename_param('snow_param_2017_Greenland_P3.xls'),'20170311_02','analysis');
% params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'));
data_types = []; idx = 0;
% idx = idx+1;
% data_types(idx).type = 'coh_noise';
% data_types(idx).dir = 'analysis';
% data_types(idx).imgs = [1];
% data_types(idx).wf_adcs = [];
idx = idx+1;
data_types(idx).type = 'echogram';
data_types(idx).dir = 'qlook';
data_types(idx).imgs = [0];
idx = idx+1;
data_types(idx).type = 'echogram';
data_types(idx).dir = 'qlook_snow';
data_types(idx).imgs = [0];
idx = idx+1;
data_types(idx).type = 'echogram';
data_types(idx).dir = 'qlook_kuband';
data_types(idx).imgs = [0];

params = ct_set_params(params,'cmd.generic',1);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20170310_01');


%% Automated Section
% ----------------------------------------------------------------------

update_mult_factor;

return;
