%% data setup  (refer to run_tracker.m)
params             = read_param_xls(ct_filename_param('rds_param_2017_Greenland_P3_random_tests.xls'),'','post');
options.name       = 'CSARP_post/mvdr';
options.debug      = false;
geotiff_fn         = 'greenland/IceMask/GimpIceMask_90m_v1.1.tif';
geotiff_fn         = ct_filename_gis([],geotiff_fn);
geotiff2_fn        = '';

% Note: If Antartica 2D data is selected, please go to line 170 in
% Viterbi_2D_kernel.m change 'if 0' to 'if 1'

global gRadar;
physical_constants;
clear('param_override');
% Input checking
if ~exist('params','var')
  error('Use run_tracker: A struct array of parameters must be passed in\n');
end
if exist('param_override','var')
  param_override = merge_structs(gRadar, param_override);
else
  param_override = gRadar;
end


%% grid setup

% Viterbi_param.smooth_weight = [90 100 110 120 150 160 170 180 190];
% Viterbi_param.smooth_var = [0 50000 100000 1000000];
% Viterbi_param.repulsion = [0 50000 100000];
% Viterbi_param.ice_b_thr = [0 5 10 15 20];

Viterbi_param.smooth_weight = [120];
Viterbi_param.smooth_var = [1000000];
Viterbi_param.repulsion = [100000];
Viterbi_param.ice_b_thr = [20];