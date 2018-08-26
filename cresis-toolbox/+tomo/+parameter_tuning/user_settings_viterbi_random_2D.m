%% setup  (refer to run_tracker.m)
params             = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_DC8_test1.xls'),'','post');
options.name       = 'CSARP_post/mvdr';
options.debug      = false;

geotiff_fn         = 'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_icemask_grounded_and_shelves.tif';
geotiff2_fn        = 'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_rockmask.tif';
% geotiff_fn         = 'greenland/IceMask/GimpIceMask_90m_v1.1.tif';
geotiff_fn         = ct_filename_gis([],geotiff_fn);
geotiff2_fn        = ct_filename_gis([],geotiff2_fn);


% Note: If Antartica 2D data is selected, please go to line 174 in
% cluster_Viterbi_2D_kernel.m change 'if 0' to 'if 1'

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


%% parameters

num_trials = 1;               % number of random vectors to generate

Viterbi_param.smooth_weight_min = 80;       % setup the lower bound for smoothness weight for random number generation
Viterbi_param.smooth_weight_max = 105;      % setup the upper bound for smoothness weight for random number generation

Viterbi_param.smooth_var_min = 140000;      % setup the lower bound for smoothness variance for random number generation
Viterbi_param.smooth_var_max = 160000;      % setup the upper bound for smoothness variance for random number generation

Viterbi_param.repulsion_min = 0;            % setup the lower bound for repulsion for random number generation
Viterbi_param.repulsion_max = 200000;       % setup the upper bound for repulsion for random number generation

Viterbi_param.ice_b_thr_min = 10;           % setup the lower bound for ice distance threshold for random number generation
Viterbi_param.ice_b_thr_max = 20;           % setup the upper bound for ice distance threshold for random number generation