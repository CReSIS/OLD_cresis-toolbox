%   Interface to run MCMC_HMM_automated_frames.m with given parameters.
%
%   See also: MCMC_HMM_automated_frames.m
%
%   You may need to recompile RJ_MCMC.cpp (MCMC) and stereo.cpp (HMM)
%     example: mex -largeArrayDims RJ_MCMC.cpp
%
%   Author: Victor Berger

%{

Parameters:
Alg: 'MCMC' or 'HMM' -- which layer tracking algorithm to use
Debug: LOGICAL -- display extra info
Display: LOGICAL -- display image on screen
PNG_Output: LOGICAL -- whether or not to save a .png file
        (this file is saved to img_output_dir)
Save_Layer_Data: LOGICAL -- whether or not to save the new layerData
        (this file is saved to output_layer_dir)

%}

MCMC_HMM_param.alg = 'MCMC';
MCMC_HMM_param.debug = false;
MCMC_HMM_param.display = false;
MCMC_HMM_param.png_output = true;
MCMC_HMM_param.save_layer_data = true;

% Directories:
MCMC_HMM_param.param_dir = '/users/victor/scripts/ct_params/rds_param_2016_Antarctica_DC8.xls';
MCMC_HMM_param.data_dir = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Antarctica_DC8/CSARP_standard';
MCMC_HMM_param.img_output_dir = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Antarctica_DC8/CSARP_standard_MCMC';
MCMC_HMM_param.original_layer_dir = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Antarctica_DC8/CSARP_post/CSARP_layerData';
MCMC_HMM_param.output_layer_dir = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Antarctica_DC8/CSARP_layerData_MCMC';

%-------------------------------%


if (~exist(MCMC_HMM_param.param_dir, 'file') || ~exist(MCMC_HMM_param.data_dir, 'dir') || ...
        (MCMC_HMM_param.png_output == true && ~exist(MCMC_HMM_param.img_output_dir, 'dir')) || ...
        ~exist(MCMC_HMM_param.original_layer_dir, 'dir') || ~exist(MCMC_HMM_param.output_layer_dir, 'dir'))
    error('Error: one or more directories were not found. Check the path and make sure they exist.')
end

if (~islogical(MCMC_HMM_param.debug) || ~islogical(MCMC_HMM_param.display) || ...
        ~islogical(MCMC_HMM_param.png_output) || ~islogical(MCMC_HMM_param.save_layer_data) || ...
        (~strcmp(MCMC_HMM_param.alg, 'MCMC') && ~strcmp(MCMC_HMM_param.alg, 'HMM')))
    error('Error with parameters. All should be logical (true/false), except algorithm (MCMC/HMM).');
end

% Actual function call
MCMC_HMM_automated_frames(MCMC_HMM_param);                          

