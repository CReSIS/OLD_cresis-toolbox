% script compress_echogram
%
% Compresses echogram files (CSARP_*). This involves three steps:
% 1. Elevation compensation
% 2. Fast-time truncation according to param spreadsheet posting settings
% 3. Quantization of 32 bit float to user specified type (e.g. uint8)
%
% Inputs are full resolution qlook files and outputs are saved in posting
% directory.  The parameter file's command worksheet is used to determine
% which segments and frames are compressed (using the generic field).
% The posting worksheet is also used.
%
% Settings for NSIDC:
% compress_type = '';
% truncate_data = true;
% guard_depth_rng = [0 0];
% min_depth_rng = [-8 5];
% elev_comp = true;
% param_post_echo_depth_override_sea = '[min(Surface_Depth)-3 max(Surface_Depth)+5]'; % Sea Ice Segments Only
% param_post_echo_depth_override_land = '[min(Surface_Depth)-10 max(Surface_Depth)+80]'; % Land Ice Segments Only
% frm_types = {-1,0,-1,-1,-1}; % Which types of frames to compress (usual is {-1,0,-1,-1,-1})
%
% Author: John Paden

% ====================================================================
% User Settings
% ====================================================================
% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2009_Greenland_P3/CSARP_qlook/';
% param_fn = ct_filename_param('snow_param_2009_Greenland_P3.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2009_Greenland_P3/CSARP_post/CSARP_qlook/';
%
% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2009_Greenland_P3/CSARP_deconv/';
% param_fn = ct_filename_param('snow_param_2009_Greenland_P3.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2009_Greenland_P3/CSARP_post/CSARP_deconv/';
% 
% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Greenland_DC8/CSARP_qlook/';
% param_fn = ct_filename_param('snow_param_2010_Greenland_DC8.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Greenland_DC8/CSARP_post/CSARP_qlook/';
% 
% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Greenland_DC8/CSARP_deconv/';
% param_fn = ct_filename_param('snow_param_2010_Greenland_DC8.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Greenland_DC8/CSARP_post/CSARP_deconv/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Greenland_P3/CSARP_qlook/';
% param_fn = ct_filename_param('snow_param_2010_Greenland_P3.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Greenland_P3/CSARP_post/CSARP_qlook/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Greenland_P3/CSARP_deconv/';
% param_fn = ct_filename_param('snow_param_2010_Greenland_P3.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Greenland_P3/CSARP_post/CSARP_deconv/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2011_Greenland_P3/CSARP_qlook/';
% param_fn = ct_filename_param('snow_param_2011_Greenland_P3.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2011_Greenland_P3/CSARP_post/CSARP_qlook/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2011_Greenland_P3/CSARP_deconv/';
% param_fn = ct_filename_param('snow_param_2011_Greenland_P3.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2011_Greenland_P3/CSARP_post/CSARP_deconv/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Greenland_P3/CSARP_qlook/';
% param_fn = ct_filename_param('snow_param_2012_Greenland_P3.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Greenland_P3/CSARP_post/CSARP_qlook/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Greenland_P3/CSARP_deconv/';
% param_fn = ct_filename_param('snow_param_2012_Greenland_P3.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Greenland_P3/CSARP_post/CSARP_deconv/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2009_Antarctica_DC8/CSARP_qlook/';
% param_fn = ct_filename_param('snow_param_2009_Antarctica_DC8.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2009_Antarctica_DC8/CSARP_post/CSARP_qlook/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2009_Antarctica_DC8/CSARP_deconv/';
% param_fn = ct_filename_param('snow_param_2009_Antarctica_DC8.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2009_Antarctica_DC8/CSARP_post/CSARP_deconv/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Antarctica_DC8/CSARP_qlook/';
% param_fn = ct_filename_param('snow_param_2010_Antarctica_DC8.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Antarctica_DC8/CSARP_post/CSARP_qlook/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Antarctica_DC8/CSARP_deconv/';
% param_fn = ct_filename_param('snow_param_2010_Antarctica_DC8.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Antarctica_DC8/CSARP_post/CSARP_deconv/';

echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2011_Antarctica_DC8/CSARP_qlook/';
param_fn = ct_filename_param('snow_param_2011_Antarctica_DC8.xls');
out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2011_Antarctica_DC8/CSARP_post/CSARP_qlook/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2011_Antarctica_DC8/CSARP_deconv/';
% param_fn = ct_filename_param('snow_param_2011_Antarctica_DC8.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2011_Antarctica_DC8/CSARP_post/CSARP_deconv/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Antarctica_DC8/CSARP_qlook/';
% param_fn = ct_filename_param('snow_param_2012_Antarctica_DC8.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Antarctica_DC8/CSARP_post/CSARP_qlook/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Antarctica_DC8/CSARP_deconv/';
% param_fn = ct_filename_param('snow_param_2012_Antarctica_DC8.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Antarctica_DC8/CSARP_post/CSARP_deconv/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/kuband/2012_Antarctica_DC8/CSARP_qlook/';
% param_fn = ct_filename_param('kuband_param_2012_Antarctica_DC8.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/kuband/2012_Antarctica_DC8/CSARP_post/CSARP_qlook/';

% echogram_dir = '/cresis/snfs1/dataproducts/ct_data/kuband/2009_Antarctica_DC8/CSARP_qlook/';
% param_fn = ct_filename_param('kuband_param_2009_Antarctica_DC8.xls');
% out_dir = '/cresis/snfs1/dataproducts/ct_data/kuband/2009_Antarctica_DC8/CSARP_post/CSARP_qlook/';

% Generally speaking do not change these parameters
compress_type = '';
truncate_data = true;
guard_depth_rng = [0 0];
min_depth_rng = [-8 5];
elev_comp = true;
param_post_echo_depth_override_sea = '[min(Surface_Depth)-3 max(Surface_Depth)+5]'; % Sea Ice Segments Only
param_post_echo_depth_override_land = '[min(Surface_Depth)-10 max(Surface_Depth)+80]'; % Land Ice Segments Only
frm_types = {-1,0,-1,-1,-1}; % Which types of frames to compress (usual is {-1,0,-1,-1,-1})

% ====================================================================
% Automated Section
% ====================================================================

compress_echogram;

return;
