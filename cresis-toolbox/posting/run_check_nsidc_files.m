% script run_check_nsidc_files
%
% Script for checking if the number of  nsidc files generated matches with expected
% Authors: Jilu Li
%
% See also: check_nsidc_files
%% User Setup
% =====================================================================
param_override = [];
params = read_param_xls(ct_filename_param('snow_param_2018_Antarctica_DC8.xls'));

% Example to run a specific segment and frame by overriding parameter spreadsheet values
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20181027_02');
param_override.nsidc_file_folder = '/cresis/snfs1/scratch/jliwestc/nsidc/2018_Antarctica_DC8/snow/IRSNO1B_Files/Data_2018_AN';
param_override.nsidc_SP_file_folder = '/cresis/snfs1/scratch/jliwestc/nsidc/2018_Antarctica_DC8/snow/IRSNO1B_Files/Spatial_Premet_Files';

%% Automated Section
% =====================================================================
% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%% Process each of the segments
% =====================================================================
all_file_num = 0; % (not include Spatial and Premet files)
all_nc_file_num = 0;
all_supplement_file_num = 0;
all_SP_file_num = 0;
all_map_file_num = 0;
all_img_file_num = 0;
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  [total_num,nc_num,supplement_num,SP_num,map_num,img_num] = check_nsidc_files(param,param_override);
  all_file_num = all_file_num + total_num;
  all_nc_file_num = all_nc_file_num + nc_num;
  all_supplement_file_num = all_supplement_file_num + supplement_num;
  all_SP_file_num = all_SP_file_num + SP_num;
  all_map_file_num = all_map_file_num + map_num;
  all_img_file_num = all_img_file_num + img_num;
end
sprintf('number of files to deliver = %d\nnumber of nc files = %d\nnumber of supplement files = %d\nnumber of Spatial and Premet files = %d\nnumber of frames or map files = %d\nnumber of echogram images = %d\n',...
  all_file_num, all_nc_file_num,all_supplement_file_num,all_SP_file_num,all_map_file_num, all_img_file_num)

