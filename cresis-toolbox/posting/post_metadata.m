% script post_metadata
%
% Posts metadata files to CSARP_post:
% 1. frames directory is copied
% 2. param spreadsheet is copied
% 3. CSARP_noise support files are moved (temporary files are deleted)
%
% Outputs:
%   gRadar.out_path/param.radar_name/param.season_name/CSARP_post/csv/
%     Data_YYYYMMDD_SS.csv
%     Browse_Data_YYYYMMDD_SS.kml
%
% Author: John Paden

%% User Settings

% params = read_param_xls(ct_filename_param('snow_param_2009_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_DC8.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenedland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2009_Antarctica_DC8.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2010_Antarctica_DC8.xls'),[],'post');
% params = read_param_xls(ct_filename_param('snow_param_2011_Antarctica_DC8.xls'),[],'post');
params = read_param_xls(ct_filename_param('snow_param_2012_Antarctica_DC8.xls'),[],'post');

copy_frames_files = true;
copy_param_spreadsheet = true;
move_analysis_files = true;

%% Automated Section

fprintf('=============================================================\n');
fprintf('post_metadata (%s)\n\n', datestr(now));
fprintf('=============================================================\n');

% =====================================================================
% Create param structure array
% =====================================================================
tic;
global gRadar;

clear('param_override');

% Input checking
if ~exist('params','var')
  error('Use run_master: A struct array of parameters must be passed in\n');
end
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%% Process each of the segments
% =====================================================================
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
 
  if copy_frames_files
    frames_dir = fullfile(ct_filename_out(param,'post','',true),'frames');
    if ~exist(frames_dir,'dir')
      mkdir(frames_dir)
    end
    frames_fn = ct_filename_support(param,'','frames');
    fprintf('Copying %s\n  %s\n', frames_fn, frames_dir);
    copyfile(frames_fn,frames_dir);
  end
  
  if copy_param_spreadsheet
    output_dir = ct_output_dir(param.radar_name);
    params_fn = ct_filename_param(sprintf('%s_param_%s.xls',output_dir,param.season_name));
    post_dir = ct_filename_out(param,'post','',true);
    if ~exist(post_dir,'dir')
      mkdir(post_dir)
    end
    fprintf('Copying %s\n  %s\n', params_fn, post_dir);
    copyfile(params_fn,post_dir);
    copy_param_spreadsheet = false;
  end
  
  if move_analysis_files
    analysis_in_dir = ct_filename_out(param,'noise','',true);
    analysis_dir = ct_filename_out(param,'CSARP_post/noise','',true);
    if ~exist(analysis_dir,'dir')
      mkdir(analysis_dir)
    end
    analysis_files{1} = {'coh_noise_%s.mat',param.day_seg};
    analysis_files{2} = {'coh_noise_simp_%s.nc',param.day_seg};
    analysis_files{3} = {'specular_%s.mat',param.day_seg};
    analysis_files{4} = {'deconv_%s.mat',param.day_seg};
    
    for analysis_file_idx = 1:length(analysis_files)
      analysis_fn = fullfile(analysis_in_dir,sprintf(analysis_files{analysis_file_idx}{:}));
      if exist(analysis_fn)
        fprintf('Moving %s\n  %s\n', analysis_fn, analysis_dir);
        movefile(analysis_fn,analysis_dir);
      else
        warning('File does not exist %s', analysis_fn);
      end
    end
  end
end

return;
