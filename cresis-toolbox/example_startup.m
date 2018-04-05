% script startup.m
%
% This is a template startup.m file.
%
% STEPS TO SETUP MATLAB ENVIRONMENT ARE ON WIKI:
%  https://wiki.cresis.ku.edu/
%
% Step 1
%  PC:
%  Rename this file to startup.m and place in:
%    C:\Users\USERNAME\Documents\MATLAB
%
%  LINUX-KU:
%  Rename this file to startup.m and place in /users/USERNAME/.matlab/
%  LINUX-IU:
%  Rename this file to startup.m and place in /N/u/USERNAME/Quarry/.matlab/
%
% Step 2
%  Modify one of the profiles below and set cur_profile
%  - Probably just need to update personal_path and tmp_file_path
%
% Step 3
%  Restart matlab making suring the this startup script runs.
%
% Authors: John Paden

if ~(~ismcc && isdeployed)
  
  %% ===> SELECT CURRENT PROFILE HERE <===
  % cur_profile: Set the current profile to load (match to pidx below)
  % KU Profile Linux (PROFILE 1)
  % IU Profile Linux (PROFILE 2)
  % Field Profile Linux (PROFILE 3)
  % KU Mobile Profile Windows (PROFILE 5)
  % KU Desktop Profile Windows (PROFILE 6)
  % AWI Profile Windows (PROFILE 7)
  % AWI Profile Linux (PROFILE 8)
  % AWI Profile Ollie (PROFILE 9)
  if ispc
    cur_profile = 6;
  else
    cur_profile = 1;
  end
  
  fprintf('Startup Script Running\n');
  
  format short; format compact;
  
  %% KU Profile Linux (PROFILE 1)
  % ----------------------------------------------------------------------
  pidx = 1; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = '/users/paden/scripts/matlab/';
  profile(pidx).ct_path                   = '/users/paden/scripts/cresis-toolbox/cresis-toolbox/';
  profile(pidx).param_path                = '/users/paden/scripts/ct_params/';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = '/cresis/snfs1/scratch/paden/mdce_tmp/';
  
  profile(pidx).data_path                 = '/cresis/snfs1/data/';
  profile(pidx).data_support_path         = '/cresis/snfs1/dataproducts/metadata/';
  profile(pidx).support_path              = '/cresis/snfs1/dataproducts/csarp_support/';
  profile(pidx).out_path                  = '/cresis/snfs1/dataproducts/ct_data/';
  profile(pidx).gis_path                  = '/cresis/snfs1/dataproducts/GIS_data/';
  profile(pidx).ct_tmp_file_path          = fullfile(profile(pidx).out_path,'ct_tmp');
  
  profile(pidx).cluster.data_location       = fullfile(profile(pidx).tmp_file_path,'cluster-temp');
  
  profile(pidx).cluster.type                  = 'torque';
  %profile(pidx).cluster.type                  = 'matlab';
  %profile(pidx).cluster.type                  = 'debug';
  profile(pidx).cluster.max_jobs_active       = 512;
  profile(pidx).cluster.max_time_per_job      = 86400;
  profile(pidx).cluster.desired_time_per_job  = 0;
  profile(pidx).cluster.max_retries           = 2;
  profile(pidx).cluster.submit_pause          = 0.2;
  profile(pidx).cluster.stat_pause            = 2;
  profile(pidx).cluster.file_check_pause      = 4;
  
  %% IU Profile Linux (PROFILE 2)
  % ----------------------------------------------------------------------
  pidx = 2; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = '/N/u/jpaden/Karst/scripts/matlab/';
  profile(pidx).ct_path                   = '/N/u/jpaden/Karst/scripts/cresis-toolbox/cresis-toolbox/';
  profile(pidx).param_path                = '/N/u/jpaden/Karst/scripts/ct_params/';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = '/N/dc2/scratch/jpaden/mdce_tmp/'; % scratch may be on dcwan or dc2
  
  profile(pidx).data_path                 = '/N/dcwan/projects/cresis/';
  profile(pidx).data_support_path         = '/N/dcwan/projects/cresis/metadata/';
  profile(pidx).support_path              = '/N/dcwan/projects/cresis/csarp_support/';
  profile(pidx).out_path                  = '/N/dcwan/projects/cresis/output/';
  profile(pidx).gis_path                  = '/N/dcwan/projects/cresis/GIS_data';
  profile(pidx).ct_tmp_file_path          = fullfile(profile(pidx).out_path,'ct_tmp');

  profile(pidx).cluster.data_location       = fullfile(profile(pidx).tmp_file_path,'cluster-temp');
  
  profile(pidx).cluster.type                  = 'torque';
  %profile(pidx).cluster.type                  = 'matlab';
  %profile(pidx).cluster.type                  = 'debug';
  profile(pidx).cluster.max_jobs_active       = 128;
  profile(pidx).cluster.max_time_per_job      = 4*86400;
  profile(pidx).cluster.desired_time_per_job  = 8*3600;
  profile(pidx).cluster.max_retries           = 2;
  profile(pidx).cluster.submit_pause          = 0.5;
  profile(pidx).cluster.stat_pause            = 2;
  profile(pidx).cluster.file_check_pause      = 4;
  
  profile(pidx).cluster.qsub_submit_arguments = '-m n -l nodes=1:ppn=%p:dcwan:dc2,pmem=%m,walltime=%t';

 
  %% Field Profile Linux (PROFILE 3)
  % ----------------------------------------------------------------------
  pidx = 3; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = '/scratch/scripts/matlab/';
  profile(pidx).ct_path                   = '/scratch/scripts/cresis-toolbox/cresis-toolbox/';
  profile(pidx).param_path                = '/scratch/scripts/ct_params/';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = '/scratch/tmp/';
  
  profile(pidx).data_path                 = '/scratch/';
  profile(pidx).data_support_path         = '/scratch/metadata/';
  profile(pidx).support_path              = '/scratch/csarp_support/';
  profile(pidx).out_path                  = '/scratch/';
  profile(pidx).gis_path                  = '/scratch/GIS_data';
  profile(pidx).ct_tmp_file_path          = fullfile(profile(pidx).out_path,'ct_tmp');
  
  profile(pidx).cluster.data_location       = fullfile(profile(pidx).tmp_file_path,'cluster-temp');
  
  profile(pidx).cluster.type                  = 'torque';
  %profile(pidx).cluster.type                  = 'matlab';
  %profile(pidx).cluster.type                  = 'debug';
  profile(pidx).cluster.max_jobs_active       = 512;
  profile(pidx).cluster.max_time_per_job      = 86400;
  profile(pidx).cluster.desired_time_per_job  = 0;
  profile(pidx).cluster.max_retries           = 2;
  profile(pidx).cluster.submit_pause          = 0.2;
  profile(pidx).cluster.stat_pause            = 2;
  profile(pidx).cluster.file_check_pause      = 4;
  profile(pidx).cluster.mem_to_ppn            = 0.9 * 131754468000 / 46;
  
  %% KU Mobile Profile Windows (PROFILE 5)
  % ----------------------------------------------------------------------
  pidx = 5; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = 'C:\Users\paden\Documents\scripts\matlab\';
  profile(pidx).ct_path                   = 'C:\Users\paden\Documents\scripts\cresis-toolbox\cresis-toolbox\';
  profile(pidx).param_path                = 'C:\Users\paden\Documents\scripts\ct_params\';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = 'C:\tmp\mdce_tmp\';
  
  profile(pidx).data_path                 = 'D:\';
  profile(pidx).data_support_path         = 'C:\metadata\';
  profile(pidx).support_path              = 'C:\csarp_support\';
  profile(pidx).out_path                  = 'D:\output\';
  profile(pidx).gis_path                  = 'C:\GIS_data\';
  profile(pidx).ct_tmp_file_path          = fullfile(profile(pidx).out_path,'ct_tmp');
  
  profile(pidx).cluster.data_location       = fullfile(profile(pidx).tmp_file_path,'cluster-temp');
  
  profile(pidx).cluster.type                  = 'matlab';
  %profile(pidx).cluster.type                  = 'debug';
  %profile(pidx).cluster.type                  = 'none';
  profile(pidx).cluster.max_jobs_active       = 4;
  profile(pidx).cluster.max_time_per_job      = 86400;
  profile(pidx).cluster.desired_time_per_job  = 0;
  profile(pidx).cluster.max_retries           = 1;
  profile(pidx).cluster.submit_pause          = 0;
  profile(pidx).cluster.stat_pause            = 1;
  
  %% KU Desktop Profile Windows (PROFILE 6)
  % ----------------------------------------------------------------------
  % H:\ --> \\emperor.cresis.ku.edu\paden\
  % S:\ --> \\cfs1.cresis.ku.edu\\
  % P:\ --> \\titan.cresis.ku.edu\projects\
  % V:\ --> \\titan.cresis.ku.edu\data3\
  % W:\ --> \\titan.cresis.ku.edu\data1\
  % X:\ --> \\titan.cresis.ku.edu\data2\
  % Y:\ --> \\titan.cresis.ku.edu\scratch1\
  % Z:\ --> \\titan.cresis.ku.edu\scratch2\
  pidx = 6; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = 'C:\Users\paden\Documents\scripts\matlab\';
  profile(pidx).ct_path                   = 'C:\Users\paden\Documents\scripts\cresis-toolbox\cresis-toolbox\';
  profile(pidx).param_path                = 'C:\Users\paden\Documents\scripts\ct_params\';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = 'Y:/paden/mdce_tmp/';
  
  profile(pidx).data_path                 = 'V:/';
  profile(pidx).data_support_path         = 'X:/metadata/';
  profile(pidx).support_path              = 'X:/csarp_support/';
  profile(pidx).out_path                  = 'X:/ct_data/';
  profile(pidx).gis_path                  = 'X:/GIS_data/';
  profile(pidx).ct_tmp_file_path          = fullfile(profile(pidx).out_path,'ct_tmp');
  
  profile(pidx).cluster.data_location       = fullfile(profile(pidx).tmp_file_path,'cluster-temp');
  
  profile(pidx).cluster.type                  = 'matlab';
  %profile(pidx).cluster.type                  = 'debug';
  %profile(pidx).cluster.type                  = 'none';
  profile(pidx).cluster.max_jobs_active       = 4;
  profile(pidx).cluster.max_time_per_job      = 86400;
  profile(pidx).cluster.desired_time_per_job  = 0;
  profile(pidx).cluster.max_retries           = 1;
  profile(pidx).cluster.submit_pause          = 0;
  profile(pidx).cluster.stat_pause            = 1;
  
  %% AWI Profile Windows (PROFILE 7)
  % ----------------------------------------------------------------------
  pidx = 7; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = 'C:\tmp\scripts\matlab\';
  profile(pidx).ct_path                   = 'C:\tmp\scripts\cresis-toolbox\cresis-toolbox\';
  profile(pidx).param_path                = 'C:\tmp\scripts\ct_params\';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = 'F:\mdce_tmp\';
  
  profile(pidx).data_path                 = 'D:\';
  profile(pidx).data_support_path         = 'F:\metadata\';
  profile(pidx).support_path              = 'F:\csarp_support\';
  profile(pidx).out_path                  = 'F:\';
  profile(pidx).gis_path                  = 'C:\tmp\GIS_data\';
  profile(pidx).ct_tmp_file_path          = fullfile(profile(pidx).out_path,'ct_tmp');
  
  profile(pidx).cluster.data_location       = fullfile(profile(pidx).tmp_file_path,'cluster-temp');

  profile(pidx).cluster.type                  = 'matlab';
  %profile(pidx).cluster.type                  = 'debug';
  %profile(pidx).cluster.type                  = 'none';
  profile(pidx).cluster.max_jobs_active       = 4;
  profile(pidx).cluster.max_time_per_job      = 86400;
  profile(pidx).cluster.desired_time_per_job  = 0;
  profile(pidx).cluster.max_retries           = 1;
  profile(pidx).cluster.submit_pause          = 0;
  profile(pidx).cluster.stat_pause            = 1;
  
  %% AWI Profile Linux (PROFILE 8)
  % ----------------------------------------------------------------------
  pidx = 8; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = '/home/administrator/scripts/matlab/';
  profile(pidx).ct_path                   = '/home/administrator/scripts/cresis-toolbox/cresis-toolbox/';
  profile(pidx).param_path                = '/home/administrator/scripts/ct_params/';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = '/home/administrator/Scratch/mdce_tmp/';
  
  profile(pidx).data_path                 = '/mnt/AWI_SSD0/';
  profile(pidx).data_support_path         = '/home/administrator/Scratch/metadata/';
  profile(pidx).support_path              = '/home/administrator/Scratch/csarp_support/';
  profile(pidx).out_path                  = '/home/administrator/Scratch/';
  profile(pidx).gis_path                  = '/home/administrator/GIS_data/';
  profile(pidx).ct_tmp_file_path          = fullfile(profile(pidx).out_path,'ct_tmp');
  
  profile(pidx).cluster.data_location       = fullfile(profile(pidx).tmp_file_path,'cluster-temp');
  
  profile(pidx).cluster.type                  = 'torque';
  %profile(pidx).cluster.type                  = 'matlab';
  %profile(pidx).cluster.type                  = 'debug';
  profile(pidx).cluster.max_jobs_active       = 512;
  profile(pidx).cluster.max_time_per_job      = 86400;
  profile(pidx).cluster.desired_time_per_job  = 0;
  profile(pidx).cluster.max_retries           = 2;
  profile(pidx).cluster.submit_pause          = 0.2;
  profile(pidx).cluster.stat_pause            = 2;
  profile(pidx).cluster.file_check_pause      = 4;

  %% AWI Profile Ollie (PROFILE 9)
  % ----------------------------------------------------------------------
  pidx = 9; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = '/home/ollie/jpaden/scripts/matlab/';
  profile(pidx).ct_path                   = '/home/ollie/jpaden/scripts/cresis-toolbox/cresis-toolbox/';
  profile(pidx).param_path                = '/home/ollie/jpaden/scripts/ct_params/';
  profile(pidx).slurm_jobs_path           = '/home/ollie/tbinder/jobs';

  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = '/work/ollie/jpaden/Scratch/mdce_tmp/';

  profile(pidx).data_path                 = '/work/ollie/jpaden/Data/';
  profile(pidx).data_support_path         = '/work/ollie/jpaden/Scratch/metadata/';
  profile(pidx).support_path              = '/work/ollie/jpaden/Scratch/csarp_support/';
  profile(pidx).out_path                  = '/work/ollie/jpaden/Scratch/';
  profile(pidx).gis_path                  = '/work/ollie/jpaden/GIS_data/';
  profile(pidx).ct_tmp_file_path          = fullfile(profile(pidx).out_path,'ct_tmp');
 
  profile(pidx).cluster.data_location       = fullfile(profile(pidx).tmp_file_path,'cluster-temp');

  %profile(pidx).cluster.type                    = 'matlab';
  profile(pidx).cluster.type                    = 'slurm';
  %profile(pidx).cluster.type                    = 'ollie';
  %profile(pidx).cluster.type                    = 'debug';
  profile(pidx).cluster.max_jobs_active         = 64;
  profile(pidx).cluster.max_time_per_job        = 86400;
  profile(pidx).cluster.desired_time_per_job    = 2*3600;
  profile(pidx).cluster.max_retries             = 2;
  profile(pidx).cluster.submit_pause            = 0.2;
  profile(pidx).cluster.stat_pause              = 2;
  profile(pidx).cluster.file_check_pause        = 4;

  profile(pidx).cluster.mcc                     = 'eval';
  
  %% Startup code (Automated Section)
  % =====================================================================
  
  fprintf('  Resetting path\n');
  path(pathdef);
  AdditionalPaths = {};
  
  if ~exist(profile(cur_profile).ct_path,'dir')
    fprintf('Cresis toolbox not found: %s\n', profile(cur_profile).ct_path);
  else
    fprintf('  Adding cresis path: %s\n',profile(cur_profile).ct_path);
    % Add get_filenames to the path
    addpath(fullfile(profile(cur_profile).ct_path));
    addpath(fullfile(profile(cur_profile).ct_path,'utility'));
    % Add each directory which is not an svn support directory to the path
    fns = get_filenames(profile(cur_profile).ct_path,'','','',struct('type','d','recursive',1));
    if ispc
      addpath(profile(cur_profile).ct_path);
      AdditionalPaths{end+1} = profile(cur_profile).ct_path;
    end
    for fn_idx = 1:length(fns)
      [fn_dir fn_name] = fileparts(fns{fn_idx});
      if ~isempty(fn_name) && fn_name(1) ~= '@' && fn_name(1) ~= '+' ...
          && isempty(strfind(fns{fn_idx},'.svn')) && isempty(strfind(fns{fn_idx},'.git'))
        % Ignore .svn directories and Matlab class and package directories
        addpath(fns{fn_idx});
        AdditionalPaths{end+1} = fns{fn_idx};
      end
    end
  end
  
  if ~exist(profile(cur_profile).personal_path,'dir')
    fprintf('Personal toolbox not found: %s\n', profile(cur_profile).personal_path);
  else
    % Add personal path after cresis-toolbox so that it overrides cresis-toolbox
    fprintf('  Adding personal path: %s\n',profile(cur_profile).personal_path);
    fns = get_filenames(profile(cur_profile).personal_path,'','','',struct('type','d','recursive',1));
    addpath(profile(cur_profile).personal_path);
    AdditionalPaths{end+1} = profile(cur_profile).personal_path;
    for fn_idx = 1:length(fns)
      [fn_dir fn_name] = fileparts(fns{fn_idx});
      if ~isempty(fn_name) && fn_name(1) ~= '@' && fn_name(1) ~= '+' ...
          && isempty(strfind(fns{fn_idx},'.svn')) && isempty(strfind(fns{fn_idx},'.git'))
        % Ignore .svn directories and Matlab class and package directories
        addpath(fns{fn_idx});
        AdditionalPaths{end+1} = fns{fn_idx};
      end
    end
  end
  
  fprintf('  Setting global preferences in global variable gRadar\n');
  global gRadar;
  
  % .cluster: structure with default cluster scheduler arguments
  gRadar.cluster = profile(cur_profile).cluster;
  % .sched.hidden_depends_funs = cell vector of strings containing
  %   all the dependent functions that should be included when the
  %   torque_compile function is called. This list should only include
  %   functions which are passed in (e.g. "@hanning" from param spreadsheet)
  %   since all directly dependent functions will automatically be included.
  %   Second cell entry is date check level in torque_compile:
  %     0: do not check for this file (system commands)
  %     1: check this file and its dependencies (function pointers)
  %     2: check this file and its dependencies only if "fun" is not
  %        specified in call to torque_compile.m (all functions called
  %        by torque_compile)
  gRadar.cluster.hidden_depend_funs = {};
  gRadar.cluster.hidden_depend_funs{end+1} = {'tomo_collate_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'create_records_accum2_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'create_records_acords_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'create_records_mcords_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'basic_load_fmcw.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'basic_load_fmcw2.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'rx_chan_equal_sar_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'rx_chan_equal_raw_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'coh_noise_tracker_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'coh_noise_tracker_combine_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'radiometric_calibration_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'get_heights_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'get_heights_combine_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'csarp_sar_coord_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'csarp_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'combine_wf_chan_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'combine_wf_chan_combine_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'nsidc_delivery_script_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'arena_packet_strip_task.m' 2};
  gRadar.cluster.hidden_depend_funs{end+1} = {'hanning.m' 0};
  gRadar.cluster.hidden_depend_funs{end+1} = {'hamming.m' 0};
  gRadar.cluster.hidden_depend_funs{end+1} = {'blackman.m' 0};
  gRadar.cluster.hidden_depend_funs{end+1} = {'tukeywin.m' 0};
  gRadar.cluster.hidden_depend_funs{end+1} = {'tukeywin_trim.m' 1};
  gRadar.cluster.hidden_depend_funs{end+1} = {'chebwin.m' 0};
  gRadar.cluster.hidden_depend_funs{end+1} = {'kaiser.m' 0};
  gRadar.cluster.hidden_depend_funs{end+1} = {'boxcar.m' 0};
  gRadar.cluster.hidden_depend_funs{end+1} = {'butter.m' 0};
  gRadar.cluster.hidden_depend_funs{end+1} = {'array_proc_sv.m' 1};
  gRadar.cluster.hidden_depend_funs{end+1} = {'lever_arm.m' 1};
  gRadar.cluster.hidden_depend_funs{end+1} = {'doa_nonlcon.m' 1};

  % .path = used by the scheduler to build the file dependency path
  %   typically this is the same at .ct_path
  gRadar.path = profile(cur_profile).code_path;
  % .path_override = used by the scheduler to build the file dependency
  %   path, however only files that also exist in the .path directory
  %   will be overwritten so that if a file only exists in .path_override
  %   it will not be included in the file dependency list
  gRadar.path_override = profile(cur_profile).code_path_override;
  % .param_path = parameter spreadsheet folder
  gRadar.param_path = profile(cur_profile).param_path;
  % .tmp_path = this is where personal temporary files will be stored (e.g.
  %   the picker should store files here)
  gRadar.tmp_path = profile(cur_profile).tmp_file_path;
  % .tmp_path = this is where global temporary files will be stored (e.g.
  %   create records, headers, etc. files should be stored here so that
  %   the process can be resumed by anyone without having to redo all
  %   the work or track down other people's temporary files)
  gRadar.ct_tmp_path = profile(cur_profile).ct_tmp_file_path;
  
  % .data_path = where the data is stored, the records and vectors files
  %   contain relative paths to the data from this base path
  gRadar.data_path = profile(cur_profile).data_path;
  % .data_support_path = where the support data is stored, the
  %   make_gps_* functions have relative paths to the data from here
  gRadar.data_support_path = profile(cur_profile).data_support_path;
  % .support_path = where the support data is stored, this includes
  %   radar configs, gps, vectors, records, frames, and GIS files
  gRadar.support_path = profile(cur_profile).support_path;
  % .out_path = where the output directories are stored (i.e. quick look,
  %   csarp, combine, pick, and post results)
  gRadar.out_path = profile(cur_profile).out_path;
  % .GIS_path = where GIS files are stored (e.g. Landsat-7 imagery)
  gRadar.gis_path = profile(cur_profile).gis_path;
  % .slurm_jobs_path = where param structures are stored for slurm scripts
  if ~isfield(profile(cur_profile),'slurm_jobs_path')
    profile(cur_profile).slurm_jobs_path = '';
  end  
  gRadar.slurm_jobs_path = profile(cur_profile).slurm_jobs_path;  
  
  clear profile cur_profile fn_dir fn_idx fn_name fns pidx;

else
  % fprintf('Compiling code: not running addpath in the compiled code\n');
end

return;

