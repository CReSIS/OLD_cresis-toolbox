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
  cur_profile = 1;
  
  fprintf('Startup Script Running\n');
  
  format short; format compact;
  
  %% KU Profile Linux (PROFILE 1)
  % ----------------------------------------------------------------------
  pidx = 1; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = '/users/paden/scripts/matlab/';
  profile(pidx).ct_path                   = '/users/paden/scripts/cresis-toolbox/cresis-toolbox/';
  profile(pidx).param_path                = '/users/paden/scripts/params/';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = '/cresis/snfs1/scratch/paden/mdce_tmp/';
  profile(pidx).ct_tmp_file_path          = '/cresis/snfs1/dataproducts/ct_data/ct_tmp';
  
  profile(pidx).data_path                 = '/cresis/snfs1/data/';
  profile(pidx).data_support_path         = '/cresis/snfs1/dataproducts/metadata/';
  profile(pidx).support_path              = '/cresis/snfs1/dataproducts/csarp_support/';
  profile(pidx).out_path                  = '/cresis/snfs1/dataproducts/ct_data/';
  profile(pidx).gis_path                  = '/cresis/snfs1/dataproducts/GIS_data/';
  
  profile(pidx).sched.type                = 'custom_torque';
  %profile(pidx).sched.type                = 'local'; % local parallel processing
  %profile(pidx).sched.type                = 'no scheduler'; % no parallel processing (DEBUG MODE)
  profile(pidx).sched.ver                 = 2; % local and jobmanager only
  %profile(pidx).sched.name                = 'kjm';
  %profile(pidx).sched.url                 = 'kraken.cluster.cresis.ku.edu';
  profile(pidx).sched.data_location       = '/cresis/snfs1/scratch/paden/mdce_tmp/torque-temp';
  %profile(pidx).sched.submit_arguments    = '-l nodes=1:ppn=1,pmem=1725mb,walltime=120:00'; % 4 processes per node
  %profile(pidx).sched.submit_arguments    = '-l nodes=1:ppn=1,pmem=2300mb,walltime=120:00'; % 3 processes per node
  profile(pidx).sched.submit_arguments    = '-l nodes=1:ppn=1,pmem=3450mb,walltime=120:00'; % 2 processes per node
  %profile(pidx).sched.submit_arguments    = '-l nodes=1:ppn=1,pmem=6500mb,walltime=120:00'; % 1 process per node
  profile(pidx).sched.max_in_queue        = 64;
  profile(pidx).sched.cluster_size        = inf;
  profile(pidx).sched.stop_on_fail        = true;
  profile(pidx).sched.max_retries         = 4;
  profile(pidx).sched.worker_fn           = '/users/paden/scripts/cresis-toolbox/cresis-toolbox-torque/worker';
  profile(pidx).sched.force_compile       = false;
  profile(pidx).sched.rerun_only          = false;
  
  %% IU Profile Linux (PROFILE 2)
  % ----------------------------------------------------------------------
  pidx = 2; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = '/N/u/paden/Quarry/scripts/matlab/';
  profile(pidx).ct_path                   = '/N/u/paden/Quarry/scripts/cresis-toolbox/cresis-toolbox/';
  profile(pidx).param_path                = '/N/u/paden/Quarry/scripts/params/';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = '/N/dc2/scratch/paden/mdce_tmp/';
  profile(pidx).ct_tmp_file_path          = '/N/dc2/projects/cresis/ct_tmp/';
  
  profile(pidx).data_path                 = '/N/dc2/projects/cresis/';
  profile(pidx).data_support_path         = '/N/dc2/projects/cresis/metadata/';
  profile(pidx).support_path              = '/N/dc2/projects/cresis/csarp_support/';
  profile(pidx).out_path                  = '/N/dc2/projects/cresis/output/';
  profile(pidx).gis_path                  = '/N/dc2/projects/cresis/GIS_data';
  
  profile(pidx).sched.type                = 'custom_torque';
  %profile(pidx).sched.type                = 'local'; % local parallel processing
  %profile(pidx).sched.type                = 'no scheduler'; % no parallel processing (DEBUG MODE)
  profile(pidx).sched.ver                 = 2; % local and jobmanager only
  profile(pidx).sched.name                = '';
  profile(pidx).sched.url                 = 'qm2.quarry.teragrid.iu.edu';
  profile(pidx).sched.data_location       = '/N/dc2/scratch/paden/matlab_torque/';
  profile(pidx).sched.submit_arguments    = '-q cresis -l nodes=1:ppn=1:dc2,pmem=2gb,walltime=40:00';
  profile(pidx).sched.max_in_queue        = 64;
  profile(pidx).sched.max_tasks_per_jobs  = 64;
  profile(pidx).sched.cluster_size        = inf;
  profile(pidx).sched.stop_on_fail        = true;
  profile(pidx).sched.max_retries         = 4;
  profile(pidx).sched.worker_fn           = '/N/u/paden/Quarry/scripts/cresis-toolbox-torque/worker';
  profile(pidx).sched.force_compile       = false;
  profile(pidx).sched.rerun_only          = false;
  
  %% Field Profile Linux (PROFILE 3)
  % ----------------------------------------------------------------------
  pidx = 3; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = '/scratch/scripts/matlab/';
  profile(pidx).ct_path                   = '/scratch/scripts/cresis-toolbox/cresis-toolbox/';
  profile(pidx).param_path                = '/scratch/scripts/params/';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = '/scratch/tmp/';
  profile(pidx).ct_tmp_file_path          = '/scratch/ct_tmp/';
  
  profile(pidx).data_path                 = '/scratch/';
  profile(pidx).data_support_path         = '/scratch/metadata/';
  profile(pidx).support_path              = '/scratch/csarp_support/';
  profile(pidx).out_path                  = '/scratch/';
  profile(pidx).gis_path                  = '/scratch/GIS_data';
  
  profile(pidx).sched.type                = 'custom_torque';
  %profile(pidx).sched.type                = 'local'; % local parallel processing
  %profile(pidx).sched.type                = 'no scheduler'; % no parallel processing (DEBUG MODE)
  profile(pidx).sched.ver                 = 2; % local and jobmanager only
  profile(pidx).sched.data_location       = '/scratch/tmp/torque-temp';
  profile(pidx).sched.submit_arguments    = '-l nodes=1:ppn=2,walltime=15:00';
  profile(pidx).sched.max_in_queue        = 64;
  profile(pidx).sched.cluster_size        = inf;
  profile(pidx).sched.stop_on_fail        = true;
  profile(pidx).sched.max_retries         = 4;
  profile(pidx).sched.worker_fn           = '/scratch/cresis-toolbox-torque/worker';
  profile(pidx).sched.force_compile       = false;
  profile(pidx).sched.rerun_only          = false;
  
  %% KU Mobile Profile Windows (PROFILE 5)
  % ----------------------------------------------------------------------
  pidx = 5; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = 'C:\Users\paden\Documents\scripts\matlab\';
  profile(pidx).ct_path                   = 'C:\Users\paden\Documents\scripts\cresis-toolbox\cresis-toolbox\';
  profile(pidx).param_path                = 'C:\Users\paden\Documents\scripts\params\';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = 'C:\tmp\mdce_tmp\';
  profile(pidx).ct_tmp_file_path          = 'D:\output\ct_tmp\';
  
  profile(pidx).data_path                 = 'D:\';
  profile(pidx).data_support_path         = 'C:\metadata\';
  profile(pidx).support_path              = 'C:\csarp_support\';
  profile(pidx).out_path                  = 'D:\output\';
  profile(pidx).gis_path                  = 'C:\GIS_data\';
  
  profile(pidx).sched.type                = 'local';
  %profile(pidx).sched.type                = 'no scheduler'; % no parallel processing (DEBUG MODE)
  profile(pidx).sched.ver                 = 2; % local and jobmanager only
  %profile(pidx).sched.name                = '';
  %profile(pidx).sched.url                 = '';
  profile(pidx).sched.data_location       = '';
  profile(pidx).sched.submit_arguments    = '';
  profile(pidx).sched.max_in_queue        = inf;
  profile(pidx).sched.cluster_size        = inf;
  profile(pidx).sched.stop_on_fail        = true;
  profile(pidx).sched.max_retries         = 4;
  profile(pidx).sched.worker_fn           = '';
  profile(pidx).sched.force_compile       = false;
  profile(pidx).sched.rerun_only          = false;
  
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
  profile(pidx).param_path                = 'C:\Users\paden\Documents\scripts\params\';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = 'Y:/paden/mdce_tmp/';
  profile(pidx).ct_tmp_file_path          = 'X:/ct_data/ct_tmp/';
  
  profile(pidx).data_path                 = 'V:/';
  profile(pidx).data_support_path         = 'X:/metadata/';
  profile(pidx).support_path              = 'X:/csarp_support/';
  profile(pidx).out_path                  = 'X:/ct_data/';
  profile(pidx).gis_path                  = 'X:/GIS_data/';
  
  %profile(pidx).sched.type                = 'local';
  profile(pidx).sched.type                = 'no scheduler'; % no parallel processing (DEBUG MODE)
  profile(pidx).sched.ver                 = 2; % local and jobmanager only
  %profile(pidx).sched.name                = '';
  %profile(pidx).sched.url                 = '';
  profile(pidx).sched.data_location       = '';
  profile(pidx).sched.submit_arguments    = '';
  profile(pidx).sched.max_in_queue        = inf;
  profile(pidx).sched.cluster_size        = inf;
  profile(pidx).sched.stop_on_fail        = true;
  profile(pidx).sched.max_retries         = 4;
  profile(pidx).sched.worker_fn           = '';
  profile(pidx).sched.force_compile       = false;
  profile(pidx).sched.rerun_only          = false;
  
  %% AWI Profile Windows (PROFILE 7)
  % ----------------------------------------------------------------------
  pidx = 7; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = 'C:\tmp\scripts\matlab\';
  profile(pidx).ct_path                   = 'C:\tmp\scripts\cresis-toolbox\cresis-toolbox\';
  profile(pidx).param_path                = 'C:\tmp\scripts\params\';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = 'F:\mdce_tmp\';
  profile(pidx).ct_tmp_file_path          = 'F:\ct_tmp\';
  
  profile(pidx).data_path                 = 'D:\';
  profile(pidx).data_support_path         = 'F:\metadata\';
  profile(pidx).support_path              = 'F:\csarp_support\';
  profile(pidx).out_path                  = 'F:\';
  profile(pidx).gis_path                  = 'C:\tmp\GIS_data\';
  
  profile(pidx).sched.type                = 'local';
  profile(pidx).sched.ver                 = 2; % local and jobmanager only
  profile(pidx).sched.data_location       = '';
  profile(pidx).sched.submit_arguments    = '';
  profile(pidx).sched.max_in_queue        = 256;
  profile(pidx).sched.max_tasks_per_jobs  = 256;
  profile(pidx).sched.cluster_size        = inf;
  profile(pidx).sched.stop_on_fail        = true;
  profile(pidx).sched.max_retries         = 4;
  profile(pidx).sched.worker_fn           = '';
  profile(pidx).sched.force_compile       = false;
  profile(pidx).sched.rerun_only          = false;
  
  %% AWI Profile Linux (PROFILE 8)
  % ----------------------------------------------------------------------
  pidx = 8; % profile index
  profile(pidx).debug_level               = 1;
  profile(pidx).personal_path             = '/home/administrator/scripts/matlab/';
  profile(pidx).ct_path                   = '/home/administrator/scripts/cresis-toolbox/';
  profile(pidx).param_path                = '/home/administrator/scripts/params/';
  
  profile(pidx).code_path                 = profile(pidx).ct_path;
  profile(pidx).code_path_override        = profile(pidx).personal_path;
  profile(pidx).tmp_file_path             = '/home/administrator/Scratch/mdce_tmp/';
  profile(pidx).ct_tmp_file_path          = '/home/administrator/Scratch/ct_tmp/';
  
  profile(pidx).data_path                 = '/mnt/AWI_SSD0/';
  profile(pidx).data_support_path         = '/home/administrator/Scratch/metadata/';
  profile(pidx).support_path              = '/home/administrator/Scratch/csarp_support/';
  profile(pidx).out_path                  = '/home/administrator/Scratch/';
  profile(pidx).gis_path                  = '/home/administrator/GIS_data/';
  
  profile(pidx).sched.type                = 'custom_torque';
  profile(pidx).sched.ver                 = 2; % local and jobmanager only
  profile(pidx).sched.data_location       = '/home/administrator/Scratch/torque-temp';
  profile(pidx).sched.submit_arguments    = '-l nodes=1:ppn=2,walltime=15:00';
  profile(pidx).sched.max_in_queue        = 64;
  profile(pidx).sched.cluster_size        = inf;
  profile(pidx).sched.stop_on_fail        = true;
  profile(pidx).sched.max_retries         = 4;
  profile(pidx).sched.worker_fn           = '/home/administrator/scripts/cresis-toolbox-torque/worker';
  profile(pidx).sched.force_compile       = false;
  profile(pidx).sched.rerun_only          = false;
    
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
  
  % .debug_level = higher means more print outs, 0 is none, 1 is minimal
  %   (the preferred default)
  gRadar.debug_level = profile(cur_profile).debug_level;
  % .sched.type = 'local', 'jobmanager', 'torque', or 'no scheduler'
  %   Currently KU runs jobmanager, IU runs torque
  %   local uses only the cores on your local machine
  gRadar.sched.type = profile(cur_profile).sched.type;
  % .sched.url = url of scheduler (e.g. qm2.quarry.teragrid.iu.edu or
  %   heimdall.cluster.cresis.ku.edu)
  if ~isfield(profile(cur_profile).sched,'url')
    profile(cur_profile).sched.url = '';
  end
  gRadar.sched.url = profile(cur_profile).sched.url;
  % .sched.name = name of scheduler (often leave blank, sometimes 'jm' or
  %   'jm2')
  if ~isfield(profile(cur_profile).sched,'name')
    profile(cur_profile).sched.name = '';
  end
  gRadar.sched.name = profile(cur_profile).sched.name;
  % .sched.cluster_size = forces the number of submitted tasks to never
  %   exceed this number (i.e. even if there are more nodes available, it
  %   will not use more than this number, set to inf to use all)
  gRadar.sched.cluster_size = profile(cur_profile).sched.cluster_size;
  % .sched.stop_on_fail = enter debug mode if failure occurs
  gRadar.sched.stop_on_fail = profile(cur_profile).sched.stop_on_fail;
  % .sched.data_location = torque scheduler only, this is where the
  %   temporary job files will be stored
  gRadar.sched.data_location = profile(cur_profile).sched.data_location;
  % .sched.submit_arguments = torque scheduler only, these are the
  %   additional arguments passed to the qsub command (qsub is the queue
  %   submit command for the torque scheduler)
  gRadar.sched.submit_arguments = profile(cur_profile).sched.submit_arguments;
  % .sched.max_in_queue = same as cluster_size, but is used by the
  %   wrapper function create_task. This function is necessary when
  %   submitting thousands of jobs at once.  Matlab does not handle
  %   large numbers of jobs (tasks) elegantly and slows down too much
  %   and can crash the jobmanager.
  gRadar.sched.max_in_queue = profile(cur_profile).sched.max_in_queue;
  % .sched.max_tasks_per_jobs = unfortunately Matlab's torque scheduler
  %   interface performs slowly when the number of tasks exceed about 50
  %   tasks per job so you can force it to submit just a few (e.g. 50)
  %   tasks per job and then submit multiple jobs
  if ~isfield(profile(cur_profile).sched,'max_tasks_per_jobs')
    profile(cur_profile).sched.max_tasks_per_jobs = '';
  end
  gRadar.sched.max_tasks_per_jobs = profile(cur_profile).sched.max_tasks_per_jobs;
  % .sched.max_retries = Matlab's torque scheduler occasionally
  %   submits tasks that don't ever execute and no error information
  %   is returned except for the fact that the process does not return any
  %   output arguments.  In this case, the toolbox will try to resubmit the
  %   job up to this many times.
  gRadar.sched.max_retries = profile(cur_profile).sched.max_retries;
  % .sched.worker_fn = worker shell script filename
  %   qsub job submissions all call this worker shell script which
  %   then calls an MCC compiled worker_task.m
  gRadar.sched.worker_fn = profile(cur_profile).sched.worker_fn;
  % .sched.force_compile = boolean which controls whether or not
  %   the torque scheduler will force a recompile before submitting
  %   a batch job (typically this should be set to false because
  %   the function checks the time stamp of dependent functions and
  %   only recompiles when it is necessary to do so)
  gRadar.sched.force_compile = profile(cur_profile).sched.force_compile;
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
  gRadar.sched.hidden_depend_funs = {};
  gRadar.sched.hidden_depend_funs{end+1} = {'create_records_accum2_task.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'create_records_acords_task.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'create_records_mcords_task.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'basic_load_fmcw.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'basic_load_fmcw2.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'rx_chan_equal_sar_task.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'rx_chan_equal_raw_task.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'coh_noise_tracker_task.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'analysis_task.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'radiometric_calibration_task.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'get_heights_task.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'csarp_task.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'combine_wf_chan_task.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'nsidc_delivery_script_task.m' 2};
  gRadar.sched.hidden_depend_funs{end+1} = {'hanning.m' 0};
  gRadar.sched.hidden_depend_funs{end+1} = {'hamming.m' 0};
  gRadar.sched.hidden_depend_funs{end+1} = {'blackman.m' 0};
  gRadar.sched.hidden_depend_funs{end+1} = {'tukeywin.m' 0};
  gRadar.sched.hidden_depend_funs{end+1} = {'tukeywin_trim.m' 1};
  gRadar.sched.hidden_depend_funs{end+1} = {'chebwin.m' 0};
  gRadar.sched.hidden_depend_funs{end+1} = {'kaiser.m' 0};
  gRadar.sched.hidden_depend_funs{end+1} = {'boxcar.m' 0};
  gRadar.sched.hidden_depend_funs{end+1} = {'butter.m' 0};
  gRadar.sched.hidden_depend_funs{end+1} = {'array_proc_sv.m' 1};
  gRadar.sched.hidden_depend_funs{end+1} = {'lever_arm.m' 1};
  gRadar.sched.hidden_depend_funs{end+1} = {'doa_nonlcon.m' 1};
  % For newer Matlab cluster access, we need a list of all the paths
  gRadar.sched.AdditionalPaths = AdditionalPaths;
  clear AdditionalPaths;
  % For newer Matlab we need to know the version of the cluster interface
  gRadar.sched.ver = profile(cur_profile).sched.ver;

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
  
  clear profile cur_profile;

else
  % fprintf('Compiling code: not running addpath in the compiled code\n');
end

return;

