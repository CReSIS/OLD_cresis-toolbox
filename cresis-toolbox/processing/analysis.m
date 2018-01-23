function analysis(param,param_override)
% analysis(param,param_override)
%
% Function separates radar data for one day_seg into blocks which are
% processed by analysis_task.m.
%
% This function has two capabilities (enabled separately):
% 1.  Break radar data specified by day_seg into slow time blocks for
%     evaluating power spectral density of noise (performed by
%     analysis_task.m).
%     -> For evaluating psd, block size is determined by number of
%        incoherent averages specified by user in param file.
% 2.  Break radar data specified by day_seg into slow time blocks for
%     evaluating noise power of each range line in block (performed by
%     analysis_task.m).
%     ->  For evaluating np, block size is determined by np.analysis_frm_size
%         field of param specified by user.
%
%
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
% Authors: Theresa Stumpf, John Paden
%
% See also: master_mcords.m, analysis_task.m
% =====================================================================
% General Setup
% =====================================================================
% clear; % Useful when running as script
%close all; % Optional
tic;
fprintf('\n\n==============================================\n\n');

% =====================================================================
% User Settings
% =====================================================================
% param = []; % Uncomment if running as a script
if ~exist('param','var') || isempty(param)
  %param = read_param_xls('E:\mcords_param_2010_Antarctica_DC8.xls','20101013_seg4');
  param = read_param_xls('/users/tstumpf/scripts/matlab/mcords_param_2010_Greenland_P3.xls','20100510_04');
  
  clear('param_override');
  param_override.sched.type = 'no scheduler';
  
  % Input checking
  if ~exist('param','var')
    error('A struct array of parameters must be passed in\n');
  end
  global gRadar;
  if exist('param_override','var')
    param_override = merge_structs(gRadar,param_override);
  else
    param_override = gRadar;
  end
  
elseif ~isstruct(param)
  % Functional form
  param();
end
param = merge_structs(param, param_override);

% =====================================================================
% Setup processing
% =====================================================================
fprintf('=====================================================================\n');
fprintf('Analysis %s (%s)\n', param.day_seg, datestr(now));
fprintf('=====================================================================\n');

physical_constants;

param.analysis.file.base_dir         = param.vectors.file.base_dirs;
param.analysis.file.adc_folder_name  = param.vectors.file.adc_folder_name;
param.analysis.file.file_prefix      = param.vectors.file.file_prefix;


% Load frames file
load(ct_filename_support(param,param.frames.frames_fn,'frames'));

if ~isfield(frames,'records_fn')
  frames.records_fn = '';
end

% Load records file
load(ct_filename_support(param,frames.records_fn,'records'));

% -------------------------------------------------------------------------
% Create analysis processing blocks
% -------------------------------------------------------------------------

% Determine total records being processed
% If no records are specified, process all records for the particular day seg
if isempty(param.analysis.records)
  param.analysis.records = [1 length(records.time)];
end

if isinf(param.analysis.records(2))
  param.analysis.records(2) = length(records.time);
end

% Determine block size
% -------------------------------------------------------------------------
% SPECIAL CASES:
%   1.  param.analysis.block_size = inf
%       or param.analysis.block_size = 0
%       param.analysis.np.block_size = inf
%       or param.analysis.np.block_size = 0
%           -->> process all records in one block
%
%   2. param.analysis.records(2) ~= num_blocks * block_size
%           -->> process maximum number of processing blocks

block_size = param.analysis.block_size;

if isinf(block_size)
  block_size = (param.analysis.records(2) - param.analysis.records(1))+1;
end

if isempty(block_size)
  block_size = (param.analysis.records(2) - param.analysis.records(1))+1;
end

if (param.analysis.records(2) - param.analysis.records(1)) < block_size
  block_size = (param.analysis.records(2) - param.analysis.records(1))+1;
end

param.analysis.block_size = block_size;
num_blks = floor(((param.analysis.records(2)-param.analysis.records(1))+1)/block_size);

if strcmpi(param.analysis.analysis_type,'psd')
  if (param.analysis.records(1) + (num_blks*block_size)) < param.analysis.records(2)
    param.analysis.records(2) = param.analysis.records(1) + (num_blks - 1)*block_size;
    blocks = param.analysis.records(1):block_size:param.analysis.records(2);
  else
    blocks = param.analysis.records(1):block_size:param.analysis.records(2);
  end
end

if strcmpi(param.analysis.analysis_type,'np')
  block_size = 1000;
  blocks = param.analysis.records(1):block_size:param.analysis.records(2);
end


global g_data;
g_data = [];

% -------------------------------------------------------------------------
% Clean up old directories, create new paths and store analysis parameters
% -------------------------------------------------------------------------

if strcmpi(param.analysis.analysis_type,'np')
  param.analysis.coh_ave = 1;
end

% Specify output directory
param.analysis.base_dir =  ct_filename_out(param, ...
  param.analysis.out_path, 'CSARP_analysis');

if strcmpi(param.analysis.analysis_type,'psd')
  for blk_idx = 1:length(blocks)
    start = blocks(blk_idx);
    stop  = start + block_size -1;
    recs = [start stop];
    
    for img_idx = 1:size(param.analysis.imgs,2)
      bins = param.analysis.ft_bins(img_idx,:);
      wf = param.analysis.imgs{img_idx}(1,1);
      
      
      coh_ave = param.analysis.psd.coh_ave;
      
      
      old_psd_dir = fullfile(param.analysis.base_dir,sprintf('psd_analysis/recs_%010d_%010d/bins_%04d_%04d/coh_ave_%05d/wf_%02d/',recs(1),recs(2),bins(1),bins(2),coh_ave,wf));
      
      if exist(old_psd_dir,'dir')
        rmdir(old_psd_dir,'s');
        fprintf('Removing path: %s\n', old_psd_dir);
      end
      
      psd_analysis_dir = fullfile(param.analysis.base_dir,sprintf('psd_analysis/recs_%010d_%010d/bins_%04d_%04d/coh_ave_%05d/wf_%02d/',recs(1),recs(2),bins(1),bins(2),coh_ave,wf));
      
      if ~exist(psd_analysis_dir,'dir')
        mkdir(psd_analysis_dir);
      end
      
    end
  end
end

% Clean up preexisting noise power directories
if strcmpi(param.analysis.analysis_type,'np')
  
  recs = [param.analysis.records(1) param.analysis.records(2)];
  old_np_dir = fullfile(param.analysis.base_dir,sprintf('np_analysis/recs_%010d_%010d',recs(1),recs(end)));
  
  if exist(old_np_dir,'dir')
    rmdir(old_np_dir,'s')
    fprintf('Removing path: %s\n', old_np_dir);
  end
  
  % Create new tmp directory for output of noise power processing blocks
  np_analysis_dir = fullfile(param.analysis.base_dir,sprintf('np_analysis/recs_%010d_%010d',recs(1),recs(2)));
  
  if ~exist(np_analysis_dir,'dir')
    mkdir(np_analysis_dir);
  end
  
  np_tmp_dir = fullfile(np_analysis_dir,sprintf('tmp/'));
  
  if ~exist(np_tmp_dir,'dir')
    mkdir(np_tmp_dir);
  end
  
end

% =====================================================================
% Setup static inputs for analysis_task
% =====================================================================
global gRadar

task_param.gRadar           = gRadar;
task_param.ft_wind_time     = false;
task_param.radar_name       = param.radar_name;
task_param.season_name      = param.season_name;
task_param.day_seg          = param.day_seg;
task_param.load.imgs        = param.analysis.imgs;
task_param.radar            = param.radar;
task_param.profile.out_path = param.out_path;

if ~isfield(param,'debug_level')
  task_param.debug_level = 1;
else
  task_param.debug_level = param.debug_level;
end

task_param.load.records_fn  = ct_filename_support(param,'','records');

% Currently every field in param.analysis is used by analysis_task
% so we just pass the whole structure
task_param.analysis     = param.analysis;
task_param.analysis     = rmfield(task_param.analysis, {'imgs'});

% =====================================================================
% Setup the scheduler
% =====================================================================
fd = [get_filenames(param.path,'','','.m',struct('recursive',1)); ...
  get_filenames(param.path,'','','.mexa64',struct('recursive',1))];
fd_override = [get_filenames(param.path_override,'','','.m',struct('recursive',1)); ...
  get_filenames(param.path_override,'','','.mexa64',struct('recursive',1))];

fd = merge_filelists(fd, fd_override);

if ~strcmpi(param.sched.type,'no scheduler')
  % Initialize submission ctrl structure
  ctrl.cmd = 'init';
  ctrl.sched = param.sched;
  ctrl.fd = fd;
  ctrl = create_task(ctrl);
  
  % Prepare submission ctrl structure for queing jobs
  ctrl.cmd = 'task';
end

% =====================================================================
% Load data and run analysis tasks
% For each frame, load block_size records at a time (code groups by
% file index, but has to watch negative offset values which imply the
% record starts in a previous file and carries over into the next
%   --> The last block can be up to 2*REC_BLOCK_SIZE
% =====================================================================
out_recs = {};
retry_fields = {};

for block_idx = 1:length(blocks)
  block = blocks(block_idx);
  
  fprintf('Processing block %d of %d\n', block_idx, length(blocks));
  
  % Find record numbers associated with frame boundaries
  cur_recs = block + [0 block_size-1];
  
  if block_idx == length(blocks) && block_idx < block + block_size - 1
    stop_rec = param.analysis.records(2);
    cur_recs = [block stop_rec];
  end
    
   
  % =====================================================================
  % Prepare task inputs
  % =====================================================================
  task_param.load.recs = cur_recs;
  
  for img_idx = 1:length(param.analysis.imgs)
    wf  = param.analysis.imgs{img_idx}(1);
    task_param.analysis.ft_bins = ...
      param.analysis.ft_bins(img_idx,:);
    for adc_idx = 1:size(param.analysis.imgs{img_idx})
      adc = param.analysis.imgs{img_idx}(adc_idx,2);
      task_param.load.imgs = {param.analysis.imgs{img_idx}(adc_idx,:)};
      
      % =================================================================
      % Execute tasks/jobs
      fh = @analysis_task;
      arg{1} = task_param;
      
      if ~strcmp(param.sched.type,'no scheduler')
        [ctrl,job_id,task_id] = create_task(ctrl,fh,1,arg);
        fprintf('  ADC %d Wf %d Records %d to %d in job,task %d,%d (%s)\n', ...
          adc, wf, cur_recs(1), cur_recs(end), job_id, task_id, datestr(now));
        retry_fields{job_id,task_id}.cur_recs = cur_recs;
        retry_fields{job_id,task_id}.adc = adc;
        retry_fields{job_id,task_id}.wf = wf;
        retry_fields{job_id,task_id}.arg = arg;
        retry_fields{job_id,task_id}.block_idx = block_idx;
        if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
          % Quit if a bad error occurred
          fprintf('Bad errors occurred, quitting (%s)\n\n', datestr(now));
          ctrl.cmd = 'done';
          ctrl = create_task(ctrl);
          return;
        end
      else
        fprintf('  ADC %d Wf %d Records %d to %d (%s)\n', ...
          adc, wf, cur_recs(1), cur_recs(end), datestr(now));
        success = fh(arg{1});
      end
      
    end
  end
  
end

% =======================================================================
% Wait for jobs to complete if a scheduler was used
% =======================================================================
if ~strcmpi(param.sched.type,'no scheduler')
  ctrl.cmd = 'done';
  ctrl = create_task(ctrl);
  if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
    % Quit if a bad error occurred
    fprintf('Bad errors occurred, quitting (%s)\n\n', datestr(now));
    return;
  end
  
  retry = 1;
  while ctrl.error_mask == 2 && retry <= param.sched.max_retries
    fprintf('Tasks failed, retry %d of max %d\n', retry, param.sched.max_retries);
    
    % Bookkeeping (move previous run info to "old_" variables)
    old_ctrl = ctrl;
    old_retry_fields = retry_fields;
    retry_fields = {};
    
    % Initialize submission ctrl structure
    ctrl = [];
    ctrl.cmd = 'init';
    ctrl.sched = param.sched;
    ctrl.fd = fd;
    ctrl = create_task(ctrl);
    
    % Prepare submission ctrl structure for queing jobs
    ctrl.cmd = 'task';
    
    for job_idx = 1:length(old_ctrl.jobs)
      for task_idx = old_ctrl.jobs{job_idx}.error_idxs
        [ctrl,job_id,task_id] = create_task(ctrl,fh,1,old_retry_fields{job_idx,task_idx}.arg);
        fprintf('  %d ADC %d Wf %d Records %d to %d in job,task %d,%d (%s)\n', ...
          old_retry_fields{job_idx,task_idx}.block_idx, ...
          old_retry_fields{job_idx,task_idx}.adc, ...
          old_retry_fields{job_idx,task_idx}.wf, ...
          old_retry_fields{job_idx,task_idx}.arg, ...
          old_retry_fields{job_idx,task_idx}.cur_recs(1), ...
          old_retry_fields{job_idx,task_idx}.cur_recs(end), ...
          job_id, task_id, datestr(now));
        retry_fields{job_id,task_id} = old_retry_fields{job_idx,task_idx};
      end
    end
    ctrl.cmd = 'done';
    ctrl = create_task(ctrl);
    retry = retry + 1;
  end
  if ctrl.error_mask ~= 0
    fprintf('Not all jobs completed, but out of retries (%s)\n', datestr(now));
    return;
  else
    fprintf('Jobs completed (%s)\n\n', datestr(now));
  end
end

% Concatenate Noise Power Blocks
% -------------------------------------------------------------------------
if strcmpi(param.analysis.analysis_type,'np')
  
  % Build param_analyze
param_radar                       = param.radar;
param_analysis                    = param.analysis;
param_radar.radar_name            = param.radar_name;
param_analysis.season_name        = param.season_name;
param_analysis.day_seg            = param.day_seg;
  
  analysis_recs = param.analysis.records(1):param.analysis.records(2);
  np_analysis_dir = fullfile(param.analysis.base_dir, ...
    sprintf('np_analysis/recs_%010d_%010d',param.analysis.records(1),param.analysis.records(2)));
  np_tmp_dir = fullfile(np_analysis_dir, sprintf('tmp/'));
  
  for img_idx = 1:length(param.analysis.imgs)
    wf = param.analysis.imgs{img_idx}(1);
    tmp_wf_dir = fullfile(np_tmp_dir,sprintf('wf_%02d/',wf));
    analysis_wf_dir = fullfile(np_analysis_dir, sprintf('wf_%02d/',wf));
    
    if ~exist(analysis_wf_dir,'dir')
      mkdir(analysis_wf_dir)
    end
    
    if param.analysis.gps.en
      gps_string = sprintf('np_gps_time');
      gps_files = get_filenames(tmp_wf_dir,gps_string,'','');
      np_gps_time = [];
      for gps_file_idx = 1:length(gps_files);
        load(gps_files{gps_file_idx});
        np_gps_time = [np_gps_time gps_time];
      end
      np_gps_time_fn = fullfile(analysis_wf_dir,sprintf('np_gps_time_%s_%010f',param.day_seg,np_gps_time(1)));
      save(np_gps_time_fn,'np_gps_time')
    end
    
    for adc_idx = 1:size(param.analysis.imgs{img_idx},1)
      adc = param.analysis.imgs{img_idx}(adc_idx,2);
      element = param_radar.wfs(wf).rx_paths(adc);
      adc_element_string = sprintf('adc_%02d_element_%02d',adc,element);
      np_files = get_filenames(tmp_wf_dir,adc_element_string,'','');
      np_dBm = [];
      for file_idx = 1:length(np_files)
        load(np_files{file_idx});
        np_dBm = [np_dBm meas_np_dBm];
        
        if file_idx == 1
          np_exp_dBm = exp_np_dBm;
          hw_presums = presums;
        end
      end
      np_fn = fullfile(analysis_wf_dir,...
        sprintf('np_adc_%02d_element_%02d_%s.mat',adc,element,param.day_seg));
      fprintf('  Saving noise power to file %s\n', np_fn);
      save(np_fn,'np_dBm','np_exp_dBm','analysis_recs','hw_presums','param_radar','param_analysis')
      
    end
    
  end
  % Clean up tmp directory
  fprintf('Removing path: %s\n', np_tmp_dir);
  rmdir(np_tmp_dir,'s');
end

return;

