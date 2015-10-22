function qlook_accum(param, param_override)
% qlook_accum(param, param_override)
%
% Qlook processing function called from master (for accumulation radar).
%
% Author: John Paden, Sam Buchanan
%
% See also: master

% =====================================================================
% General Setup
% =====================================================================
%clear; % Useful when running as script
%close all; % Optional
fprintf('\n\n==============================================\n\n');

% =====================================================================
% User Settings
% =====================================================================
%param = []; % Uncomment if running as a script
if ~exist('param','var') || isempty(param)
  param = read_param_xls('/mnt/scratch2/csarp_support/documents/accum_param_2011_Greenland_P3.xls','20110412_01');
  
  clear('param_override');
  %param_override.sched.type = 'no scheduler';
  
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
fprintf('\n=============================================================\n');
fprintf('Accum qlook processing for %s (%s)\n', param.day_seg, datestr(now));
fprintf('=============================================================\n');

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
  param.qlook.plot.en = 0; % Turn qlook plotting off
  
  % Prepare submission ctrl structure for queing jobs
  ctrl.cmd = 'task';
end

% =======================================================================
% Loop through and process each file, store processed output
%  (one output file per input file, one echogram jpg per input file, one
%   map jpg per input file, create pdfs and kml browse files)
%   auto-track surface for optimal echogram jpg
% =======================================================================

base_dir = fullfile(param.vectors.file.base_dir,param.vectors.file.adc_folder_name);

fprintf('Getting files for %s (%s)\n', base_dir, datestr(now));
fns = get_filenames(base_dir,param.vectors.file.file_prefix,'','.dat');
fprintf('  Found %d files in %s\n', length(fns), base_dir);

if param.qlook.file.stop_idx > length(fns)
  stop_idx = length(fns);
else
  stop_idx = param.qlook.file.stop_idx;
end
file_idxs = param.qlook.file.start_idx:stop_idx;

if isempty(param.cmd.frms)
  % If frames to process is empty, process all frames
  for file_num = 1:length(file_idxs)
    fn = fns{file_idxs(file_num)};
    finfo = fname_info_accum(fn);
    param.cmd.frms(end+1) = finfo.file_idx;
  end
end

% Get small records file info for start/stop idxs
records_fn = ct_filename_support(param,'','records');
[path name] = fileparts(records_fn);
records_small_fn = fullfile(path, sprintf('small_%s.mat', name));
if exist(records_small_fn, 'file')
  load(records_small_fn);
else
  records = struct([]);
end

if isfield(records, 'ver') && records.ver >= 2 ...
    && isfield(records, 'radar_name') && strcmpi(records.radar_name, 'accum')
  param.qlook.records_en = 1;
else
  fprintf('  Records file not available\n');
  param.qlook.records_en = 0;
end

% Get small records file info for start/stop idxs
frames_fn = ct_filename_support(param,'','frames');
if exist(frames_fn, 'file')
  frames = load(frames_fn);
else
  frames.proc_mode = [];
end

for file_num = 1:length(file_idxs)
  fn = fns{file_idxs(file_num)};
  finfo = fname_info_accum(fn);
  frm = finfo.file_idx;
  
  if param.qlook.records_en
    param.qlook.load_idxs = [records.relative_rec_num(file_num) ...
      records.relative_rec_num(file_num+1)-1];
  end
  
  if ~isempty(find(param.cmd.frms == frm, 1)) ...
      && (isempty(frames.proc_mode) || mod(floor(frames.proc_mode(frm+1)/10),10) ~= 2)
    % =================================================================
    % Execute tasks/jobs
    fh = @qlook_accum_task;
    arg{1} = fn;
    arg{2} = param;
    
    if ~strcmp(param.sched.type,'no scheduler')
      [ctrl,job_id,task_id] = create_task(ctrl,fh,1,arg);
      fprintf('  File %d of %d, idx %d in job,task %d,%d (%s)\n', ...
        file_num, length(file_idxs), file_idxs(file_num), job_id, task_id, datestr(now));
      retry_fields{job_id,task_id}.file_num = file_num;
      retry_fields{job_id,task_id}.arg = arg;
    else
      fprintf('  File %d of %d, idx %d (%s)\n', ...
        file_num, length(file_idxs), file_idxs(file_num), datestr(now));
      [success] = fh(arg{1},arg{2});
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
        fprintf('  File %d of %d, idx %d in job,task %d,%d (%s)\n', ...
          old_retry_fields{job_idx,task_idx}.file_num, length(file_idxs), ...
          file_idxs(old_retry_fields{job_idx,task_idx}.file_num), ...
          job_id, task_id, datestr(now));
        retry_fields{job_id,task_id}.file_num = file_num;
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

return;
