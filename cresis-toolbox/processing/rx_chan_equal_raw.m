function rx_chan_equal_raw(param,param_override)
% rx_chan_equal_raw(param,param_override)
%
% This function has two capabilities (enabled separately)
% 1. Automated or manual pick: Gets the air/ice surface height using the
%    radar depth sounder data like an altimeter (manual mode is optional)
% 2. Generate quick look outputs
%
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Authors: John Paden
%
% See also: master.m, rx_chan_equal_raw_task.m

% =====================================================================
% General Setup
% =====================================================================

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  param = read_param_xls(ct_filename_param('rds_param_2009_Greenland_TO.xls'),'20090331_01','equal');
  
  clear('param_override');
  param_override.sched.type = 'no scheduler';
  param_override.sched.rerun_only = false;
  %param_override.cmd.frms = 7;

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

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

% =====================================================================
% Setup processing
% =====================================================================

% Get WGS84 ellipsoid parameters
physical_constants;

% Load frames file
load(ct_filename_support(param,param.records.frames_fn,'frames'));
% Load records file
records_fn = ct_filename_support(param,param.records.records_fn,'records');
records = load(records_fn);

global g_data;
g_data = [];

equal_out_path = ct_filename_out(param, param.equal.out_path, 'CSARP_equal');

if isempty(param.cmd.frms)
  param.cmd.frms = 1:length(frames.frame_idxs);
end
% Remove frames that do not exist from param.cmd.frms list
[valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
if length(valid_frms) ~= length(param.cmd.frms)
  bad_mask = ones(size(param.cmd.frms));
  bad_mask(keep_idxs) = 0;
  warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
    param.cmd.frms(find(bad_mask,1)));
  param.cmd.frms = valid_frms;
end

if ~isfield(param.sched,'rerun_only') || isempty(param.sched.rerun_only)
  param.sched.rerun_only = false;
end

% Cleanup folders
if ~param.sched.rerun_only
  if exist(equal_out_path,'dir')
    for frm = param.cmd.frms
      del_paths = get_filenames(equal_out_path,sprintf('equal_data_%03d',frm),'','',struct('type','d'));
      for idx = 1:length(del_paths)
        fprintf('Removing path: %s\n', del_paths{idx});
        rmdir(del_paths{idx},'s');
      end
    end
  end
end

% =====================================================================
% Setup static inputs for equal_task
% =====================================================================
task_param = param;
task_param.load.imgs = param.equal.imgs;

% =====================================================================
% Setup the scheduler
% =====================================================================

if strcmpi(param.sched.type,'custom_torque')
  global ctrl; % Make this global for convenience in debugging
  ctrl = torque_new_batch(param);
  fprintf('Torque batch: %s\n', ctrl.batch_dir);
  torque_compile('rx_chan_equal_raw_task.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
  
elseif ~strcmpi(param.sched.type,'no scheduler')
  fd = [get_filenames(param.path,'','','.m',struct('recursive',1)); ...
    get_filenames(param.path,'','','.mexa64',struct('recursive',1))];
  fd_override = [get_filenames(param.path_override,'','','.m',struct('recursive',1)); ...
    get_filenames(param.path_override,'','','.mexa64',struct('recursive',1))];
  
  fd = merge_filelists(fd, fd_override);
  
  % Remove SVN files from path
  non_svn_path = zeros(size(fd));
  for fn_idx = 1:length(fd)
    if isempty(strfind(fd{fn_idx},'.svn'))
      non_svn_path(fn_idx) = 1;
    end
  end
  fd = fd(non_svn_path);
  
  % Initialize submission ctrl structure
  global ctrl; % Make this global for convenience in debugging
  ctrl = [];
  ctrl.cmd = 'init';
  ctrl.sched = param.sched;
  ctrl.fd = fd;
  ctrl = create_task(ctrl);
  
  % Prepare submission ctrl structure for queing jobs
  ctrl.cmd = 'task';
end

% =====================================================================
% Load data and run equal tasks
%
% For each frame load REC_BLOCK_SIZE records at a time (code groups
% by file index, but has to watch negative offset values which imply
% the record starts in a previous file and carries over into the next)
%    --> The last block can be up to 2*REC_BLOCK_SIZE
% =====================================================================
out_recs = {};
retry_fields = {};
for frame_idx = 1:length(param.cmd.frms);
  frame = param.cmd.frms(frame_idx);
  if frames.proc_mode(frame) ~= 2
    fprintf('equal frame %s_%03i (%i of %i) %s\n', param.day_seg, frame, frame_idx, length(param.cmd.frms), datestr(now));
  else
    fprintf('Skipping frame %d (no process frame)\n', frame);
    continue;
  end
  
  equal_dir = fullfile(equal_out_path, sprintf('equal_data_%03d',frame));
  % Create the array_proc output directories
  if ~exist(equal_dir,'dir')
    mkdir(equal_dir);
  end
  
  if frame < length(frames.frame_idxs)
    recs = frames.frame_idxs(frame):frames.frame_idxs(frame+1);
  else
    recs = frames.frame_idxs(frame):length(records.lat);
  end
  
  % Determine where breaks in processing are going to occur
  rec = recs(1);
  if isempty(param.equal.block_size)
    REC_BLOCK_SIZE = 10000;
    block_overlap = 0;
  else
    if numel(param.equal.block_size) == 2
      REC_BLOCK_SIZE = param.equal.block_size(1);
      block_overlap = param.equal.block_size(2);
    else
      REC_BLOCK_SIZE = param.equal.block_size;
      block_overlap = 0;
    end
  end
  if length(recs) < 2*REC_BLOCK_SIZE
    breaks = 1;
  else
    breaks = 1:REC_BLOCK_SIZE:length(recs)-REC_BLOCK_SIZE;
  end
  
  task_param.proc.frm = frame;
  
  % Begin loading data
  for break_idx = 1:length(breaks)
    task_param.proc.block = break_idx;
    % Determine the current records being processed
    if break_idx < length(breaks)
      cur_recs_keep = [recs(breaks(break_idx)) recs(breaks(break_idx+1)-1)];
      cur_recs = [max(1,recs(breaks(break_idx))-block_overlap) ...
        recs(breaks(break_idx+1)-1)+block_overlap];
    else
      cur_recs_keep = [recs(breaks(break_idx)) recs(end)];
      cur_recs = [max(1,recs(breaks(break_idx))-block_overlap) recs(end)];
    end
    
    % =====================================================================
    % Prepare task inputs
    % =====================================================================
    task_param.load.recs = cur_recs;
    task_param.load.recs_keep = cur_recs_keep;
    
    % =================================================================
    % Rerun only mode: Test to see if we need to run this task
    if param.sched.rerun_only
      % If we are in rerun only mode AND all the get heights task output files
      % already exists, then we do not run the task
      file_exists = true;
      sub_apt_shift_idx = 1;
      sub_band_idx = 1;
      out_path = fullfile(ct_filename_out(param, ...
        param.equal.out_path, 'CSARP_equal'), ...
        sprintf('equal_%03d',frame));
      start_time_for_fn = records.gps_time(cur_recs(1));
      for img = 1:length(param.equal.imgs)
        
        out_fn = fullfile(ct_filename_out(param, ...
          param.equal.out_path, 'CSARP_equal'), ...
          sprintf('equal_%03d',frame), ...
          sprintf('blk_%d_img_%02d.mat', break_idx, img));
        if ~exist(out_fn,'file')
          file_exists = false;
        end
      end
      if file_exists
        fprintf('  %d: Already exists records %d to %d [rerun_only skipping] (%s)\n', ...
          break_idx, cur_recs(1), cur_recs(end), datestr(now));
        continue;
      end
    end
    
    % =================================================================
    % Execute tasks/jobs
    fh = @rx_chan_equal_raw_task;
    if isfield(frames,'nyquist_zone') && ~isnan(frames.nyquist_zone(frame))
      task_param.radar.wfs(1).nyquist_zone = frames.nyquist_zone(frame);
    end
    arg{1} = task_param;
    
    if strcmp(param.sched.type,'custom_torque')
      create_task_param.conforming = true;
      create_task_param.notes = sprintf('%d/%d records %d-%d', ...
        frame, break_idx, cur_recs(1), cur_recs(end));
      ctrl = torque_create_task(ctrl,fh,1,arg,create_task_param);
      
    elseif ~strcmp(param.sched.type,'no scheduler')
      [ctrl,job_id,task_id] = create_task(ctrl,fh,1,arg);
      fprintf('  %d/%d: records %d to %d in job,task %d,%d (%s)\n', ...
        frame, break_idx, cur_recs(1), cur_recs(end), job_id, task_id, datestr(now));
      retry_fields{job_id,task_id}.frm = frame;
      retry_fields{job_id,task_id}.break_idx = break_idx;
      retry_fields{job_id,task_id}.arg = arg;
      out_recs{end + 1} = cur_recs;
      retry_fields{job_id,task_id}.out_idx = length(out_recs);
    else
      fprintf('  %d/%d: records %d to %d (%s)\n', ...
        frame, break_idx, cur_recs(1), cur_recs(end), datestr(now));
      [success] = fh(arg{1});
    end
    
  end
end

% =======================================================================
% Wait for jobs to complete if a scheduler was used
% =======================================================================
if strcmpi(param.sched.type,'custom_torque')
  % Wait until all submitted jobs to complete
  ctrl = torque_rerun(ctrl);
  if ~all(ctrl.error_mask == 0)
    if ctrl.sched.stop_on_fail
      torque_cleanup(ctrl);
      error('Not all jobs completed, but out of retries (%s)', datestr(now));
    else
      warning('Not all jobs completed, but out of retries (%s)', datestr(now));
      keyboard;
    end
  else
    fprintf('Jobs completed (%s)\n\n', datestr(now));
  end
  torque_cleanup(ctrl);
  
elseif ~strcmpi(param.sched.type,'no scheduler')
  % ======================================================================
  % Wait for jobs to finish and clean up
  ctrl.cmd = 'done';
  ctrl = create_task(ctrl);
  if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
    % Quit if a bad error occurred
    fprintf('Bad errors occurred, quitting (%s)\n\n', datestr(now));
    if strcmp(ctrl.sched.type,'torque')
      fprintf('Often on the Torque scheduler, these are not bad errors\n');
      fprintf('because of system instabilities (e.g. file IO failure)\n');
      fprintf('and the task simply needs to be resubmitted. If this is the case,\n');
      fprintf('run "ctrl.error_mask = 2" and then run "dbcont".\n');
      keyboard
      if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
        return;
      end
    else
      return
    end
  end
  
  % ======================================================================
  % Retry jobs that failed
  retry = 1;
  while ctrl.error_mask == 2 && retry <= param.sched.max_retries
    fprintf('Tasks failed, retry %d of max %d\n', retry, param.sched.max_retries);
    
    % Bookkeeping (move previous run information to "old_" variables)
    old_ctrl = ctrl;
    old_retry_fields = retry_fields;
    retry_fields = {};
    old_out_recs = out_recs;
    out_recs = {};
    
    % Initialize submission ctrl structure
    ctrl = [];
    ctrl.cmd = 'init';
    ctrl.sched = param.sched;
    ctrl.fd = fd;
    ctrl = create_task(ctrl);
    
    % Prepare submission ctrl structure for queing jobs
    ctrl.cmd = 'task';
    
    % Submit failed tasks, but keep track of these in case they fail again
    for job_idx = 1:length(old_ctrl.jobs)
      for task_idx = old_ctrl.jobs{job_idx}.error_idxs
        [ctrl,job_id,task_id] = create_task(ctrl,fh,2,old_retry_fields{job_idx,task_idx}.arg);
        out_idx = old_retry_fields{job_idx,task_idx}.out_idx;
        fprintf('  %d/%d: Processing records %d to %d in job,task %d,%d (%s)\n', ...
          old_retry_fields{job_idx,task_idx}.frm, old_retry_fields{job_idx,task_idx}.break_idx, ...
          old_out_recs{out_idx}(1), old_out_recs{out_idx}(end), ...
          job_id, task_id, datestr(now));
        retry_fields{job_id,task_id} = old_retry_fields{job_idx,task_idx};
        out_recs{end + 1} = old_out_recs{out_idx};
        retry_fields{job_id,task_id}.out_idx = length(out_recs);
      end
    end
    
    % Wait for tasks to complete and then cleanup
    ctrl.cmd = 'done';
    ctrl = create_task(ctrl);
    retry = retry + 1;
    
    if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
      % Quit if a bad error occurred
      fprintf('Bad errors occurred, quitting (%s)\n\n', datestr(now));
      if strcmp(ctrl.sched.type,'torque')
        fprintf('Often on the Torque scheduler, these are not bad errors\n');
        fprintf('because of system instabilities (e.g. file IO failure)\n');
        fprintf('and the task simply needs to be resubmitted. If this is the case,\n');
        fprintf('run "ctrl.error_mask = 2" and then run "dbcont".\n');
        keyboard
        if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
          return;
        end
      else
        return
      end
    end
  end
  if ctrl.error_mask ~= 0
    fprintf('Not all jobs completed, but out of retries (%s)\n', datestr(now));
    return;
  else
    fprintf('Jobs completed (%s)\n\n', datestr(now));
  end
end

% =====================================================================
% Loop through all the array_path directories and combine
% =====================================================================
for frame_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frame_idx);
  
  out_fn_dir = fullfile(ct_filename_out(param,param.equal.out_path,'CSARP_equal'), ...
    sprintf('equal_%03d', frm));
  
  for img = 1:length(param.equal.imgs)
    fns = get_filenames(out_fn_dir,'blk_','',sprintf('_img_%02d.mat', img));
    
    num_blocks_per_task = 1;
    layer_bins = [];
    layer_mask = [];
    layer_bins_vals = [];
    peak_val_corr = [];
    peak_offset= [];
    peak_val = [];
    roll_est_theta = [];
    roll_est_val = [];
    gps_time = [];
    roll = [];
    for block = 1:num_blocks_per_task:length(fns)
      equal_fn_name = sprintf('blk_%d_img_%02d.mat',block,img);
      equal_fn = fullfile(out_fn_dir, equal_fn_name);
      fprintf('  Combining %s (%s)\n', equal_fn_name, datestr(now));
      equal = load(equal_fn);
      
      layer_bins = [layer_bins equal.layer_bins];
      layer_mask = [layer_mask equal.layer_mask];
      layer_bins_vals = [layer_bins_vals equal.layer_bins_vals];
      peak_val_corr = [peak_val_corr equal.peak_val_corr];
      peak_offset = [peak_offset equal.peak_offset];
      peak_val = [peak_val equal.peak_val];
      roll_est_theta = [roll_est_theta equal.roll_est_theta];
      roll_est_val = [roll_est_val equal.roll_est_val];
      gps_time = [gps_time equal.gps_time];
      roll = [roll equal.roll];
      
      %delete(equal_fn);
    end
    wfs = equal.wfs;
    gps_source = equal.gps_source;
    param_equal = equal.param_equal;
    param_records = equal.param_records;
    
    out_fn = fullfile(ct_filename_out(param,param.equal.out_path,'CSARP_equal'), ...
      sprintf('Equal_%s_%03d_%02d.mat', param.day_seg, frm, img));
    fprintf('  Writing output to %s\n', out_fn);
    save(out_fn, 'layer_bins', 'layer_mask', 'layer_bins_vals','peak_val_corr','peak_offset','peak_val',...
      'roll_est_theta','roll_est_val','gps_time','gps_source','roll',...
      'param_equal', 'param_records','wfs');
  end
end

return;


for frm = 1:27
  frm
  load(sprintf('/cresis/scratch2/mdce/cr1/rds/2008_Greenland_TO/CSARP_feedthru/20080630_01/Equal_20080630_01_%03d_01.mat',frm))

  [B,A] = butter(4,0.05);
  figure(1); clf;
  h = plot(angle(filtfilt(B,A,(peak_val ./ repmat(peak_val(param_equal.equal.ref_wf_adc_idx,:),[size(peak_val,1) 1])).'))*180/pi,'.');
  legend(h,{'1','2','3','4','5','6'})
  grid on;

%   hold on;
%   h = plot(angle(filtfilt(B,A,(peak_val_corr ./ repmat(peak_val_corr(param_equal.equal.ref_wf_adc_idx,:),[size(peak_val,1) 1])).'))*180/pi,'.');
%   hold off;
%   set(h,'Marker','x');
%   legend(h,{'1','2','3','4','5','6'})
%   grid on;

  figure(2); clf;
  h = plot(medfilt1(peak_offset.',11));
  legend(h,{'1','2','3','4','5','6'})
  grid on;
  
  wf = param_equal.equal.imgs{1}(1);
  time = wfs(wf).time;
  dt = time(2)-time(1);
  
  fprintf('%.2f ',round((param_equal.radar.wfs(wf).Tsys*1e9 + mean(medfilt1(peak_offset.',11))*dt*1e9)*100)/100)
  fprintf('\n');
  
  fprintf('%.2f ', round((param_equal.radar.wfs(wf).chan_equal_deg + angle(mean(peak_val.',1))*180/pi)*100)/100)
  fprintf('\n');
  
  pause
end


