function get_heights(param,param_override)
% get_heights(param,param_override)
%
% This function has two capabilities (enabled separately)
% 1. Automated or manual pick: Gets the air/ice surface height using the
%    data like an altimeter (manual mode is optional)
% 2. Generate quick look outputs (default location: CSARP_qlook)
%
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Example:
%  See run_get_heights.m for how to run this function directly.
%  Normally this function is called from master.m using the param spreadsheet.
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_get_heights.m, get_heights.m,
%   get_heights_task.m

% =====================================================================
% General Setup
% =====================================================================

if ~isstruct(param)
  % Functional form
  param();
end
param = merge_structs(param, param_override);

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

% =====================================================================
% Setup processing
% =====================================================================

% Get WGS84 ellipsoid parameters
physical_constants;

if ~isfield(param.get_heights,'combine_only') || isempty(param.get_heights.combine_only)
  param.get_heights.combine_only = false;
end

% Load frames file
load(ct_filename_support(param,'','frames'));
% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);

global g_data;
g_data = [];

qlook_out_path = ct_filename_out(param, param.get_heights.qlook.out_path, 'CSARP_qlook');

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

if ~isfield(param.get_heights,'ground_based')
  param.get_heights.ground_based = [];
end

if ~isfield(param.get_heights.qlook,'save_format') || isempty(param.get_heights.qlook.save_format)
  param.get_heights.qlook.save_format = '7.3';
end
save_format = sprintf('-v%s',param.get_heights.qlook.save_format);

if isfield(param.get_heights,'deconvolution') ...
    && ~isempty(param.get_heights.deconvolution) ...
    && param.get_heights.deconvolution == 3
  %% Get version information out of the deconvolution file
  out_fn_dir = ct_filename_out(param,'', 'CSARP_noise');
  out_segment_fn_dir = fileparts(out_fn_dir);
  out_segment_fn = fullfile(out_segment_fn_dir,sprintf('deconv_%s.mat', param.day_seg));
  spec = load(out_segment_fn,'param_collate');
  
  param.get_heights.deconvolution_sw_version = spec.param_collate.sw_version;
  param.get_heights.deconvolution_params = spec.param_collate.analysis.specular;
end

if isfield(param.get_heights,'coh_noise_method') ...
    && ~isempty(param.get_heights.coh_noise_method) ...
    && any(param.get_heights.coh_noise_method == [17 19])
  %% Get version information out of the coherent noise file
  
  cdf_fn_dir = fileparts(ct_filename_out(param,param.get_heights.coh_noise_arg{4}, ''));
  cdf_fn = fullfile(cdf_fn_dir,sprintf('coh_noise_simp_%s.nc', param.day_seg));
  
  tmp = netcdf_to_mat(cdf_fn,[],'^sw_version.*');
  param.get_heights.coh_noise_version = tmp.sw_version;
  tmp = netcdf_to_mat(cdf_fn,[],'^param_collate.*');
  param.get_heights.coh_noise_params = tmp.param_collate;
end

% Cleanup folders
if ~param.sched.rerun_only
  if exist(qlook_out_path,'dir')
    for frm = param.cmd.frms
      del_paths = get_filenames(qlook_out_path,sprintf('ql_data_%03d',frm),'','',struct('type','d'));
      for idx = 1:length(del_paths)
        fprintf('If required, manually remove path: %s\n', del_paths{idx});
        %rmdir(del_paths{idx},'s');
      end
    end
  end
end

% =====================================================================
% Setup static inputs for get_heights_task
% =====================================================================
task_param = param;
task_param.load.imgs = param.get_heights.imgs;
if isempty(param.get_heights.imgs)
  error('No images specified in param.get_heights.imgs.');
end

if ~param.get_heights.combine_only
  % =====================================================================
  % Setup the scheduler
  % =====================================================================
  
  if strcmpi(param.sched.type,'custom_torque')
    global ctrl; % Make this global for convenience in debugging
    ctrl = torque_new_batch(param);
    fprintf('Torque batch: %s\n', ctrl.batch_dir);
    torque_compile('get_heights_task.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
    
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
  
  if ~strcmpi(param.sched.type,'no scheduler')
    param.get_heights.surf.manual = 0; % Turn manual pick off
  end
  
  % =====================================================================
  % Load data and run get_heights tasks
  %
  % For each frame load REC_BLOCK_SIZE records at a time (code groups
  % by file index, but has to watch negative offset values which imply
  % the record starts in a previous file and carries over into the next)
  %    --> The last block can be up to 2*REC_BLOCK_SIZE
  % =====================================================================
  out_recs = {};
  retry_fields = {};
  for frm_idx = 1:length(param.cmd.frms)
    frm = param.cmd.frms(frm_idx);
    
    % Check digits of proc_mode from frames file and make sure the user has
    % specified to process this frame type
    if ct_proc_frame(frames.proc_mode(frm),param.get_heights.frm_types)
      fprintf('get_heights %s_%03i (%i of %i) %s\n', param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now,'HH:MM:SS'));
    else
      fprintf('Skipping %s_%03i (no process frame)\n', param.day_seg, frm);
      continue;
    end
    
    ql_path = fullfile(qlook_out_path, sprintf('ql_data_%03d_01_01',frm));
    % Create the array_proc output directories
    if ~exist(ql_path,'dir')
      mkdir(ql_path);
    end
    
    if frm < length(frames.frame_idxs)
      recs = frames.frame_idxs(frm):frames.frame_idxs(frm+1)-1;
    else
      recs = frames.frame_idxs(frm):length(records.lat);
    end
    
    % Determine where breaks in processing are going to occur
    rec = recs(1);
    if isempty(param.get_heights.block_size)
      REC_BLOCK_SIZE = 10000;
    else
      if numel(param.get_heights.block_size) == 2
        REC_BLOCK_SIZE = param.get_heights.block_size(1);
        block_overlap = param.get_heights.block_size(2);
      else
        REC_BLOCK_SIZE = param.get_heights.block_size;
        block_overlap = 0;
      end
    end
    if length(recs) < 2*REC_BLOCK_SIZE
      breaks = 1;
    else
      breaks = 1:REC_BLOCK_SIZE:length(recs)-REC_BLOCK_SIZE;
    end
    
    task_param.proc.frm = frm;
    
    % Begin loading data
    for break_idx = 1:length(breaks)
      % Determine the current records being processed
      if break_idx < length(breaks)
        cur_recs_keep = [recs(breaks(break_idx)) recs(breaks(break_idx+1)-1)];
        cur_recs = [max(1,recs(breaks(break_idx))-block_overlap) ...
          recs(breaks(break_idx+1)-1)+block_overlap];
      else
        cur_recs_keep = [recs(breaks(break_idx)) recs(end)];
        cur_recs = [max(1,recs(breaks(break_idx))-block_overlap) min(length(records.lat),recs(end)+block_overlap)];
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
          param.get_heights.qlook.out_path, 'CSARP_qlook'), ...
          sprintf('ql_data_%03d_%02d_%02d',frm, ...
          sub_apt_shift_idx, sub_band_idx));
        start_time_for_fn = records.gps_time(cur_recs(1));
        for imgs_list_idx = 1:length(param.get_heights.imgs)
          out_fn = sprintf('%s_img_%02d', ...
            datestr(epoch_to_datenum(start_time_for_fn),'yyyymmdd_HHMMSS'), ...
            imgs_list_idx);
          out_full_fn = fullfile(out_path,[out_fn '.mat']);
          lower_out_fn = sprintf('%s_img_%02d', ...
            datestr(epoch_to_datenum(start_time_for_fn-1),'yyyymmdd_HHMMSS'), ...
            imgs_list_idx);
          lower_out_full_fn = fullfile(out_path,[lower_out_fn '.mat']);
          upper_out_fn = sprintf('%s_img_%02d', ...
            datestr(epoch_to_datenum(start_time_for_fn+1),'yyyymmdd_HHMMSS'), ...
            imgs_list_idx);
          upper_out_full_fn = fullfile(out_path,[upper_out_fn '.mat']);
          if ~exist(out_full_fn,'file') && ~exist(lower_out_full_fn,'file') && ~exist(upper_out_full_fn,'file')
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
      fh = @get_heights_task;
      if isfield(frames,'nyquist_zone') && ~isnan(frames.nyquist_zone(frm))
        task_param.radar.wfs(1).nyquist_zone = frames.nyquist_zone(frm);
      elseif isfield(param.radar.wfs(1),'nyquist_zone') ...
          && ~isempty(param.radar.wfs(1).nyquist_zone) ...
          && ~isnan(param.radar.wfs(1).nyquist_zone)
        task_param.radar.wfs(1).nyquist_zone = param.radar.wfs(1).nyquist_zone;
      end
      arg{1} = task_param;
      
      if strcmp(param.sched.type,'custom_torque')
        create_task_param.conforming = true;
        create_task_param.notes = sprintf('%s %s_%03d (%d of %d)/%d of %d records %d-%d', ...
          param.radar_name, param.day_seg, frm, frm_idx, length(param.cmd.frms), break_idx, length(breaks), cur_recs(1), cur_recs(end));
        ctrl = torque_create_task(ctrl,fh,1,arg,create_task_param);
        
      elseif strcmp(param.sched.type,'ollie')
        dynamic_param.frms.(['frm',num2str(frm)]).frm_id = frm;
        dynamic_param.frms.(['frm',num2str(frm)]).breaks.(['break',num2str(break_idx)]).break_id = break_idx;
        dynamic_param.frms.(['frm',num2str(frm)]).breaks.(['break',num2str(break_idx)]).recs = task_param.load.recs;
        dynamic_param.frms.(['frm',num2str(frm)]).breaks.(['break',num2str(break_idx)]).recs_keep = task_param.load.recs_keep;
        
      elseif ~strcmp(param.sched.type,'no scheduler')
        [ctrl,job_id,task_id] = create_task(ctrl,fh,1,arg);
        fprintf('  %d/%d: records %d to %d in job,task %d,%d (%s)\n', ...
          frm, break_idx, cur_recs(1), cur_recs(end), job_id, task_id, datestr(now));
        retry_fields{job_id,task_id}.frm = frm;
        retry_fields{job_id,task_id}.break_idx = break_idx;
        retry_fields{job_id,task_id}.arg = arg;
        out_recs{end + 1} = cur_recs;
        retry_fields{job_id,task_id}.out_idx = length(out_recs);
      else
        fprintf('  %s_%03d (%d of %d)/%d of %d: records %d-%d (%s)\n', ...
          param.day_seg, frm, frm_idx, length(param.cmd.frms), break_idx, length(breaks), cur_recs(1), cur_recs(end), datestr(now));
        [success] = fh(arg{1});
      end
    end
  end
  
  % Export parameter structs in case of Schedule Type Ollie
  if strcmp(param.sched.type,'ollie')
    dynamic_param.day_seg = param.day_seg;
    steady_param = task_param;
    dynamic_param_file_name = sprintf('%s/qlook_%s_dynamic_param.mat', param.slurm_jobs_path, param.day_seg);
    save(dynamic_param_file_name,'dynamic_param');
    steady_param_file_name = sprintf('%s/qlook_%s_steady_param.mat', param.slurm_jobs_path, param.day_seg);
    save(steady_param_file_name,'steady_param');
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
    get_heights_check_cluster_files; % Function call temporarily added to track down compute system problem
    torque_cleanup(ctrl);
    
  elseif strcmp(param.sched.type,'ollie')
    return
    
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
  
end

%% Check img_comb
if length(param.get_heights.imgs) == 1 || isempty(param.get_heights.qlook.img_comb)
  num_imgs = 1;
else
  num_imgs = length(param.get_heights.imgs);
  if length(param.get_heights.qlook.img_comb) ~= 3*(num_imgs-1)
    warning('param.get_heights.qlook.img_comb not the right length. Since it is not empty, there should be 3 entries for each image combination interface ([Tpd second image for surface saturation, -inf for second image blank, Tpd first image to avoid roll off] is typical). Set correctly here and update param spreadsheet before dbcont.');
    keyboard
  end
end

% =====================================================================
%% Loop through all the frames: combine and surface track
[output_dir,radar_type] = ct_output_dir(param.radar_name);
for frm_idx = 1:length(param.cmd.frms);
  frm = param.cmd.frms(frm_idx);
  
  % Check digits of proc_mode from frames file and make sure the user has
  % specified to process this frame type
  if ct_proc_frame(frames.proc_mode(frm),param.get_heights.frm_types)
    fprintf('get_heights combine frame %s_%03i (%i of %i) %s\n', param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now));
  else
    fprintf('Skipping frame %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end
  
  %% Output directory
  in_path = fullfile(qlook_out_path, sprintf('ql_data_%03d_01_01',frm));
  
  %% Concatenate blocks for each of the images
  for img = 1:length(param.get_heights.imgs)
    
    filenames = get_filenames(in_path,'','',sprintf('img_%02d.mat', img));
    
    if length(param.get_heights.imgs) == 1
      out_fn = fullfile(qlook_out_path, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
    else
      out_fn = fullfile(qlook_out_path, sprintf('Data_img_%02d_%s_%03d.mat', ...
        img, param.day_seg, frm));
    end
    
    %% Concatenate all the ql block outputs for this image
    Time = [];
    Latitude = [];
    Longitude = [];
    Elevation = [];
    Roll = [];
    Pitch = [];
    Heading = [];
    GPS_time = [];
    Data = [];
    custom = [];
    if ~param.get_heights.combine_only
      for block_idx = 1:length(filenames)
        qlook_fn = filenames{block_idx};
        tmp = load(qlook_fn);
        time_vector_changed = false;
        if isempty(Time)
          Time = tmp.Time;
        elseif any(size(Time) ~= size(tmp.Time)) || any(Time ~= tmp.Time)
          % Determine the new time axis
          time_vector_changed = true;
          old_time = Time;
          dt = Time(2) - Time(1);
          start_time_diff = round((Time(1) - tmp.Time(1))/dt);
          end_time_diff = round((tmp.Time(end) - Time(end))/dt);
          if start_time_diff > 0
            Time = [Time(1)+dt*(-start_time_diff:-1)'; Time];
          end
          if end_time_diff > 0
            Time = [Time; Time(end)+dt*(1:end_time_diff)'];
          end
        end
        Latitude = [Latitude double(tmp.Latitude)];
        Longitude = [Longitude double(tmp.Longitude)];
        Elevation = [Elevation double(tmp.Elevation)];
        Roll = [Roll double(tmp.Roll)];
        Pitch = [Pitch double(tmp.Pitch)];
        Heading = [Heading double(tmp.Heading)];
        GPS_time = [GPS_time tmp.GPS_time];
        param_records = tmp.param_records;
        param_get_heights = tmp.param_get_heights;
        
        if isfield(tmp,'custom') && ~isempty(tmp.custom)
          % Custom fields are present, so concatenate them on the second dimension
          % - Used for deconvolution waveform index
          fields = fieldnames(tmp.custom);
          if block_idx == 1
            for field_idx = 1:length(fields)
              custom.(fields{field_idx}) = tmp.custom.(fields{field_idx});
            end
          else
            for field_idx = 1:length(fields)
              max_dim = length(size(tmp.custom.(fields{field_idx})));
              custom.(fields{field_idx}) = cat(max_dim,custom.(fields{field_idx}),tmp.custom.(fields{field_idx}));
            end
          end
          
        end
        
        if time_vector_changed
          if strcmpi(radar_type,'fmcw')
            Data = [interp1(old_time,Data,Time,'nearest',0) interp1(tmp.Time,tmp.Data,Time,'nearest',0)];
          else
            Data = [interp1(old_time,Data,Time,'linear',0) interp1(tmp.Time,tmp.Data,Time,'linear',0)];
          end
        else
          Data = [Data tmp.Data];
        end
      end
      %% Save output
      fprintf('  Writing output to %s\n', out_fn);
      Data = single(Data);
      if isempty(custom)
        save(save_format,out_fn,'Time','Latitude','Longitude', ...
          'Elevation','Roll','Pitch','Heading','GPS_time','Data', ...
          'param_get_heights','param_records');
      else
        save(save_format,out_fn,'Time','Latitude','Longitude', ...
          'Elevation','Roll','Pitch','Heading','GPS_time','Data', ...
          'param_get_heights','param_records','custom');
      end
    else
      fprintf('  Reading output %s\n', out_fn);
      load(out_fn);
    end
    
    %% Create temporary output for surface tracker
    if img == 1
      Time_Surface = Time;
      Data_Surface = Data;
    elseif ~isempty(param.get_heights.qlook.img_comb)
      %% Combine image with previous
      
      if Time(end) > Time_Surface(end)
        combine.idx                   = img;
        combine.Data                  = Data_Surface;
        combine.Time                  = Time_Surface;
        combine.appendData            = Data;
        combine.appendTime            = Time;
        combine.imb_comb_surf         = 0;
        combine.img_comb_weights      = [];
        combine.img_comb_mult         = inf;
        combine.img_comb              = param.get_heights.qlook.img_comb;
        combine.img_comb_bins         = 0;
        combine.img_comb_weights_mode = 'get_heights';
        
        % Call img_combine
        [Data, Time]                  = img_combine(combine);
      end
    end
  end
  
  %% Run surface tracker
  surf = param.get_heights.surf;
  if isfield(param.get_heights.surf,'min_bin')
    % Convert time min_bin into range bins
    surf.min_bin = find(Time > param.get_heights.surf.min_bin, 1);
  end
  if isfield(param.get_heights.surf,'max_bin') && ~isempty(param.get_heights.surf.max_bin)
    % Convert time max_bin into range bins
    surf.max_bin = find(Time > param.get_heights.surf.max_bin, 1);
  end
  if isfield(param.get_heights.surf,'max_diff')
    % Convert time max_diff into range bins
    dt = Time(2) - Time(1);
    surf.max_diff = param.get_heights.surf.max_diff/dt;
  end
  
  if ~isfield(surf,'manual')
    surf.manual = false;
  end
  
  
  if isfield(surf,'feedthru')
    %% Optional feed through removal
    
    % Interpolate feed through power levels on to data time axis
    feedthru_threshold = interp1(surf.feedthru.time,surf.feedthru.power_dB,Time);
    feedthru_threshold = interp_finite(feedthru_threshold);
    
    % Remove all data not exceeding feed through threshold power
    for rline=1:size(Data,2)
      Data(:,rline) = Data(:,rline) .* (lp(Data(:,rline)) > feedthru_threshold);
    end
  end
  
  if ~isempty(param.get_heights.ground_based)
    % Hack for ground based radar, surface time is zero
    Surface = param.get_heights.ground_based * ones(1,size(Data,2));
    
  else
    if surf.manual
      [new_surface,pnt] = tracker_snake_simple(Data,surf);
      fprintf('  Press F1 for help\n');
      layer = tracker_snake_manual_gui(lp(Data),pnt);
      
    elseif strcmpi(surf.method,'threshold')
      new_surface = tracker_threshold(Data,surf);
    elseif strcmpi(surf.method,'max')
      new_surface = tracker_max(Data,surf);
    elseif strcmpi(surf.method,'snake')
      new_surface = tracker_snake_simple(Data,surf);
    else
      error('Not a supported surface tracking method.');
    end
    
    %% Apply optional median filter
    if isfield(surf,'medfilt') && ~isempty(surf.medfilt)
      new_surface = medfilt1(new_surface,surf.medfilt);
    end
    
    %% Convert from range bins to two way travel time
    Surface = interp1(1:length(Time), Time, new_surface);
    
    Surface = reshape(Surface, [1 length(Surface)]);
  end
  
  % Reset the "Data" variable in case it was modified during surface
  % tracking
  Data = Data_Surface;
  
  %% Combine images into a single image (also trim time<0 values)
  % Load each image and then combine with previous image
  for img = 1:num_imgs
    if length(param.get_heights.imgs) == 1
      out_fn = fullfile(qlook_out_path, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
    else
      out_fn = fullfile(qlook_out_path, sprintf('Data_img_%02d_%s_%03d.mat', ...
        img, param.day_seg, frm));
    end
    if img == 1
      load(out_fn);
      first_idx = find(Time <= 0,1,'last');
      if ~isempty(first_idx)
        Time = Time(first_idx:end);
        Data = Data(first_idx:end,:,:);
      end
    else
      append = load(out_fn,'Time','Data');
      combine.idx                   = img;
      combine.Data                  = Data;
      combine.Time                  = Time;
      combine.appendData            = append.Data;
      combine.appendTime            = append.Time;
      combine.imb_comb_surf         = Surface;
      combine.img_comb_weights      = [];
      combine.img_comb_mult         = inf;
      combine.img_comb              = param.get_heights.qlook.img_comb;
      combine.img_comb_bins         = 0;
      combine.img_comb_weights_mode = 'get_heights';
      
      % Call img_combine
      [Data, Time]                  = img_combine(combine);
    end
  end
  
  %% Save combined image output
  out_fn = fullfile(qlook_out_path, sprintf('Data_%s_%03d.mat', ...
    param.day_seg, frm));
  fprintf('  Writing output to %s\n', out_fn);
  Data = single(Data);
  if isempty(custom)
    save(save_format,out_fn,'Time','Latitude','Longitude', ...
      'Elevation','Roll','Pitch','Heading','GPS_time','Data','Surface', ...
      'param_get_heights','param_records');
  else
    save(save_format,out_fn,'Time','Latitude','Longitude', ...
      'Elevation','Roll','Pitch','Heading','GPS_time','Data','Surface', ...
      'param_get_heights','param_records','custom');
  end
  
  %% Remove the temporary block files
  if ~param.sched.rerun_only && exist(in_path,'dir')
    rmdir(in_path,'s');
  end
end

if param.get_heights.surf.en
  % Read the "Surface" variable from all the frames that were created
  % by this particular run of get_heights
  
  if isfield(records,'lat')
    num_recs = length(records.lat);
  end
  
  if ~isfield(records,'surface')
    records.surface = zeros(1,num_recs);
  end
  
  if length(records.surface) ~= num_recs
    warning('Debug catch incase records.surface is wrong length... hand-fix to be correct length');
    keyboard
  end
  
  for frm = param.cmd.frms
    if ~ct_proc_frame(frames.proc_mode(frm),param.get_heights.frm_types)
      continue;
    end
    out_fn = fullfile(qlook_out_path, sprintf('Data_%s_%03d.mat', ...
      param.day_seg, frm));
    load(out_fn,'GPS_time','Surface');
    
    if frm < length(frames.frame_idxs)
      recs = frames.frame_idxs(frm):frames.frame_idxs(frm+1);
    else
      recs = frames.frame_idxs(frm):num_recs;
    end
    
    records.surface(recs) = interp1(GPS_time,Surface,records.gps_time(recs),'linear','extrap');
  end
  
  % Store surface information to the records file
  save(records_fn,'-APPEND','-struct','records','surface');
  create_records_aux_files(records_fn);
end

fprintf('Done %s\n', datestr(now));

return;

