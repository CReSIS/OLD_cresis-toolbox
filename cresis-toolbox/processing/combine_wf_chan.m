function combine_wf_chan(param,param_override)
% combine_wf_chan(param,param_override)
%
% This script combines the receive channels and outputs the result
% for each waveform. It also combines the waveforms. It takes in
% f-k files one directory at a time and:
%  1. combines the receive channels
%  2. concatenates the results
%  3. square-law detects the data, abs()^2
%  4. takes incoherent averages (multilooks data)
%  5. saves the result in a new directory
%  6. Combines the waveforms
%
% The assumption is that the directories in the input_path are named
% using the following convention:
%   PROC-TYPE-STRING_data_#{_SUBAPERTURE-STRING}
% where
%   PROC-TYPE-STRING can be 'fk','tdbp', or 'pc' for f-k migrated,time domain
%   back projected,and pulse compressed respectively ('fk' and tdbp supported)
%   _data_ is always present
%   #, \d+: one or more numbers
%   _SUBAPERTURE-STRING, {_[mp]\d\.\d}: optional subaperture string
% Examples:
%   fk_data_01_01: f-k migrated, frame 1, subaperture 1
%   fk_data_04_02: f-k migrated, frame 4, subaperture 2
%   fk_data_01_03: f-k migrated, frame 1, subaperture 3
%   pc_data_01: pulse compressed only, frame 1
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
% See also: run_master.m, master.m, run_combine_wf_chan.m, combine_wf_chan.m,
%   combine_wf_chan_task.m

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

if ~isfield(param,'debug_level')
  param.debug_level = 1;
end

if ~isfield(param.sched,'rerun_only') || isempty(param.sched.rerun_only)
  param.sched.rerun_only = false;
end

% Handles multilooking syntax:
%  {{[1 1],[1 2],[1 3],[1 4],[1 5]},{[2 1],[2 2],[2 3],[2 4],[2 5]}}
%  If the image is a cell array it describes multilooking across apertures
if ~iscell(param.combine.imgs{1})
  % No special multilooking, reformat old syntax to new multilooking syntax
  for img = 1:length(param.combine.imgs)
    param.combine.imgs{img} = {param.combine.imgs{img}};
  end
end

for img = 1:length(param.combine.imgs)
  for ml_idx = 1:length(param.combine.imgs{img})
    % Imaginary image indices is for IQ combining during raw data load
    % which we do not need here.
    param.combine.imgs{img}{ml_idx} = abs(param.combine.imgs{img}{ml_idx});
  end
end

img_list = param.combine.imgs;
  
in_path = ct_filename_out(param, ...
  param.combine.in_path, 'CSARP_out');

array_path = ct_filename_out(param, ...
  param.combine.array_path, 'CSARP_out');

out_path = ct_filename_out(param, ...
  param.combine.out_path, sprintf('CSARP_%s', ...
  param.combine.method));

% Create the output directory
if ~exist(out_path,'dir')
  mkdir(out_path);
end

% Load frames file
if ~isfield(param.records,'frames_fn')
  param.records.frames_fn = '';
end
load(ct_filename_support(param,param.records.frames_fn,'frames'));

% Check frames to process variable
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

% =====================================================================
% Setup the scheduler
% =====================================================================
fd = [get_filenames(param.path,'','','.m',struct('recursive',1)); ...
  get_filenames(param.path,'','','.mexa64',struct('recursive',1))];
fd_override = [get_filenames(param.path_override,'','','.m',struct('recursive',1)); ...
  get_filenames(param.path_override,'','','.mexa64',struct('recursive',1))];

fd = merge_filelists(fd, fd_override);

if strcmpi(param.sched.type,'custom_torque')
  global ctrl; % Make this global for convenience in debugging
  ctrl = torque_new_batch(param);
  fprintf('Torque batch: %s\n', ctrl.batch_dir);
  torque_compile('combine_wf_chan_task.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
elseif ~strcmpi(param.sched.type,'no scheduler')
  % Initialize submission ctrl structure
  global ctrl
  ctrl = [];
  ctrl.cmd = 'init';
  ctrl.sched = param.sched;
  ctrl.fd = fd;
  ctrl = create_task(ctrl);
  param.surf.manual = 0; % Turn manual pick off
  
  % Prepare submission ctrl structure for queing jobs
  ctrl.cmd = 'task';
end

%% Loop through all the frame directories and process the fk
% chunks in those directories
% =====================================================================
retry_fields = {};
for frm_idx = 1:length(param.cmd.frms);
  frm = param.cmd.frms(frm_idx);
  
  if ct_proc_frame(frames.proc_mode(frm),param.csarp.frm_types)
    fprintf('%s combine %s_%03i (%i of %i) %s\n', param.radar_name, param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now,'HH:MM:SS'));
  else
    fprintf('Skipping frame %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end
  
  %% Input directory for this frame (only look at the first subaperture
  % "_01_01" since combine_wf_chan_task will know which subapertures to load)
  if strcmpi(param.csarp.sar_type,'f-k')
    sar_type = 'fk';
  elseif strcmpi(param.csarp.sar_type,'tdbp')
    sar_type = 'tdbp';
  elseif strcmpi(param.csarp.sar_type,'mltdp')
    sar_type = 'mltdp';
  end
  param.combine.in_path = fullfile(in_path,sprintf('%s_data_%03d_01_01', sar_type, frm));
  
  %% Output directory
  param.combine.out_path = fullfile(array_path,sprintf('array_%03d', frm));
  
  %% Get the filenames for each chunk of data processed by csarp
  % DEBUG:
  %  - Get all the chunk files
  %  - Normal operation is 1 to +inf (i.e. all files)
  %  - The start/stop funtionality is not used except for debugging
  %  - Set start_chunk and stop_chunk to restrict which chunks get processed
  start_chunk = 1;
  stop_chunk = inf;
  % Get all time stamps in the directory: the assumption is that csarp.m
  % has created all the necessary files for each wf/adc pair required.  So
  % to get the time stamps, we just search for all the files for a particular
  % pair (in this case the first one in the list).
  img = 1;
  wf_adc_idx = 1;
  filenames = get_filenames(param.combine.in_path,'', ...
    sprintf('wf_%02d_adc_%02d',img_list{img}{1}(wf_adc_idx,1), ...
    img_list{img}{1}(wf_adc_idx,2)),'.mat');
  if isempty(filenames)
    error('No filenames found in %s', param.combine.in_path);
  end
  chunk_ids = {};
  for idx = 1:length(filenames)
    [path,name,ext] = fileparts(filenames{idx});
    chunk_idx = str2double(name(end-2:end));
    if chunk_idx >= start_chunk && chunk_idx <=stop_chunk
      chunk_ids{end+1} = name(end-2:end);
    end
  end
  
  %% Create and clean the array_proc output directories
  if exist(param.combine.out_path,'dir') && ~param.sched.rerun_only
    % If folders do exist, clear them out
    fprintf('  Cleaning array directory %s\n', param.combine.out_path);
    rmdir(param.combine.out_path,'s');
  end
  mkdir(param.combine.out_path);
  
  %% Combine Channels: Standard beam-forming, MUSIC, MVDR, etc)
  % - This is setup so that multiple fk chunks can be processed in the
  %   same task/job.
  
  load(filenames{1},'param_csarp');
  if strcmpi(param_csarp.csarp.sar_type,'f-k')
    if max(param.combine.rline_rng) - min(param.combine.rline_rng) > param_csarp.csarp.chunk_overlap
      error('SAR processing chunks will not align properly, chunk_overlap too small');
    end
  end
  param.combine.sar_type = param_csarp.csarp.sar_type;
  num_chunks_per_task = 1;
    
  for chunk_idx = 1:num_chunks_per_task:length(chunk_ids)
    %% To make the SAR processed chunks fit together seamlessly without
    % having to resample, we determine the start range line output for
    % each chunk.
    % chunk_idxs: SAR chunks that will be processed by this task,
    %   plus one additional one for calculating rlines(2)
    chunk_idxs = chunk_idx + (0:num_chunks_per_task);
    % chunk_Nx: the number of non-overlapping SAR chunk outputs
    chunk_Nx = floor(param_csarp.csarp.chunk_len / param_csarp.csarp.sigma_x);
    % min_offset: the minimum offset into the SAR chunk which array_proc can
    %   output a full support estimate (since the output uses a neighborhood
    %   of points around the pixel in question, the first output line generally
    %   be from the first input line)
    min_offset = -min(param.combine.rline_rng);
    % rlines(1,:): this will be the first range line output by array_proc
    %   for each SAR chunk this task is array processing
    rlines = [];
    rlines(1,:) = 1+mod(min_offset+param.combine.dline-(chunk_idxs-1)*chunk_Nx, param.combine.dline);
    rlines(rlines<min_offset) = rlines(rlines<min_offset) + ceil(param.combine.dline/(1+min_offset)) * param.combine.dline;
    rlines(2,1:end-1) = chunk_Nx + rlines(1,2:end) - param.combine.dline;
    rlines = rlines(:,1:end-1);
    
    % Check if this is the last chunk. This last chunk could have variable
    % length and we want to return all of the data from this chunk. To tell
    % combine_task to do this, we set rlines(2) to infinity for
    % this chunk
    if chunk_idx+num_chunks_per_task-1 >= length(chunk_ids)
      rlines(2,end) = inf;
    end
    
    %% Get the chunk ids that this task will process
    chunk_idx_last = min(chunk_idx+num_chunks_per_task-1, length(chunk_ids));
    param.combine.chunk_ids = chunk_ids(chunk_idx:chunk_idx_last);
    param.combine.rlines = rlines;
    
    %% Pass in the frame
    param.combine.frm = frm;
    
    %% Rerun only mode: Test to see if we need to run this task
    if param.sched.rerun_only
      % If we are in rerun only mode AND all the combine_wf_chan task output files
      % already exists, then we do not run the task
      file_exists = true;
      for img = 1:length(param.combine.imgs)
%         array_fn = fullfile(param.combine.out_path, ...
%           sprintf('chk_%03d_img_%02d.mat', chunk_ids{chunk_idx}, img));
        array_fn = fullfile(param.combine.out_path, ...
          sprintf('chk_%s_img_%02d.mat', chunk_ids{chunk_idx}, img));
        if ~exist(array_fn,'file')
          file_exists = false;
        end
      end
      if file_exists
        fprintf('  %s already exists [rerun_only skipping] (%s)\n', ...
          param.combine.chunk_ids{1}, datestr(now));
        continue;
      end
    end
    
    %% Execute tasks/jobs
    fh = @combine_wf_chan_task;
    arg{1} = param;
    
    if strcmp(param.sched.type,'custom_torque')
      create_task_param.conforming = true;
      create_task_param.notes = sprintf('%s', ...
        param.combine.chunk_ids{1});
      ctrl = torque_create_task(ctrl,fh,1,arg,create_task_param);
      
    elseif ~strcmp(param.sched.type,'no scheduler')
      [ctrl,job_id,task_id] = create_task(ctrl,fh,1,arg);
      fprintf('  %s task %d,%d (%s)\n', ...
        param.combine.chunk_ids{1}, job_id, task_id, datestr(now));
      retry_fields{job_id,task_id}.file_time = param.combine.chunk_ids{1};
      retry_fields{job_id,task_id}.arg = arg;
      retry_fields{job_id,task_id}.frm_dir_name = sprintf('Frm %d', frm);
    else
      success = fh(arg{1});
    end
  end
end

%% Wait for jobs to complete if a scheduler was used
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
        fprintf('  %s/%s task %d,%d (%s)\n', ...
          old_retry_fields{job_idx,task_idx}.frm_dir_name, ...
          old_retry_fields{job_idx,task_idx}.file_time, job_id, task_id, datestr(now));
        retry_fields{job_id,task_id} = old_retry_fields{job_idx,task_idx};
      end
    end
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

%% Loop through all the frames
% =====================================================================
for frm_idx = 1:length(param.cmd.frms);
  frm = param.cmd.frms(frm_idx);
  
  if ct_proc_frame(frames.proc_mode(frm),param.csarp.frm_types)
    fprintf('combine %s_%03i (%i of %i) %s\n', param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now,'HH:MM:SS'));
  else
    fprintf('Skipping frame %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end
  
  %% Output directory
  param.combine.out_path = fullfile(array_path,sprintf('array_%03d', frm));
  
  %% Loop through all the images
  for img = 1:length(param.combine.imgs)
    
    %% Loop through all the chunks and combine
    Latitude = [];
    Longitude = [];
    Elevation = [];
    Roll = [];
    Pitch = [];
    Heading = [];
    GPS_time = [];
    Surface = [];
    Bottom = [];
    Data = [];
    Topography = [];
    chunk_fns = get_filenames(param.combine.out_path,'chk','',sprintf('img_%02d.mat',img));
    for chunk_idxs = 1:length(chunk_fns)
      tmp = load(chunk_fns{chunk_idxs});
      Time = tmp.Time;
      Latitude = [Latitude double(tmp.Latitude)];
      Longitude = [Longitude double(tmp.Longitude)];
      Elevation = [Elevation double(tmp.Elevation)];
      Roll = [Roll double(tmp.Roll)];
      Pitch = [Pitch double(tmp.Pitch)];
      Heading = [Heading double(tmp.Heading)];
      GPS_time = [GPS_time tmp.GPS_time];
      Surface = [Surface double(tmp.Surface)];
      Bottom = [Bottom double(tmp.Bottom)];
      Data = [Data tmp.Data];
      param_records = tmp.param_records;
      param_csarp = tmp.param_csarp;
      if chunk_idxs == 1
      param_combine = tmp.param_combine;
        param_combine.array_param.fcs{1}{1}.x = tmp.param_combine.array_param.fcs{1}{1}.x(:,tmp.param_combine.array_param.lines);
        param_combine.array_param.fcs{1}{1}.y = tmp.param_combine.array_param.fcs{1}{1}.y(:,tmp.param_combine.array_param.lines);
        param_combine.array_param.fcs{1}{1}.z = tmp.param_combine.array_param.fcs{1}{1}.z(:,tmp.param_combine.array_param.lines);
        param_combine.array_param.fcs{1}{1}.origin = tmp.param_combine.array_param.fcs{1}{1}.origin(:,tmp.param_combine.array_param.lines);
      else
        % Concatenate the fcs field
        param_combine.array_param.fcs{1}{1}.x = [param_combine.array_param.fcs{1}{1}.x tmp.param_combine.array_param.fcs{1}{1}.x(:,tmp.param_combine.array_param.lines)];
        param_combine.array_param.fcs{1}{1}.y = [param_combine.array_param.fcs{1}{1}.y tmp.param_combine.array_param.fcs{1}{1}.y(:,tmp.param_combine.array_param.lines)];
        param_combine.array_param.fcs{1}{1}.z = [param_combine.array_param.fcs{1}{1}.z tmp.param_combine.array_param.fcs{1}{1}.z(:,tmp.param_combine.array_param.lines)];
        param_combine.array_param.fcs{1}{1}.origin = [param_combine.array_param.fcs{1}{1}.origin tmp.param_combine.array_param.fcs{1}{1}.origin(:,tmp.param_combine.array_param.lines)];
      end
      if isfield(tmp,'Topography')
%         3D-surface is present so concatenate it too
%         Topography = cat(3,Topography,tmp.Topography);
%         Concatenate all the fields under struct Topography: valR, bins, val, freq
%         and img.
        fields = fieldnames(tmp.Topography);      
        if chunk_idxs == 1
          for field_idx = 1:length(fields)
            Topography.(fields{field_idx}) = tmp.Topography.(fields{field_idx});
          end       
        else        
          for field_idx = 1:length(fields)
            max_dim = length(size(tmp.Topography.(fields{field_idx})));
            Topography.(fields{field_idx}) = cat(max_dim,Topography.(fields{field_idx}),tmp.Topography.(fields{field_idx}));
          end       
        end
        
      end
    end
    
    % =====================================================================
    % Save output
    if length(param.combine.imgs) == 1
      out_fn = fullfile(out_path, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
    else
      out_fn = fullfile(out_path, sprintf('Data_img_%02d_%s_%03d.mat', ...
        img, param.day_seg, frm));
    end
    fprintf('  Writing output to %s\n', out_fn);
    if isempty(Topography)
      % Do not save 3D surface
      save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
        'Elevation','GPS_time','Data','Surface','Bottom', ...
        'param_combine','param_records','param_csarp', ...
        'Roll', 'Pitch', 'Heading');
    else
      % Save 3D surface
      save('-v7.3',out_fn,'Topography','Time','Latitude', ...
        'Longitude','Elevation','GPS_time','Data','Surface','Bottom', ...
        'param_combine','param_records','param_csarp', ...
        'Roll', 'Pitch', 'Heading');
    end
  end
  
  if isempty(param.combine.img_comb)
    % No image combining is required
    continue;
  end
  
  if length(param.combine.img_comb) ~= 3*(length(param.combine.imgs)-1)
    warning('param.combine.img_comb not the right length. There should be 3 entries for each image combination interface ([Tpd second image for surface saturation, -inf for second image blank, Tpd first image to avoid roll off] is typical). Set correctly here and update param spreadsheet before dbcont.');
    keyboard
  end
  
  %% Load each image and then combine with previous image (also trim time<0 values)
  for img = 1:length(param.combine.imgs)
    
    if length(param.combine.imgs) == 1
      out_fn = fullfile(out_path, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
    else
      out_fn = fullfile(out_path, sprintf('Data_img_%02d_%s_%03d.mat', ...
        img, param.day_seg, frm));
    end
    if img == 1
      load(out_fn);
      first_idx = find(Time <= 0,1,'last');
      if ~isempty(first_idx)
        Time = Time(first_idx:end);
        Data = Data(first_idx:end,:);
      end
    else
      append = load(out_fn,'Time','Data');
      %% Combine images
      % Data,Time => already loaded data
      % append.Data, append.Time => new data to append
      % New_Time, New_Data => Combined result
      
      % Interpolate image N onto already loaded data (assumption is that image
      % N-1 always comes before image N)
      dt = Time(2)-Time(1);
      New_Time = (Time(1) : dt : append.Time(end)).';
      append.Data = interp1(append.Time,append.Data,New_Time,'linear',0);
      
      % Surface tracking image combine
      %  param.combine.img_comb(1): time after surface return where
      %    combine will happen
      %  param.combine.img_comb(2): minimum time that combine will occur
      %  param.combine.img_comb(3): guard time which specifies how
      %    many seconds at the end of img1 will not be used... this is
      %    important because the last samples of img1 will have low signal
      %    power and blurred because they will only have captured a portion
      %    of the chirp energy (typically this will be set to something
      %    close to the pulse duration for img1)
      %  param.combine.img_comb(4-6, 7-9, etc.): same fields as above
      %    except between images 2 and 3, 3 and 4, etc.
      
      Surface = interp_finite(Surface,0);
      % First row of img_bins indicates the start of the blend-region
      img_bins = round(interp1(New_Time, 1:length(New_Time), ...
        max(Surface+param.combine.img_comb((img-2)*3+1),param.combine.img_comb((img-2)*3+2)), 'linear','extrap'));
      
      % Determine guard at end of image 1 that will not be used
      guard_bins = 1 + round(param.combine.img_comb((img-2)*3+3)/dt);
      
      % Check to make sure requested time is inside window and just
      % force the combination bin to occur at the second to last bin
      %   img_bins outside the img1 time window will be NaN due to interp1
      %   img_bins inside the img1 time window may still be larger than
      %     the guard allows
      max_good_time = length(Time)*ones(1,size(Data,2));
      invalid_rlines = find(isnan(img_bins) ...
        | img_bins > max_good_time-guard_bins);
      img_bins(invalid_rlines) = max_good_time(invalid_rlines)-guard_bins;
      
      % Second row of img_bins indicates the end of the blend-region
      img_bins(2,:) = img_bins(1,:) + 1;
      
      difference = 10^(-0/10);
      
      % Combine images
      New_Data = zeros(size(append.Data),'single');
      for rline = 1:size(New_Data,2)
        trans_bins = img_bins(1,rline)+1:img_bins(2,rline);
        weights = 0.5+0.5*cos(pi*linspace(0,1,length(trans_bins)).');
        if trans_bins <= size(append.Data,1)
          New_Data(:,rline) = [Data(1:img_bins(1,rline),rline); ...
            weights.*Data(trans_bins,rline) ...
            + difference*(1-weights).*append.Data(trans_bins,rline); ...
            difference*append.Data(img_bins(2,rline)+1:end,rline)];
        else
          New_Data(:,rline) = Data(1:size(New_Data,1),rline);
        end        
      end
      Time = New_Time;
      Data = New_Data;
    end
  end
  
  %% Save output
  out_fn = fullfile(out_path, sprintf('Data_%s_%03d.mat', ...
    param.day_seg, frm));
  fprintf('  Writing output to %s\n', out_fn);
  save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
    'Elevation','GPS_time','Data','Surface','Bottom', ...
    'param_combine','param_records','param_csarp', ...
    'Roll', 'Pitch', 'Heading');
end

fprintf('Done %s\n', datestr(now));

return;
