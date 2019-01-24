function ctrl_chain = array(param,param_override)
% ctrl_chain = array(param,param_override)
%
% This script combines the receive channels and outputs the result
% for each waveform. It also combines the waveforms. It takes in
% sar files one directory at a time and:
%  1. combines the receive channels
%  2. concatenates the results
%  3. square-law detects the data, abs()^2
%  4. takes incoherent averages (multilooks data)
%  5. saves the result in a new directory
%
% The assumption is that the directories in the input_path are named
% using the following convention:
%   PROC-TYPE-STRING_data_#{_SUBAPERTURE-STRING}
% where
%   PROC-TYPE-STRING can be 'fk' or 'tdbp'for f-k migrated or time domain
%   back projected respectively.
%   _data_ is always present
%   #, \d+: one or more numbers
%   _SUBAPERTURE-STRING, {_[mp]\d\.\d}: optional subaperture string
% Examples:
%   fk_data_001_01_01: f-k migrated, frame 1, subaperture 1, subband 1
%   fk_data_004_02_01: f-k migrated, frame 4, subaperture 2, subband 1
%   fk_data_001_03_01: f-k migrated, frame 1, subaperture 3, subband 1
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
% See also: run_master.m, master.m, run_array.m, array.m, load_sar_data.m,
% array_proc.m, array_task.m, array_combine_task.m


%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks
% =====================================================================

% Remove frames that do not exist from param.cmd.frms list
load(ct_filename_support(param,'','frames')); % Load "frames" variable
if ~isfield(param.cmd,'frms') || isempty(param.cmd.frms)
  param.cmd.frms = 1:length(frames.frame_idxs);
end
[valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
if length(valid_frms) ~= length(param.cmd.frms)
  bad_mask = ones(size(param.cmd.frms));
  bad_mask(keep_idxs) = 0;
  warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
    param.cmd.frms(find(bad_mask,1)));
  param.cmd.frms = valid_frms;
end

% Handles multilooking syntax:
%  {{[1 1],[1 2],[1 3],[1 4],[1 5]},{[2 1],[2 2],[2 3],[2 4],[2 5]}}
%  If the image is a cell array it describes multilooking across apertures
if ~iscell(param.array.imgs{1})
  % No special multilooking, reformat old syntax to new multilooking syntax
  for img = 1:length(param.array.imgs)
    param.array.imgs{img} = {param.array.imgs{img}};
  end
end

% param.sar gets used to create the defaults of several fields so we make
% sure the field is available
if ~isfield(param,'sar') || isempty(param.sar)
  param.sar = [];
end
if ~isfield(param.sar,'sigma_x') || isempty(param.sar.sigma_x)
  error('The param.sar.sigma_x field must be set to calculate cpu time and memory requirements.');
end

% param.array.* fields
% -------------------------------------------------------------------------

if ~isfield(param.array,'chunk_len') || isempty(param.array.chunk_len)
  if ~isfield(param.sar,'chunk_len') || isempty(param.sar.chunk_len)
    error('param.array.chunk_len or param.sar.chunk_len must be defined');
  else
    param.array.chunk_len = param.sar.chunk_len;
  end
end

if ~isfield(param.array,'frm_types') || isempty(param.array.frm_types)
  param.array.frm_types = {-1,-1,-1,-1,-1};
end

if ~isfield(param.array,'in_path') || isempty(param.array.in_path)
  param.array.in_path = 'sar';
end

if ~isfield(param.array,'out_path') || isempty(param.array.out_path)
  param.array.out_path = param.array.method;
end

if ~isfield(param.array,'presums') || isempty(param.array.presums)
  if ~isfield(param.sar,'presums') || isempty(param.sar.presums)
    param.sar.presums = 1;
  end
  if ~isfield(param.sar,'presums') || isempty(param.sar.presums)
    param.array.presums = 1;
  else
    param.array.presums = param.sar.presums;
  end
end

if ~isfield(param.array,'sar_type') || isempty(param.array.sar_type)
  if ~isfield(param.sar,'sar_type') || isempty(param.sar.sar_type)
    param.array.sar_type = 'fk';
  else
    param.array.sar_type = param.sar.sar_type;
  end
end
if strcmpi(param.array.sar_type,'f-k')
  error('Deprecated sar_type name. Change param.array.sar_type from ''f-k'' to ''fk'' in  your parameters (or remove parameter since ''fk'' is the default mode).');
end

if ~isfield(param.array,'subaps') || isempty(param.array.subaps)
  % If SAR sub-apertures not set, we assume that there is just one
  % subaperture to be passed in for each multilook input
  for img = 1:length(param.array.imgs)
    for ml_idx = 1:length(param.array.imgs{img})
      param.array.subaps{img}{ml_idx} = [1];
    end
  end
end

if ~isfield(param.array,'subbnds') || isempty(param.array.subbnds)
  % If subbands not set, we assume that there is just one
  % subaperture to be passed in for each multilook input
  for img = 1:length(param.array.imgs)
    for ml_idx = 1:length(param.array.imgs{img})
      param.array.subbnds{img}{ml_idx} = [1];
    end
  end
end

if ~isfield(param.array,'surf_layer') || isempty(param.array.surf_layer)
  param.array.surf_layer.name = 'surface';
  param.array.surf_layer.source = 'layerData';
end
% Never check for the existence of layers
param.array.surf_layer.existence_check = false;

% param.array.* fields used by array_proc.m
% -------------------------------------------------------------------------
param = array_proc(param);

%% Setup processing
% =====================================================================

% Get the standard radar name
[~,~,radar_name] = ct_output_dir(param.radar_name);

% Create output directory path
array_out_dir = ct_filename_out(param, param.array.out_path);

% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);
% Apply presumming
if param.sar.presums > 1
  records.lat = fir_dec(records.lat,param.sar.presums);
  records.lon = fir_dec(records.lon,param.sar.presums);
  records.elev = fir_dec(records.elev,param.sar.presums);
  records.roll = fir_dec(records.roll,param.sar.presums);
  records.pitch = fir_dec(records.pitch,param.sar.presums);
  records.heading = fir_dec(records.heading,param.sar.presums);
  records.gps_time = fir_dec(records.gps_time,param.sar.presums);
  records.surface = fir_dec(records.surface,param.sar.presums);
end
% Along-track
along_track_approx = geodetic_to_along_track(records.lat,records.lon,records.elev);

%% Collect waveform information into one structure
%  - This is used to break the frame up into chunks
% =====================================================================
[wfs,~] = data_load_wfs(setfield(param,'load',struct('imgs',{param.sar.imgs})),records);
param.radar.wfs = merge_structs(param.radar.wfs,wfs);

%% Create and setup the cluster batch
% =====================================================================
ctrl = cluster_new_batch(param);
cluster_compile({'array_task.m','array_combine_task.m'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);

total_num_sam = zeros(size(param.array.imgs));
total_num_sam_combine = zeros(size(param.array.imgs));
if param.array.three_dim.en
  Nsv = param.array.Nsv;
else
  Nsv = 1;
end
if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','mcrds','seaice','accum2'}))
  for img = 1:length(param.array.imgs)
    wf = param.array.imgs{img}{1}(1,1);
    total_num_sam_combine(img) = total_num_sam_combine(img) + wfs(wf).Nt_pc/param.array.dbin*Nsv;
    for ml_idx = 1:length(param.array.imgs{img})
      wf = param.array.imgs{img}{ml_idx}(1,1);
      total_num_sam(img) = total_num_sam(img) + wfs(wf).Nt_pc/param.array.dbin ...
        * size(param.array.imgs{img}{ml_idx},1) * numel(param.array.subaps{img}{ml_idx});
    end
  end
  cpu_time_mult = 3e-6;
  mem_mult = 16;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  for img = 1:length(param.array.imgs)
    wf = param.array.imgs{img}{1}(1,1);
    total_num_sam(img) = 32000;
  end
  cpu_time_mult = 8e-8;
  mem_mult = 64;
  
else
  error('radar_name %s not supported yet.', radar_name);
  
end

array_proc_methods; % This script assigns the integer values for each method
switch (param.array.method)
  case STANDARD_METHOD
  case MVDR_METHOD
    cpu_time_mult = cpu_time_mult*4;
  case MUSIC_METHOD
    cpu_time_mult = cpu_time_mult*4;
  case MLE_METHOD
    cpu_time_mult = cpu_time_mult*480;
end

%% Loop through all the frame directories and process the SAR chunks
% =====================================================================
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'array_task';
sparam.num_args_out = 1;
prev_frm_num_chunks = [];
for frm_idx = 1:length(param.cmd.frms);
  frm = param.cmd.frms(frm_idx);
  
  if ct_proc_frame(frames.proc_mode(frm),param.array.frm_types)
    fprintf('%s %s_%03i (%i of %i) (%s)\n', sparam.task_function, param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now));
    skip_frame = false;
  else
    fprintf('Skipping frame %s_%03i (no process frame)\n', param.day_seg, frm);
    skip_frame = true;
  end
  
  % Create combine_file_success for this frame only which is used for
  % rerun_only==true checks
  if ctrl.cluster.rerun_only
    combine_file_success = {};
    if length(param.array.imgs) > 1
      for img = 1:length(param.array.imgs)
        out_fn = fullfile(array_out_dir, sprintf('Data_img_%02d_%s_%03d.mat', ...
          img, param.day_seg, frm));
        combine_file_success{end+1} = out_fn;
      end
    end
    if length(param.array.imgs) == 1 || ~isempty(param.array.img_comb)
      % A combined file should be created
      out_fn = fullfile(array_out_dir, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
      combine_file_success{end+1} = out_fn;
    end
  end
  
  % Temporary output directory
  array_tmp_dir = fullfile(ct_filename_out(param, param.array.out_path, 'array_tmp'), ...
    sprintf('array_%03d', frm));
  
  % Current frame goes from the start record specified in the frames file
  % to the record just before the start record of the next frame.  For
  % the last frame, the stop record is just the last record in the segment.
  start_rec = ceil(frames.frame_idxs(frm)/param.sar.presums);
  if frm < length(frames.frame_idxs)
    stop_rec = ceil((frames.frame_idxs(frm+1)-1)/param.sar.presums);
  else
    stop_rec = length(records.gps_time);
  end
  
  % Determine length of the frame
  frm_dist = along_track_approx(stop_rec) - along_track_approx(start_rec);
  
  % Determine number of chunks and range lines per chunk
  num_chunks = round(frm_dist / param.array.chunk_len);
  
  %% Process each chunk (unless it is a skip frame)
  for chunk_idx = 1:num_chunks*~skip_frame
    % Prepare task inputs
    % =================================================================
    dparam = [];
    dparam.argsin{1}.load.frm = frm;
    dparam.argsin{1}.load.chunk_idx = chunk_idx;
    dparam.argsin{1}.load.num_chunks = num_chunks;
    dparam.argsin{1}.load.prev_frm_num_chunks = prev_frm_num_chunks;
  
    % Create success condition
    % =================================================================
    dparam.file_success = {};
    for img = 1:length(param.array.imgs)
      out_fn_name = sprintf('img_%02d_chk_%03d.mat',img,chunk_idx);
      out_fn = fullfile(array_tmp_dir,out_fn_name);
      dparam.file_success{end+1} = out_fn;
      if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
        delete(out_fn);
      end
    end
    
    % Rerun only mode: Test to see if we need to run this task
    % =================================================================
    dparam.notes = sprintf('%s:%s:%s %s_%03d (%d of %d)/%d of %d', ...
      sparam.task_function, param.radar_name, param.season_name, param.day_seg, frm, frm_idx, length(param.cmd.frms), ...
      chunk_idx, num_chunks);
    if ctrl.cluster.rerun_only
      % If we are in rerun only mode AND the array task file success
      % condition passes without error (or the combined file passes),
      % then we do not run the task.
      if ~cluster_file_success(dparam.file_success) || ~cluster_file_success(combine_file_success)
        fprintf('  Already exists [rerun_only skipping]: %s (%s)\n', ...
          dparam.notes, datestr(now));
        continue;
      end
    end
    
    % Create task
    % =================================================================
    
    % CPU Time and Memory estimates:
    %  Nx*total_num_sam*K where K is some manually determined multiplier.
    Nx = round(frm_dist / num_chunks / param.sar.sigma_x);
    dparam.cpu_time = 0;
    dparam.mem = 0;
    for img = 1:length(param.array.imgs)
      dparam.cpu_time = dparam.cpu_time + 10 + Nx*total_num_sam(img)*cpu_time_mult;
      % Take the max of the input data size and the output data size
      dparam.mem = max(dparam.mem,250e6 + max( Nx*total_num_sam(img)*mem_mult, ...
        Nx/param.array.dline*total_num_sam_combine(img)*mem_mult ));
    end
    
    ctrl = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);
    
  end
  prev_frm_num_chunks = num_chunks;
end

ctrl = cluster_save_dparam(ctrl);

ctrl_chain = {ctrl};

%% Create and setup the combine batch
% =====================================================================
ctrl = cluster_new_batch(param);

if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','mcrds','seaice','accum2'}))
  cpu_time_mult = 1e-6;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  cpu_time_mult = 100e-8;
  mem_mult = 24;
end

sparam = [];
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'array_combine_task';
sparam.num_args_out = 1;
sparam.cpu_time = 60;
sparam.mem = 0;
fprintf('%s %s (%s)\n', sparam.task_function, param.day_seg, datestr(now));
% Add up all records being processed and find the most records in a frame
Nx = 0;
Nx_max = 0;
for frm = param.cmd.frms
  % recs: Determine the records for this frame
  if frm < length(frames.frame_idxs)
    Nx_frm = (along_track_approx(frames.frame_idxs(frm+1)) - along_track_approx(frames.frame_idxs(frm))) / param.sar.sigma_x / param.array.dline;
  else
    Nx_frm = (along_track_approx(length(records.gps_time)) - along_track_approx(frames.frame_idxs(frm) + 1)) / param.sar.sigma_x / param.array.dline;
  end
  Nx_frm = round(Nx_frm);
  if Nx_frm > Nx_max
    Nx_max = Nx_frm;
  end
  Nx = Nx + Nx_frm;
end
% Account for averaging
for img = 1:length(param.array.imgs)
  sparam.cpu_time = sparam.cpu_time + (Nx*total_num_sam_combine(img)*cpu_time_mult) * (1 + (Nsv-1)*0.2);
  if isempty(param.array.img_comb)
    % Individual images, so need enough memory to hold the largest image
    sparam.mem = max(sparam.mem,250e6 + Nx_max*total_num_sam_combine(img)*mem_mult);
  else
    % Images combined into one so need enough memory to hold all images
    sparam.mem = max(sparam.mem,250e6 + Nx_max*sum(total_num_sam_combine)*mem_mult);
  end
end
sparam.notes = sprintf('%s:%s:%s %s combine frames', ...
  sparam.task_function, param.radar_name, param.season_name, param.day_seg);

% Create success condition
sparam.file_success = {};
for frm = param.cmd.frms
  if length(param.array.imgs) > 1
    for img = 1:length(param.array.imgs)
      out_fn = fullfile(array_out_dir, sprintf('Data_img_%02d_%s_%03d.mat', ...
        img, param.day_seg, frm));
      sparam.file_success{end+1} = out_fn;
      if ~ctrl.cluster.rerun_only
        % Mark file for deletion
        ct_file_lock_check(out_fn,3);
      end
    end
  end
  if length(param.array.imgs) == 1 || ~isempty(param.array.img_comb)
    % A combined file should be created
    out_fn = fullfile(array_out_dir, sprintf('Data_%s_%03d.mat', ...
      param.day_seg, frm));
    sparam.file_success{end+1} = out_fn;
    if ~ctrl.cluster.rerun_only
      % Mark file for deletion
      ct_file_lock_check(out_fn,3);
    end
  end

end

ctrl = cluster_new_task(ctrl,sparam,[]);

ctrl_chain{end+1} = ctrl;
    
fprintf('Done %s\n', datestr(now));

return;
