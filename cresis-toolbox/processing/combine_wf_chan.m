function ctrl_chain = combine_wf_chan(param,param_override)
% ctrl_chain = combine_wf_chan(param,param_override)
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

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks
% =====================================================================

if ~isfield(param.combine,'frm_types') || isempty(param.combine.frm_types)
  param.combine.frm_types = {-1,-1,-1,-1,-1};
end

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

if ~isfield(param.combine,'img_comb_layer_params')
  param.combine.img_comb_layer_params = [];
end

if ~isfield(param.combine,'trim_time')
  param.combine.trim_time = true;
end

if ~isfield(param.csarp,'pulse_comp') || isempty(param.csarp.pulse_comp)
  param.csarp.pulse_comp = 1;
end

if ~isfield(param.csarp,'presums') || isempty(param.csarp.presums)
  param.csarp.presums = 1;
end

if ~isfield(param.combine,'in_path') || isempty(param.combine.in_path)
  param.combine.in_path = 'out';
end

if ~isfield(param.combine,'array_path') || isempty(param.combine.array_path)
  param.combine.array_path = 'out';
end

if ~isfield(param.combine,'out_path') || isempty(param.combine.out_path)
  param.combine.out_path = param.combine.method;
end

if ~isfield(param.csarp,'out_path') || isempty(param.csarp.out_path)
  param.csarp.out_path = 'out';
end

if ~isfield(param.combine,'presums') || isempty(param.combine.presums)
  if ~isfield(param.csarp,'presums') || isempty(param.csarp.presums)
    param.combine.presums = 1;
  else
    param.combine.presums = param.csarp.presums;
  end
end

if ~isfield(param.combine,'sar_type') || isempty(param.combine.sar_type)
  if ~isfield(param.csarp,'sar_type') || isempty(param.csarp.sar_type)
    param.combine.sar_type = 'fk';
  else
    param.combine.sar_type = param.csarp.sar_type;
  end
end

if strcmpi(param.combine.sar_type,'f-k')
  error('Deprecated sar_type name. Change param.combine.sar_type from ''f-k'' to ''fk'' in  your parameters (or remove parameter since ''fk'' is the default mode).');
end

if ~isfield(param.combine,'chunk_len') || isempty(param.combine.chunk_len)
  if ~isfield(param.csarp,'chunk_len') || isempty(param.csarp.chunk_len)
    error('param.combine.chunk_len or param.csarp.chunk_len must be defined');
  else
    param.combine.chunk_len = param.csarp.chunk_len;
  end
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

if ~isfield(param.combine,'out_path') || isempty(param.combine.out_path)
  param.combine.out_path = param.combine.method;
end

%% Setup processing
% =====================================================================

% Get the standard radar name
[~,~,radar_name] = ct_output_dir(param.radar_name);

% Create output directory path
combine_out_dir = ct_filename_out(param, param.combine.out_path);

% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);
% Apply presumming
if param.csarp.presums > 1
  records.lat = fir_dec(records.lat,param.csarp.presums);
  records.lon = fir_dec(records.lon,param.csarp.presums);
  records.elev = fir_dec(records.elev,param.csarp.presums);
  records.roll = fir_dec(records.roll,param.csarp.presums);
  records.pitch = fir_dec(records.pitch,param.csarp.presums);
  records.heading = fir_dec(records.heading,param.csarp.presums);
  records.gps_time = fir_dec(records.gps_time,param.csarp.presums);
  records.surface = fir_dec(records.surface,param.csarp.presums);
end
% Along-track
along_track_approx = geodetic_to_along_track(records.lat,records.lon,records.elev);

% Preload layer for image combine if it is specified
if isempty(param.combine.img_comb_layer_params)
  layers = [];
else
  param_load_layers = param;
  param_load_layers.cmd.frms = param.cmd.frms;
  layers = opsLoadLayers(param_load_layers,param.combine.img_comb_layer_params);
end

%% Collect waveform information into one structure
%  - This is used to break the frame up into chunks
% =====================================================================
if strcmpi(radar_name,'mcrds')
  wfs = load_mcrds_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
elseif any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','mcrds','seaice','accum2'}))
  wfs = load_mcords_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
elseif any(strcmpi(radar_name,{'icards'}))% add icards---qishi
  wfs = load_icards_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5'}))
  error('Not supported');
  wfs = load_fmcw_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
  for wf=1:length(wfs)
    wfs(wf).time = param.csarp.time_of_full_support;
    wfs(wf).freq = 1;
  end
end

%% Create and setup the cluster batch
% =====================================================================
ctrl = cluster_new_batch(param);
cluster_compile({'combine_wf_chan_task.m','combine_wf_chan_combine_task.m'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);

total_num_sam = [];
if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','mcrds','seaice','accum2'}))
  for img = 1:length(param.combine.imgs)
    wf = abs(param.combine.imgs{img}{1}(1,1));
    total_num_sam(img) = wfs(wf).Nt_raw;
  end
  cpu_time_mult = 6e-8;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  for img = 1:length(param.combine.imgs)
    wf = abs(param.combine.imgs{img}{1}(1,1));
    total_num_sam(img) = 32000;
  end
  cpu_time_mult = 8e-8;
  mem_mult = 64;
  
else
  error('radar_name %s not supported yet.', radar_name);
  
end

%% Loop through all the frame directories and process the SAR chunks
% =====================================================================
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'combine_wf_chan_task';
sparam.num_args_out = 1;
prev_frm_num_chunks = [];
for frm_idx = 1:length(param.cmd.frms);
  frm = param.cmd.frms(frm_idx);
  
  if ct_proc_frame(frames.proc_mode(frm),param.combine.frm_types)
    fprintf('%s %s_%03i (%i of %i) (%s)\n', sparam.task_function, param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now));
    skip_frame = false;
  else
    fprintf('Skipping frame %s_%03i (no process frame)\n', param.day_seg, frm);
    skip_frame = true;
  end
  
  % Temporary output directory
  combine_tmp_dir = fullfile(ct_filename_out(param, param.combine.array_path), ...
    sprintf('%s_%03d', param.combine.method, frm));
  
  % Current frame goes from the start record specified in the frames file
  % to the record just before the start record of the next frame.  For
  % the last frame, the stop record is just the last record in the segment.
  start_rec = ceil(frames.frame_idxs(frm)/param.csarp.presums);
  if frm < length(frames.frame_idxs)
    stop_rec = ceil((frames.frame_idxs(frm+1)-1)/param.csarp.presums);
  else
    stop_rec = length(records.gps_time);
  end
  
  % Determine length of the frame
  frm_dist = along_track_approx(stop_rec) - along_track_approx(start_rec);
  
  % Determine number of chunks and range lines per chunk
  num_chunks = round(frm_dist / param.combine.chunk_len);
  
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
    dparam.success = '';
    for img = 1:length(param.combine.imgs)
      out_fn_name = sprintf('img_%02d_chk_%03d.mat',img,chunk_idx);
      out_fn{img} = fullfile(combine_tmp_dir,out_fn_name);
      if img == 1
        dparam.success = cat(2,dparam.success, ...
          sprintf('if ~exist(''%s'',''file'')', out_fn{img}));
      else
        dparam.success = cat(2,dparam.success, ...
          sprintf(' || ~exist(''%s'',''file'')', out_fn{img}));
      end
      if ~ctrl.cluster.rerun_only && exist(out_fn{img},'file')
        delete(out_fn{img});
      end
    end
    dparam.success = cat(2,dparam.success,sprintf('\n'));
    if 0
      % Enable this check if you want to open each output file to make
      % sure it is not corrupt.
      for img = 1:length(param.combine.imgs)
        out_fn_name = sprintf('img_%02d_chk_%03d.mat',img,chunk_idx);
        out_fn{img} = fullfile(combine_tmp_dir,out_fn_name);
        dparam.success = cat(2,dparam.success, ...
          sprintf('  load(''%s'');\n', out_fn{img}));
      end
    end
    success_error = 64;
    dparam.success = cat(2,dparam.success, ...
      sprintf('  error_mask = bitor(error_mask,%d);\n', success_error));
    dparam.success = cat(2,dparam.success,sprintf('end;\n'));
    
    % Rerun only mode: Test to see if we need to run this task
    % =================================================================
    dparam.notes = sprintf('%s:%s:%s %s_%03d (%d of %d)/%d of %d', ...
      sparam.task_function, param.radar_name, param.season_name, param.day_seg, frm, frm_idx, length(param.cmd.frms), ...
      chunk_idx, num_chunks);
    if ctrl.cluster.rerun_only
      % If we are in rerun only mode AND the get heights task success
      % condition passes without error, then we do not run the task.
      error_mask = 0;
      eval(dparam.success);
      if ~error_mask
        fprintf('  Already exists [rerun_only skipping]: %s (%s)\n', ...
          dparam.notes, datestr(now));
        continue;
      end
    end
    
    % Create task
    % =================================================================
    
    % CPU Time and Memory estimates:
    %  Nx*total_num_sam*K where K is some manually determined multiplier.
    Nx = stop_rec - start_rec + 1;
    dparam.cpu_time = 0;
    dparam.mem = 0;
    for img = 1:length(param.combine.imgs)
      dparam.cpu_time = dparam.cpu_time + 10 + size(param.combine.imgs{img},1)*Nx*total_num_sam(img)*cpu_time_mult;
      dparam.mem = max(dparam.mem,250e6 + size(param.combine.imgs{img},1)*Nx*total_num_sam(img)*mem_mult);
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
  cpu_time_mult = 6e-8;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  cpu_time_mult = 100e-8;
  mem_mult = 24;
end

sparam = [];
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'combine_wf_chan_combine_task';
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
    Nx_frm = (along_track_approx(frames.frame_idxs(frm+1)) - along_track_approx(frames.frame_idxs(frm))) / param.csarp.sigma_x;
  else
    Nx_frm = (along_track_approx(length(records.gps_time)) - along_track_approx(frames.frame_idxs(frm) + 1)) / param.csarp.sigma_x;
  end
  Nx_frm = round(Nx_frm);
  if Nx_frm > Nx_max
    Nx_max = Nx_frm;
  end
  Nx = Nx + Nx_frm;
end
% Account for averaging
for img = 1:length(param.combine.imgs)
  sparam.cpu_time = sparam.cpu_time + (Nx*total_num_sam(img)*cpu_time_mult);
  if isempty(param.combine.img_comb)
    % Individual images, so need enough memory to hold the largest image
    sparam.mem = max(sparam.mem,250e6 + Nx_max*total_num_sam(img)*mem_mult);
  else
    % Images combined into one so need enough memory to hold all images
    sparam.mem = 250e6 + Nx*sum(total_num_sam)*mem_mult;
  end
end
sparam.notes = sprintf('%s:%s:%s %s combine frames', ...
  sparam.task_function, param.radar_name, param.season_name, param.day_seg);

% Create success condition
success_error = 64;
sparam.success = '';
for frm = param.cmd.frms
  if numel(param.combine.imgs) > 1
    % More than one image: Check for individual image files
    for img = 1:length(param.combine.imgs)
      out_fn_name = sprintf('Data_img_%02d_%s_%03d.mat',img,param.day_seg,frm);
      out_fn = fullfile(combine_out_dir,out_fn_name);
      sparam.success = cat(2,sparam.success, ...
        sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
      if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
        delete(out_fn);
      end
    end
    % If combining set, then check for combined file
    if ~isempty(param.combine.img_comb)
      out_fn_name = sprintf('Data_%s_%03d.mat',param.day_seg,frm);
      out_fn = fullfile(combine_out_dir,out_fn_name);
      sparam.success = cat(2,sparam.success, ...
        sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
      if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
        delete(out_fn);
      end
    end
  else
    % Only one image: Check for combined file
    out_fn_name = sprintf('Data_%s_%03d.mat',param.day_seg,frm);
    out_fn = fullfile(combine_out_dir,out_fn_name);
    sparam.success = cat(2,sparam.success, ...
      sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
    if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
      delete(out_fn);
    end
  end
end

ctrl = cluster_new_task(ctrl,sparam,[]);

ctrl_chain{end+1} = ctrl;
    
fprintf('Done %s\n', datestr(now));

return;
