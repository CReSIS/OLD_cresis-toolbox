function ctrl_chain = collate(param, param_override)
% tomo.collate(param, param_override)
%
% Usually this function is called from tomo.run_collate.
% Calls tomo_collate_task.
%
% Inputs:
%   param = struct with processing parameters
%   param_override = parameters in this struct will override parameters
%     in param.
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfdata,
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks
% =====================================================================

% Remove frames that do not exist from param.cmd.frms list
frames = frames_load(param);
param.cmd.frms = frames_param_cmd_frms(param,frames);

% sar.* fields
% -------------------------------------------------------------------------

% param.sar gets used to create the defaults of several fields so we make
% sure the field is available
if ~isfield(param,'sar') || isempty(param.sar)
  param.sar = [];
end
if ~isfield(param.sar,'sigma_x') || isempty(param.sar.sigma_x)
  error('The param.sar.sigma_x field must be set to calculate cpu time and memory requirements.');
end

% array.* fields
% -------------------------------------------------------------------------

if ~isfield(param.array,'dline') || isempty(param.array.dline)
  error('param.array.dline must be specified.');
end

if ~isfield(param.array,'method') || isempty(param.array.method)
  param.array.method = 'music';
end

% tomo_collate.* fields
% -------------------------------------------------------------------------
if ~isfield(param.tomo_collate,'frm_types') || isempty(param.tomo_collate.frm_types)
  param.tomo_collate.frm_types = {-1,-1,-1,-1,-1};
end

if ~isfield(param.tomo_collate,'ground_based_flag') || isempty(param.tomo_collate.ground_based_flag)
  param.tomo_collate.ground_based_flag = false;
end

if ~isfield(param.tomo_collate,'in_path') || isempty(param.tomo_collate.in_path)
  param.tomo_collate.in_path = 'music3D';
end
  
if ~isfield(param.tomo_collate,'surf_out_path') || isempty(param.tomo_collate.surf_out_path)
  param.tomo_collate.surf_out_path = 'surfData';
end
  
if ~isfield(param.tomo_collate,'out_path') || isempty(param.tomo_collate.out_path)
  param.tomo_collate.out_path = param.tomo_collate.in_path;
end
  
%% Setup Processing
% =====================================================================

% Get the standard radar name
[~,~,radar_name] = ct_output_dir(param.radar_name);

% Load records file
records = records_load(param);

% Along-track
along_track_approx = geodetic_to_along_track(records.lat,records.lon,records.elev);

% Tomo collate radar echogram output directory
out_path_dir = ct_filename_out(param, param.tomo_collate.out_path);
surf_out_path_dir = ct_filename_out(param, param.tomo_collate.surf_out_path);

%% Compile C++ functions
if 0
  % If you get a C++11 option error, you may be using pre-G++ 4.7. You can
  % check the g++ version with system('g++ --version');
  % To fix this, add -v option to mex function and look for a line like this:
  %   Options file: ~/.matlab/R2015b/mex_C++_glnxa64.xml
  % Replace -std=c++11 with -std=c++0x (should occur in two places)
  % Reference: http://stackoverflow.com/questions/14674597/cc1plus-error-unrecognized-command-line-option-std-c11-with-g
  mex -largeArrayDims train_model.cpp
  mex -largeArrayDims detect.cpp
  mex -largeArrayDims extract.cpp
  mex -largeArrayDims viterbi.cpp
  mex -largeArrayDims trws.cpp
end

%% Setup cluster
% =====================================================================
ctrl = cluster_new_batch(param);
cluster_compile({'tomo_collate_task.m'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);

total_num_sam = 0;
total_img = 0;
if param.array.tomo_en
  if any(strcmpi(param.array.method,{'music'}))
    Nsv = param.array.Nsv;
  elseif any(strcmpi(param.array.method,{'mle'}))
    Nsv = 2*param.array.Nsrc;
  else
    error('Unsupported param.array.method %s\n', param.array.method);
  end
else
  error('param.array.tomo_en is false, tomography should be enabled when running tomo.collate.');
end
[wfs,~] = data_load_wfs(setfield(param,'load',struct('imgs',{param.array.imgs})),records);
if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','mcords6','mcrds','seaice','accum2','accum3'}))
  for v_img = 1:length(param.tomo_collate.imgs)
    for h_img = 1:length(param.tomo_collate.imgs{v_img})
      img = param.tomo_collate.imgs{v_img}(h_img);
      wf = param.tomo_collate.imgs{v_img}(1,1);
      if h_img == 1
        total_num_sam = total_num_sam + wfs(wf).Nt_pc;
      end
      total_img = total_img + 1;
    end
  end
  cpu_time_mult = 20e-6;
  mem_mult = 14;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  total_num_sam = 32000 * ones(size(param.tomo_collate.imgs));
  cpu_time_mult = 8e-8;
  mem_mult = 64;
  
else
  error('radar_name %s not supported yet.', radar_name);
  
end

ctrl_chain = {};

%% Create Tasks
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'tomo_collate_task';
sparam.num_args_out = 1;
sparam.argsin{1}.load.imgs = param.tomo_collate.imgs;
for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  
  % Check proc_mode from frames file that contains this frames type and
  % make sure the user has specified to process this frame type
  if ct_proc_frame(frames.proc_mode(frm),param.tomo_collate.frm_types)
    fprintf('%s %s_%03i (%i of %i) (%s)\n', sparam.task_function, param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now));
  else
    fprintf('Skipping %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end
  
  % Current frame goes from the start record specified in the frames file
  % to the record just before the start record of the next frame.  For
  % the last frame, the stop record is just the last record in the segment.
  start_rec = frames.frame_idxs(frm);
  if frm < length(frames.frame_idxs)
    stop_rec = frames.frame_idxs(frm+1)-1;
  else
    stop_rec = length(records.gps_time);
  end
  
  % Determine length of the frame
  frm_dist = along_track_approx(stop_rec) - along_track_approx(start_rec);
  
  % Prepare task inputs
  % =================================================================
  dparam = [];
  dparam.argsin{1}.load.frm = frm;
  
  % Create success condition
  % =================================================================
  dparam.file_success = {};
  if param.tomo_collate.fuse_images_flag || param.tomo_collate.add_icemask_surfacedem_flag
    out_fn = fullfile(out_path_dir, sprintf('Data_%s_%03d.mat', ...
      param.day_seg, frm));
    dparam.file_success{end+1} = out_fn;
  end
  if param.tomo_collate.create_surfData_flag
    out_fn = fullfile(surf_out_path_dir, sprintf('Data_%s_%03d.mat', ...
      param.day_seg, frm));
    dparam.file_success{end+1} = out_fn;
  end
  if ~ctrl.cluster.rerun_only
    % Mark file for deletion
    ct_file_lock_check(out_fn,3);
  end
  
  % Rerun only mode: Test to see if we need to run this task
  % =================================================================
  dparam.notes = sprintf('%s:%s:%s:%s %s_%03d (%d of %d)', ...
    sparam.task_function, param.radar_name, param.season_name, out_path_dir, param.day_seg, frm, frm_idx, length(param.cmd.frms));
  if ctrl.cluster.rerun_only
    % If we are in rerun only mode AND the tomo_collate task file success
    % condition passes without error then we do not run the task.
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
  Nx = round(frm_dist / param.sar.sigma_x / param.array.dline);
  dparam.cpu_time = 10 + Nx*Nsv*total_num_sam*cpu_time_mult;
  dparam.mem = 250e6 + Nx*Nsv*total_num_sam*mem_mult;
  ctrl = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);
end

ctrl = cluster_save_dparam(ctrl);

ctrl_chain{end+1} = ctrl;

fprintf('Done %s\n', datestr(now));
