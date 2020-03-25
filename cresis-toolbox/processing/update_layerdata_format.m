function ctrl_chain = update_layerdata_format(param,param_override)
% update_layerdata_format(param,param_override)
%
% This function updates layer data files of an old version (0) to a new version (1), and save
% the updated layer data files to a layer data destination (default is CSARP_layer).
%
% param = struct with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Example:
%  See run_update_layerdata_format.m for how to run this function directly.
%  This function may be called from master.m using the param spreadsheet and 
%  the cmd.generic column.
%
% Authors: Jilu Li
%
% See also: run_master.m, master.m, run_update_layerdata_format.m, update_layerdata_formatk.m,
% update_layerdata_format_task.m

%% General Setup
% =====================================================================

param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input checks
% =====================================================================
if ~isfield(param,'frame_overlap_removal') || isempty(param.frame_overlap_removal)
  param.frame_overlap_removal = true;
end

if ~isfield(param,'layerdata_source') || isempty(param.layerdata_source)
  param.layerdata_source = 'layerData';
end

if ~isfield(param,'file_version_old') || isempty(param.file_version_old)
  param.file_version_old = '0';
end

if ~isfield(param,'layerdata_dest') || isempty(param.layerdata_dest)
  param.layerdata_dest = 'layer';
end

if ~isfield(param,'file_version_new') || isempty(param.file_version_new)
  param.file_version_new = '1';
end


% Load frames file
load(ct_filename_support(param,'','frames'));

% If no frames specified, then do all frames
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

% Set up output path
output_path = ct_filename_out(param,param.layerdata_dest,'');
if ~exist(output_path,'dir')
  mkdir(output_path);
end

%% Setup Processing
% =====================================================================

% Get the standard radar name
[~,~,radar_name] = ct_output_dir(param.radar_name);

% Load records file
records_fn = ct_filename_support(param,'','records');
if ~exist(records_fn)
  error('You must run create the records file before running anything else:\n  %s', records_fn);
end
records = load(records_fn);

out_fn_dir = ct_filename_out(param, param.layerdata_dest);

%% Create and setup the cluster batch
% =====================================================================
ctrl = cluster_new_batch(param);
cluster_compile({'update_layerdata_format_task.m'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);
cpu_time_mult = 8e-8;
mem_mult = 64;

ctrl_chain = {};

%% Create task
% =====================================================================
%
% For each frame load REC_BLOCK_SIZE records at a time (code groups
% by file index, but has to watch negative offset values which imply
% the record starts in a previous file and carries over into the next)
%    --> The last block can range from 0.5 to 1.5 * REC_BLOCK_SIZE
% =====================================================================

% Set the task's static parameters
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'update_layerdata_format_task';
sparam.num_args_out = 1;
for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  
  % rerun_only==true checks
  if ctrl.cluster.rerun_only
    update_file_success = {};
    out_fn = fullfile(out_fn_dir, sprintf('Data_%s_%03d.mat',param.day_seg, frm));
    update_file_success{end+1} = out_fn;
  end
  
  % recs: Determine the records for this frame
  if frm < length(frames.frame_idxs)
    recs = frames.frame_idxs(frm):frames.frame_idxs(frm+1)-1;
    recs = [frames.frame_idxs(frm),frames.frame_idxs(frm+1)-1];
  else
    recs = frames.frame_idxs(frm):length(records.gps_time);
    recs = [frames.frame_idxs(frm),length(records.gps_time)];
  end
    
  % Set the task's dynamic parameters
  % =================================================================
  dparam = [];
  dparam.argsin{1}.load.frm = frm;
  dparam.argsin{1}.load.recs = recs;
  
  % Create success condition
  % =================================================================
  dparam.file_success = {};
  out_fn_name = sprintf('Data_%s_%03d.mat', param.day_seg, frm);
  out_fn = fullfile(out_fn_dir,out_fn_name);
  dparam.file_success{end+1} = out_fn;
  if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
    delete(out_fn);
  end

  % Rerun only mode: Test to see if we need to run this task
  % =================================================================
  dparam.notes = sprintf('%s:%s:%s:%s %s_%03d (%d of %d)/%d of %d', ...
    sparam.task_function, param.radar_name, param.season_name, param.layerdata_dest, param.day_seg, frm, frm_idx, length(param.cmd.frms));
  if ctrl.cluster.rerun_only
    % If we are in rerun only mode AND the layer data file success
    % condition passes without error then we do not run the task.
    if ~cluster_file_success(dparam.file_success) || ~cluster_file_success(update_file_success)
      fprintf('  Already exists [rerun_only skipping]: %s (%s)\n', ...
        dparam.notes, datestr(now));
      continue;
    end
  end
    
  % Create a cluster task for each frame
  % =================================================================
  
  % CPU Time and Memory estimates:
  %  Nx*total_num_sam*K where K is some manually determined multiplier.
  Nx = recs(end)-recs(1)+1;
  dparam.cpu_time = 0;
  dparam.mem = 250e6;
  dparam.cpu_time = dparam.cpu_time + 10 + Nx*1000*log2(1000)*cpu_time_mult;
  dparam.mem = dparam.mem + Nx*1000*mem_mult;
  ctrl = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);
end

ctrl = cluster_save_dparam(ctrl);

ctrl_chain{end+1} = ctrl;

