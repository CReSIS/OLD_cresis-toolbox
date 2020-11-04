function [ctrl_chain,param] = layer_tracker(param,param_override)
% [ctrl_chain,param] = layer_tracker(param,param_override)
%
% Check input parameters and create tracking tasks for running on a cluster
% with layer_tracker. See run_layer_tracker for an example of routine
% tracking of echograms. See run_layer_tracker_tune for an example of
% hyperparameter tuning to improve tracking parameters. The function
% "layer_tracker_task" does the actual tracking and
% "layer_tracker_combine_task" combines the tracking results and stores them
% into standard layer storage locations (either the OPS database or
% layerdata files).
%
% First stage temporary outputs stored in:
% /cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_layer_tracker_tmp/CSARP_layer/20140313_08/
% Second stage output stored in any format supported by opsCopyLayers and
% input parameters control where the final combiend output goes.
%
% Comparing four different methods for example might store the outputs like this:
%   layer_tracker_001/t001_lsm.mat, ..., layer_tracker_00N/t001_lsm.mat
%   layer_tracker_001/t002_mcmc.mat, ..., layer_tracker_00N/t002_mcmc.mat
%   layer_tracker_001/t003_stereo.mat, ..., layer_tracker_00N/t003_stereo.mat
%   layer_tracker_001/t004_viterbi.mat, ..., layer_tracker_00N/t004_viterbi.mat
% Layers in the files (all combined into one file during combine):
%   t001_lsm_surface_001, ..., t001_lsm_surface_016, t001_lsm_bottom_001, ..., t001_lsm_bottom_016
%   t002_mcmc_surface, t002_mcmc_bottom
%   t003_stereo_surface, t003_stereo_bottom
%   t004_viterbi_bottom
%
% Comparing the same method with four different sets of parameters:
%   layer_tracker_001/t001_lsm.mat, ..., layer_tracker_00N/t001_lsm.mat
%   layer_tracker_001/t002_lsm.mat, ..., layer_tracker_00N/t001_lsm.mat
%   layer_tracker_001/t003_lsm.mat, ..., layer_tracker_00N/t001_lsm.mat
%   layer_tracker_001/t004_lsm.mat, ..., layer_tracker_00N/t001_lsm.mat
% Layers in the files (all combined into one file during combine):
%   t001_lsm_surface_001, ..., t001_lsm_surface_016, t001_lsm_bottom_001, ..., t001_lsm_bottom_016
%   t002_lsm_surface_001, ..., t002_lsm_surface_016, t002_lsm_bottom_001, ..., t002_lsm_bottom_016
%   t003_lsm_surface_001, ..., t003_lsm_surface_016, t003_lsm_bottom_001, ..., t003_lsm_bottom_016
%   t004_lsm_surface_001, ..., t004_lsm_surface_016, t004_lsm_bottom_001, ..., t004_lsm_bottom_016
%
% Authors: Anjali Pare, John Paden
%
% See also: layer_tracker.m, layer_tracker_combine_task.m,
% layer_tracker_task.m, layer_tracker_profile.m, run_layer_tracker.m,
% run_layer_tracker_tune.m

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks: cmd
% =====================================================================

% Remove frames that do not exist from param.cmd.frms list
frames = frames_load(param);
param.cmd.frms = frames_param_cmd_frms(param,frames);

%% Input Checks: layer_tracker
% =====================================================================

layer_tracker_input_check;

%% Set up Cluster
% ===================================================================

ctrl = cluster_new_batch(param);
%cluster_compile({'layer_tracker_task'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);
cluster_compile({'layer_tracker_task','layer_tracker_combine_task'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);
ctrl_chain = {};

%% layer_tracker
% =========================================================================
% =========================================================================

% Cluster setup
% -------------------------------------------------------------------------

sparam.argsin{1} = param;
sparam.task_function = 'layer_tracker_task';
sparam.num_args_out = 1;
sparam.argsin{1}.load.echogram_img = param.layer_tracker.echogram_img;
sparam.cpu_time = 60;
sparam.mem = 500e6;
sparam.notes = '';

cpu_time_mult = zeros(size(param.layer_tracker.track));
mem_mult = zeros(size(param.layer_tracker.track));
for track_idx = 1:length(param.layer_tracker.track)
  switch param.layer_tracker.track{track_idx}.method
    case 'viterbi'
      cpu_time_mult(track_idx) = 11e-6;
      mem_mult(track_idx) = 64;
      
    case 'lsm'
      cpu_time_mult(track_idx) = 5.5e-7*max(param.layer_tracker.track{track_idx}.lsm.storeIter);
      mem_mult(track_idx) = 80;
      
    otherwise
      cpu_time_mult(track_idx) = 11e-6;
      mem_mult(track_idx) = 64;
  end
end

%% layer_tracker: Loop to create tasks
% -------------------------------------------------------------------------
in_fn_dir = ct_filename_out(param,param.layer_tracker.echogram_source,'');
if strcmp(param.layer_tracker.layer_params.source,'ops')
  tmp_out_fn_dir_dir = ct_filename_out(param,'ops','layer_tracker_tmp');
  param.layer_tracker.layer_params.layerdata_source = 'ops'; % Only used for stdout
else
  tmp_out_fn_dir_dir = ct_filename_out(param,param.layer_tracker.layer_params.layerdata_source,'layer_tracker_tmp');
end
mem_combine = 0;
cputime_combine = 0;
frm_idx = 1;
while frm_idx <= length(param.cmd.frms)
  Nx = 0;
  Nt = 0;
  
  start_frm_idx = frm_idx;
  frms = [];
  for subblock_idx = 1:param.layer_tracker.block_size_frms
    if frm_idx > param.cmd.frms
      break;
    end
    frm = param.cmd.frms(frm_idx);
    if ~any(frm == param.cmd.frms)
      break;
    end
    % Add frame to this block
    frm_idx = frm_idx + 1;
    frms(end+1) = frm;
    
    % Compute matrix size
    % ---------------------------------------------------------------------
    if param.layer_tracker.echogram_img == 0
      data_fn = fullfile(in_fn_dir,sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    else
      data_fn = fullfile(in_fn_dir,sprintf('Data_img_%02d_%s_%03d.mat',param.layer_tracker.echogram_img,param.day_seg,frm));
    end
    try
      mdata = load(data_fn, 'GPS_time','Time');
      if (subblock_idx==1)
        max_time = mdata.Time(end);
        min_time = mdata.Time(1);
      else
        if(max_time <= mdata.Time(end))
          max_time = mdata.Time(end);
        end
        if(min_time >= mdata.Time(1))
          min_time = mdata.Time(1);
        end
      end
      Nx = Nx + length(mdata.GPS_time);
    catch ME
      warning('Failed to load %s!!!!!!:\n  %s', data_fn, ME.getReport);
      mdata.Time = [0 1e-6];
      min_time = 0;
      max_time = 0;
      % keyboard % Uncomment for debugging why file loading failed
    end
  end
  dt = mdata.Time(2) - mdata.Time(1);
  Nt = 1 + (max_time-min_time)/dt;
  
  for track_idx = 1:param.layer_tracker.track_per_task:length(param.layer_tracker.track)
    dparam = [];
    dparam.file_success = {};
    dparam.argsin{1}.layer_tracker.frms = frms;
    
    tracks_in_task = track_idx:min(track_idx-1+param.layer_tracker.track_per_task,length(param.layer_tracker.track));
    
    dparam.argsin{1}.layer_tracker.tracks_in_task = tracks_in_task;
    
    % File Success
    % ---------------------------------------------------------------------
    for track_idx = tracks_in_task
      tmp_out_fn_name = sprintf('%s_%s.mat', param.layer_tracker.track{track_idx}.name, param.layer_tracker.track{track_idx}.method);
      tmp_out_fn = fullfile(tmp_out_fn_dir_dir,sprintf('layer_tracker_%03d', frm),tmp_out_fn_name);
      dparam.file_success{end+1} = tmp_out_fn;
      if ~ctrl.cluster.rerun_only && exist(tmp_out_fn,'file')
        delete(tmp_out_fn);
      end
    end
    
    % Rerun only check
    % ---------------------------------------------------------------------
    if ~ctrl.cluster.rerun_only
      if ~cluster_file_success(dparam.file_success)
        fprintf('  Already exists [rerun_only skipping]: %s (%s)\n', ...
          dparam.notes, datestr(now));
        continue;
      end
    end
    
    % CPU time and memory
    % ---------------------------------------------------------------------
    dparam.cpu_time = sum(cpu_time_mult(tracks_in_task)) * Nx * Nt;
    dparam.mem = 800e6 + max(mem_mult(tracks_in_task)) * Nx * Nt;
    if strcmp(track.init.method,'dem')
      % Add extra time and memory for DEM
      dparam.cpu_time = dparam.cpu_time + 120;
      dparam.mem = dparam.mem + 2e9;
    end
    mem_combine = mem_combine + 256*Nx*length(tracks_in_task);
    cputime_combine = cputime_combine + 1e-1*Nx*length(tracks_in_task);
    
    % Notes
    % ---------------------------------------------------------------------
    dparam.notes = sprintf('%s %s:%s:%s %s %s:%d-%d %s %d-%d (%d of %d)', ...
      sparam.task_function, param.radar_name, param.season_name, ...
      param.layer_tracker.echogram_source, param.layer_tracker.layer_params.layerdata_source, ...
      param.layer_tracker.track{tracks_in_task(1)}.method, tracks_in_task([1 end]), param.day_seg, ...
      dparam.argsin{1}.layer_tracker.frms([1 end]), start_frm_idx, length(param.cmd.frms));
    
    % Create task
    % ---------------------------------------------------------------------
    ctrl = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);
  end
end

ctrl = cluster_save_dparam(ctrl);
ctrl_chain{end+1} = ctrl;
fprintf('Done %s\n',datestr(now));

%% layer_tracker_combine
% =========================================================================
% =========================================================================
ctrl = cluster_new_batch(param);

sparam = [];
sparam.argsin{1} = param;
sparam.task_function = 'layer_tracker_combine_task';
sparam.num_args_out = 1;

sparam.cpu_time = 30 + cputime_combine;
sparam.mem = 500e6 + mem_combine;
sparam.notes = '';

if strcmp(param.layer_tracker.layer_params.source,'ops')
  sparam.file_success = {};
else
  sparam.file_success = {};
  out_fn_dir = ct_filename_out(param,'',param.layer_tracker.layer_params.layerdata_source);
  for frm = param.cmd.frms
    out_fn = fullfile(out_fn_dir,sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    sparam.file_success{end+1} = out_fn;
  end
end

ctrl = cluster_new_task(ctrl,sparam,[]);
ctrl_chain{end+1} = ctrl;
fprintf('Done %s\n',datestr(now));

