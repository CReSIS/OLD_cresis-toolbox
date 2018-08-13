function ctrl_chain = coh_noise_tracker(param,param_override)
% ctrl_chain = coh_noise_tracker(param,param_override)
%
% This function performs analyses on the data. A set of analyses commands
% can be passed in.
%
% param = struct with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Example:
%  See run_coh_noise_tracker.m for how to run this function directly.
%  Normally this function is called from master.m using the param spreadsheet.
%
% Authors: John Paden
%
% See also: master.m, run_coh_noise_tracker.m coh_noise_tracker.m,
%   coh_noise_tracker_task.m

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks
% =====================================================================

if ~isempty(param.cmd.frms)
  warning('All frames are processed with analysis, setting param.cmd.frms to do all frames.');
  param.cmd.frms = []; % All frames
end

if ~isfield(param.analysis,'specular') || isempty(param.analysis.specular.en)
  param.analysis.specular.en = 0;
end

if ~isfield(param.analysis.specular,'threshold_max') || isempty(param.analysis.specular.threshold_max)
  param.analysis.specular.threshold_max = 10;
end

if ~isfield(param.analysis,'coh_ave') || isempty(param.analysis.coh_ave.en)
  param.analysis.coh_ave.en = 0;
end

if ~isfield(param.analysis.coh_ave,'nz_list') || isempty(param.analysis.coh_ave.nz_list)
  param.analysis.coh_ave.nz_list = [1];
end

if ~isfield(param.analysis,'surf') || isempty(param.analysis.surf.en)
  param.analysis.surf.en = 0;
end

if ~isfield(param.analysis,'power') || isempty(param.analysis.power.en)
  param.analysis.power.en = 0;
end

if ~isfield(param.analysis,'psd') || isempty(param.analysis.psd.en)
  param.analysis.psd.en = 0;
end

if param.analysis.specular.en && param.analysis.coh_ave.en
  error('Run coh ave first, then specular');
end

if ~isfield(param.analysis,'out_dir') || isempty(param.analysis.out_dir)
  param.analysis.out_dir = 'analysis';
end

if ~isfield(param.get_heights,'trim_vals') || isempty(param.get_heights.trim_vals)
  param.get_heights.trim_vals = [0 0];
end

if ~isfield(param.get_heights,'coh_noise_method') || isempty(param.get_heights.coh_noise_method)
  param.get_heights.coh_noise_method = 0;
end

if ~isfield(param.get_heights,'coh_noise_arg')
  param.get_heights.coh_noise_arg = [];
end

if ~isfield(param.get_heights,'deconvolution') || isempty(param.get_heights.deconvolution)
  param.get_heights.deconvolution = 0;
end
if ~isfield(param.get_heights,'deconv_enforce_wf_idx') 
  param.get_heights.deconv_enforce_wf_idx = [];
end
if ~isfield(param.get_heights,'deconv_same_twtt_bin') 
  param.get_heights.deconv_same_twtt_bin = [];
end

if ~isfield(param.get_heights,'psd_smooth') || isempty(param.get_heights.psd_smooth)
  param.get_heights.psd_smooth = 0;
end

if ~isfield(param.get_heights,'ft_oversample') || isempty(param.get_heights.ft_oversample)
  param.get_heights.ft_oversample = 1;
end

if ~isfield(param.get_heights,'pulse_rfi') || isempty(param.get_heights.pulse_rfi)
  param.get_heights.pulse_rfi.en = 0;
end

if ~isfield(param.get_heights,'ft_dec') || isempty(param.get_heights.ft_dec)
  param.get_heights.ft_dec = 1;
end

if ~isfield(param.get_heights,'ft_wind_time') || isempty(param.get_heights.ft_wind_time)
  param.get_heights.ft_wind_time = 0;
end

if ~isfield(param.get_heights,'trim_vals') || isempty(param.get_heights.trim_vals)
  param.get_heights.trim_vals = 1;
end

if ~isfield(param.get_heights,'pulse_comp') || isempty(param.get_heights.pulse_comp)
  param.get_heights.pulse_comp = 1;
end

if ~isfield(param.get_heights,'raw_data') || isempty(param.get_heights.raw_data)
  param.get_heights.raw_data = 0;
end

if ~isfield(param.get_heights,'elev_correction') || isempty(param.get_heights.elev_correction)
  param.get_heights.elev_correction = false;
end

if ~isfield(param.get_heights,'roll_correction') || isempty(param.get_heights.roll_correction)
  param.get_heights.roll_correction = 0;
end

if abs(sum(param.get_heights.B_filter)-1) > 1e4*eps
  %warning('B_filter weights are not normalized. They must be normalized so normalizing to one now.')
  param.get_heights.B_filter = param.get_heights.B_filter / sum(param.get_heights.B_filter);
end

if ~isfield(param.get_heights,'inc_B_filter') || isempty(param.get_heights.inc_B_filter)
  param.get_heights.inc_B_filter = 1;
end
if abs(sum(param.get_heights.inc_B_filter)-1) > 1e4*eps
  %warning('inc_B_filter weights are not normalized. They must be normalized so normalizing to one now.')
  param.get_heights.inc_B_filter = param.get_heights.inc_B_filter / sum(param.get_heights.inc_B_filter);
end

if param.analysis.surf.en
  % Pulse compression should be enabled for surface waveform extraction
  param.get_heights.pulse_comp = 1;
end

if param.analysis.coh_ave.en
  % Coherent noise removal must be turned off to measure it
  param.get_heights.coh_noise_method    = 0;
  param.get_heights.coh_noise_arg       = [];
  % PSD smoothing must be turned off
  param.get_heights.psd_smooth = 0;
end

if param.analysis.specular.en
  % PSD smoothing must be turned off
  param.get_heights.psd_smooth = 0;
end

%% Setup processing
% =====================================================================

% Get the standard radar name
[~,~,radar_name] = ct_output_dir(param.radar_name);

% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);
  
%% Create and setup the cluster batch
% =====================================================================
ctrl = cluster_new_batch(param);
cluster_compile({'coh_noise_tracker_task.m','coh_noise_tracker_combine_task.m'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);

if any(strcmpi(radar_name,{'acords','hfrds','mcords','mcords2','mcords3','mcords4','mcords5','mcrds','seaice','accum2'}))
  [wfs,~] = load_mcords_wfs(records.settings, param, ...
    1:max(records.param_records.records.file.adcs), param.get_heights);
  for img = 1:length(param.analysis.imgs)
    wf = abs(param.analysis.imgs{img}(1,1));
    total_num_sam(img) = wfs(wf).Nt_raw;
  end
  cpu_time_mult = 140e-9;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  total_num_sam = 32000 * ones(size(param.analysis.imgs));
  cpu_time_mult = 8e-8;
  mem_mult = 33;
  
else
  error('radar_name %s not supported yet.', radar_name);
  
end

%% Split up data into blocks and run equal tasks
% =====================================================================
% Load param.analysis.block_size records at a time
%    --> The last block can be up to 1.5*param.analysis.block_size
out_recs = {};
retry_fields = {};

% Break records in segment into blocks
breaks = 1:param.analysis.block_size:length(records.gps_time);

% If the last block is less than half the desired block size, then combine
% with earlier block if possible
if length(records.gps_time)-breaks(end) < param.analysis.block_size/2 ...
    && length(breaks) > 1
  breaks = breaks(1:end-1);
end

% Create output directory string
out_fn_dir = ct_filename_out(param,param.analysis.out_dir);

sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'coh_noise_tracker_task';
sparam.num_args_out = 1;
sparam.argsin{1}.load.imgs = param.analysis.imgs;
for break_idx = 1:length(breaks)
  % Determine the start/stop record for this block
  rec_load_start = breaks(break_idx);
  if break_idx == length(breaks)
    rec_load_stop = length(records.gps_time);
  else
    rec_load_stop = rec_load_start+param.analysis.block_size-1;
  end
  cur_recs = [rec_load_start rec_load_stop];
  
  % Prepare task inputs
  % =================================================================
  dparam = [];
  dparam.argsin{1}.load.recs = cur_recs;
  
  % Create success condition
  % =================================================================
  dparam.success = '';
  success_error = 64;
  for img = 1:length(param.analysis.imgs)
    if param.analysis.coh_ave.en
      out_fn = fullfile(out_fn_dir,sprintf('coh_noise_img_%02d_%d_%d.mat',img,cur_recs(1),cur_recs(end)));
      dparam.success = cat(2,dparam.success, ...
        sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
      if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
        delete(out_fn);
      end
    end
    if param.analysis.specular.en
      out_fn = fullfile(out_fn_dir,sprintf('specular_img_%02d_%d_%d.mat',img,cur_recs(1),cur_recs(end)));
      dparam.success = cat(2,dparam.success, ...
        sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
      if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
        delete(out_fn);
      end
    end
    if param.analysis.surf.en
      out_fn = fullfile(out_fn_dir,sprintf('surf_img_%02d_%d_%d.mat',cur_recs(1),cur_recs(end)));
      dparam.success = cat(2,dparam.success, ...
        sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
      if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
        delete(out_fn);
      end
    end
  end
  
  % Rerun only mode: Test to see if we need to run this task
  % =================================================================
  dparam.notes = sprintf('%s:%s:%s %s %d of %d recs %d-%d', ...
    mfilename, param.radar_name, param.season_name, param.day_seg, ...
    break_idx, length(breaks), cur_recs(1), cur_recs(end));
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
  Nx = cur_recs(end)-cur_recs(1)+1;
  dparam.cpu_time = 0;
  dparam.mem = 0;
  for img = 1:length(param.analysis.imgs)
    if param.analysis.coh_ave.en
      dparam.cpu_time = dparam.cpu_time + 10 + size(param.analysis.imgs{img},1)*Nx*total_num_sam(img)*log2(Nx)*cpu_time_mult;
      dparam.mem = max(dparam.mem,250e6 + Nx*total_num_sam(img)*mem_mult);
    end
    if param.analysis.specular.en
      dparam.cpu_time = dparam.cpu_time + 10 + size(param.analysis.imgs{img},1)*Nx*total_num_sam(img)*log2(total_num_sam(img))*cpu_time_mult;
      dparam.mem = max(dparam.mem,250e6 + Nx*total_num_sam(img)*mem_mult);
    end
    if param.analysis.surf.en
      dparam.cpu_time = dparam.cpu_time + 10 + size(param.analysis.imgs{img},1)*Nx*total_num_sam(img)*log2(total_num_sam(img))*cpu_time_mult;
      dparam.mem = max(dparam.mem,250e6 + Nx*total_num_sam(img)*mem_mult);
    end
  end
  
  ctrl = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);
  
end

ctrl = cluster_save_dparam(ctrl);

ctrl_chain = {ctrl};

%% Create and setup the combine batch
% =====================================================================
ctrl = cluster_new_batch(param);

if any(strcmpi(radar_name,{'acords','hfrds','mcords','mcords2','mcords3','mcords4','mcords5','mcrds','seaice','accum2'}))
  cpu_time_mult = 6e-6;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  cpu_time_mult = 1000e-8;
  mem_mult = 24;
end

sparam = [];
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'coh_noise_tracker_combine_task';
sparam.num_args_out = 1;
sparam.cpu_time = 60;
sparam.mem = 0;
% Add up all records being processed and find the most records in a block
Nx = length(records.gps_time);
for img = 1:length(param.analysis.imgs)
  Nt = total_num_sam(img);
  if param.analysis.coh_ave.en
    Nx_cmd = Nx / param.analysis.coh_ave.block_ave;
    sparam.cpu_time = sparam.cpu_time + size(param.analysis.imgs{img},1)*Nx_cmd*Nt*cpu_time_mult;
    sparam.mem = max(sparam.mem,250e6 + size(param.analysis.imgs{img},1)*Nx_cmd*Nt*mem_mult);
  end
  if param.analysis.specular.en
    Nx_cmd = Nx / param.analysis.block_size * param.analysis.specular.threshold_max;
    sparam.cpu_time = sparam.cpu_time + size(param.analysis.imgs{img},1)*Nx_cmd*Nt*cpu_time_mult;
    sparam.mem = max(sparam.mem,250e6 + size(param.analysis.imgs{img},1)*Nx_cmd*Nt*mem_mult);
  end
  if param.analysis.surf.en
    Nx_cmd = Nx / param.get_heights.decimate_factor;
    if isfinite(param.analysis.surf.Nt)
      Nt = param.analysis.surf.Nt;
    end
    sparam.cpu_time = sparam.cpu_time + size(param.analysis.imgs{img},1)*Nx_cmd*Nt*cpu_time_mult;
    sparam.mem = max(sparam.mem,250e6 + size(param.analysis.imgs{img},1)*Nx_cmd*Nt*mem_mult);
  end
end
sparam.notes = sprintf('%s:%s:%s %s combine', ...
  mfilename, param.radar_name, param.season_name, param.day_seg);

% Create success condition
success_error = 64;
sparam.success = '';

out_fn_dir_dir = fileparts(out_fn_dir);
for img = 1:length(param.analysis.imgs)
  if param.analysis.coh_ave.en
    out_fn = fullfile(out_fn_dir_dir,sprintf('coh_noise_img_%02d.mat',img));
    dparam.success = cat(2,dparam.success, ...
      sprintf('  error_mask = bitor(error_mask,%d*~exist(%s,''file''));\n', success_error, out_fn));
  end
  if param.analysis.specular.en
    out_fn = fullfile(out_fn_dir_dir,sprintf('specular_img_%02d.mat',img));
    dparam.success = cat(2,dparam.success, ...
      sprintf('  error_mask = bitor(error_mask,%d*~exist(%s,''file''));\n', success_error, out_fn));
  end
  if param.analysis.surf.en
    out_fn = fullfile(out_fn_dir_dir,sprintf('surf_img_%02d.mat',img));
    dparam.success = cat(2,dparam.success, ...
      sprintf('  error_mask = bitor(error_mask,%d*~exist(%s,''file''));\n', success_error, out_fn));
  end
end

ctrl = cluster_new_task(ctrl,sparam,[]);

ctrl_chain{end+1} = ctrl;
    
fprintf('Done %s\n', datestr(now));
