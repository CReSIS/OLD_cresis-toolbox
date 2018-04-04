function ctrl_chain = csarp(param,param_override)
% ctrl_chain = csarp(param,param_override)
%
% SAR processor function which breaks up the frame into chunks
% which are processed by csarp_task.m
%
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Authors: William Blake, John Paden
%
% Example:
%  See run_csarp.m for how to run this function directly.
%  Normally this function is called from master.m using the param spreadsheet.
%
% See also: run_master.m, master.m, run_csarp.m, csarp.m,
%   csarp_task.m

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks
% =====================================================================

% Get speed of light, dielectric of ice constants
physical_constants;

if ~isfield(param.csarp,'frm_types') || isempty(param.csarp.frm_types)
  param.csarp.frm_types = {-1,-1,-1,-1,-1};
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

if ~isfield(param.csarp,'out_path') || isempty(param.csarp.out_path)
  param.csarp.out_path = 'out';
end

if ~isfield(param.csarp,'coord_path') || isempty(param.csarp.coord_path)
  param.csarp.coord_path = param.csarp.out_path;
end

if ~isfield(param.csarp,'pulse_comp') || isempty(param.csarp.pulse_comp)
  param.csarp.pulse_comp = 1;
end

if ~isfield(param.csarp,'ground_based') || isempty(param.csarp.ground_based)
  param.csarp.ground_based = false;
end

if ~isfield(param.csarp,'sar_type') || isempty(param.csarp.sar_type)
  param.csarp.sar_type = 'fk';
end

if strcmpi(param.csarp.sar_type,'f-k')
  error('Deprecated sar_type name. Change param.csarp.sar_type from ''f-k'' to ''fk'' in  your parameters (or remove parameter since ''fk'' is the default mode).');
end

if ~isfield(param.csarp,'pulse_comp') || isempty(param.csarp.pulse_comp)
  param.csarp.pulse_comp = true;
end

if ~isfield(param.csarp,'ft_dec') || isempty(param.csarp.ft_dec)
  param.csarp.ft_dec = true;
end

if ~isfield(param.csarp,'ft_wind_time') || isempty(param.csarp.ft_wind_time)
  param.csarp.ft_wind_time = [];
end

if ~isfield(param.csarp,'start_eps') || isempty(param.csarp.start_eps)
  param.csarp.start_eps = er_ice;
end

if ~isfield(param.csarp,'time_of_full_support') || isempty(param.csarp.time_of_full_support)
  param.csarp.time_of_full_support = inf;
end

if ~isfield(param.csarp,'force_one_wf_adc_pair_per_job') || isempty(param.csarp.force_one_wf_adc_pair_per_job)
  param.csarp.force_one_wf_adc_pair_per_job = false;
end

if ~isfield(param.csarp,'combine_rx') || isempty(param.csarp.combine_rx)
  param.csarp.combine_rx = false;
end

if ~isfield(param.csarp,'Lsar')
  param.csarp.Lsar = [];
end

if ~isfield(param.csarp.Lsar,'agl')
  param.csarp.Lsar.agl = 500;
end

if ~isfield(param.csarp.Lsar,'thick')
  param.csarp.Lsar.thick = 1000;
end

if ~isfield(param.csarp,'trim_vals') || isempty(param.csarp.trim_vals)
  param.csarp.trim_vals = [0 0];
end

if ~isfield(param.csarp,'coh_noise_removal') || isempty(param.csarp.coh_noise_removal) ...
    || ~param.csarp.coh_noise_removal
  param.csarp.coh_noise_removal = 0;
  default_coh_noise_method = 0;
elseif param.csarp.coh_noise_removal
  default_coh_noise_method = 1;
end

if ~isfield(param.csarp,'coh_noise_method') || isempty(param.csarp.coh_noise_method)
  param.csarp.coh_noise_method = default_coh_noise_method;
end

if ~isfield(param.csarp,'coh_noise_arg')
  param.csarp.coh_noise_arg = [];
end

if ~isfield(param.csarp,'pulse_rfi') || isempty(param.csarp.pulse_rfi)
  param.csarp.pulse_rfi.en = 0;
end

if ~isfield(param.csarp.mocomp,'filter')
  param.csarp.mocomp.filter = '';
end

if ~isfield(param.csarp,'surf_filt_dist') || isempty(param.csarp.surf_filt_dist)
  param.csarp.surf_filt_dist = 3e3;
  warning('Surface filtering is not set. Using default value: %.0f m.',param.csarp.surf_filt_dist);
end

if ~isfield(param.csarp,'presums')
  param.csarp.presums = 1;
end

if ~isfield(param.csarp,'deconvolution') || isempty(param.csarp.deconvolution)
  param.csarp.deconvolution = 0;
end

if ~isfield(param.csarp,'psd_smooth') || isempty(param.csarp.psd_smooth)
  param.csarp.psd_smooth = 0;
end

% Do not apply channel equalization during csarp combine unless receivers
% are being combined at this stage (csarp-combined method)
if ~param.csarp.combine_rx
  for wf = 1:length(param.radar.wfs)
    param.radar.wfs(wf).chan_equal_dB(:) = 0;
    param.radar.wfs(wf).chan_equal_deg(:) = 0;
  end
end

%% Setup processing
% =====================================================================

% Get the standard radar name
[~,~,radar_name] = ct_output_dir(param.radar_name);

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

% SAR output directory
csarp_out_dir = ct_filename_out(param, param.csarp.out_path);
csarp_coord_dir = ct_filename_out(param, param.csarp.coord_path);

ctrl_chain = {};

%% Collect waveform information into one structure
%  - This is used to break the frame up into chunks
% =====================================================================
if strcmpi(radar_name,'mcrds')
  wfs = load_mcrds_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
elseif any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
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

%% Create the SAR coordinate system (used for structuring processing)
% =====================================================================

sar_fn = fullfile(csarp_coord_dir,'sar_coord.mat');
if exist(sar_fn,'file')
  sar = load(sar_fn,'Lsar','gps_source','type','sigma_x','presums','version');
end
Lsar = c/wfs(1).fc*(param.csarp.Lsar.agl+param.csarp.Lsar.thick/sqrt(er_ice))/(2*param.csarp.sigma_x);
if ~exist(sar_fn,'file') ...
    || sar.Lsar ~= Lsar ...
    || ~strcmpi(sar.gps_source,records.gps_source) ...
    || sar.type ~= param.csarp.mocomp.type ...
    || sar.sigma_x ~= param.csarp.sigma_x ...
    || sar.presums ~= param.csarp.presums ...
    || sar.version ~= 1.0
  
  ctrl = cluster_new_batch(param);
  
  if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
    cpu_time_mult = 6e-3;
    mem_mult = 64;
    
  elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
    cpu_time_mult = 100e-8;
    mem_mult = 64;
  end
  
  sparam = [];
  sparam.argsin{1} = param; % Static parameters
  sparam.task_function = 'csarp_sar_coord_task';
  sparam.num_args_out = 1;
  Nx = numel(records.gps_time);
  sparam.cpu_time = 60 + Nx*cpu_time_mult;
  sparam.mem = 250e6 + Nx*mem_mult;
  sparam.notes = sprintf('%s:%s:%s %s sar coordinate system', ...
    mfilename, param.radar_name, param.season_name, param.day_seg);
    
  % Create success condition
  success_error = 64;
  sparam.success = ...
    sprintf('error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, sar_fn);
  
  ctrl = cluster_new_task(ctrl,sparam,[]);
  
  ctrl_chain{end+1} = ctrl;
end

%% Create imgs_list which divides up images into tasks
% =====================================================================
% imgs_list: cell array of images
%  imgs_list{img}: cell array of wf-adc pair lists PER task
%   imgs_list{img}{task}: array of wf-adc pairs to be processed by a task
imgs_list = {};
if param.csarp.combine_rx
  % One SAR image with all wf-adc pairs
  for img = 1:length(param.csarp.imgs)
    imgs_list{1}{img} = param.csarp.imgs{img};
  end
  
else
  % One SAR image per wf-adc pair

  if ~param.csarp.force_one_wf_adc_pair_per_job
    % All SAR images from the same data files per task
    for img = 1:length(param.csarp.imgs)
      for wf_adc = 1:size(param.csarp.imgs{img},1)
        adc = abs(param.csarp.imgs{img}(wf_adc,2));
        [board,board_idx] = adc_to_board(radar_name,adc);
        if length(imgs_list) < board_idx
          imgs_list{board_idx} = [];
        end
        if length(imgs_list{board_idx}) < img
          imgs_list{board_idx}{img} = [];
        end
        imgs_list{board_idx}{img}(end+1,:) = param.csarp.imgs{img}(wf_adc,:);
      end
    end
    
  else
    % One SAR image per task
    for img = 1:length(param.csarp.imgs)
      for wf_adc = 1:length(param.csarp.imgs{img})
        imgs_list{end+1}{1}(1,:) = param.csarp.imgs{img}(wf_adc,:);
      end
    end
  end
end

%% Create and setup the cluster batch
% =====================================================================
ctrl = cluster_new_batch(param);
cluster_compile({'csarp_task.m','csarp_sar_coord_task'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);

total_num_sam = {};
if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  for imgs_idx = 1:length(imgs_list)
    for img = 1:length(imgs_list{imgs_idx})
      wf = abs(imgs_list{imgs_idx}{img}(1,1));
      total_num_sam{imgs_idx}{img} = wfs(wf).Nt_raw;
    end
  end
  cpu_time_mult = 150e-8;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  for imgs_idx = 1:length(imgs_list)
    for img = 1:length(imgs_list{imgs_idx})
      wf = abs(imgs_list{imgs_idx}{img}(1,1));
      total_num_sam{imgs_idx}{img} = 32000;
    end
  end
  cpu_time_mult = 8e-8;
  mem_mult = 64;
  
else
  error('radar_name %s not supported yet.', radar_name);
  
end

%% Setup tasks to SAR process each frame
% =====================================================================
retry_fields = {};
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'csarp_task';
sparam.num_args_out = 1;
for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  if ct_proc_frame(frames.proc_mode(frm),param.csarp.frm_types)
    fprintf('%s %s_%03i (%i of %i) %s\n', mfilename, param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now,'HH:MM:SS'));
  else
    fprintf('Skipping frame %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end

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
  num_chunks = round(frm_dist / param.csarp.chunk_len);
  
  % Estimate number of input range lines per chunk
  num_rlines_per_chunk = round((stop_rec-start_rec) / num_chunks);
  
  for chunk_idx = 1:num_chunks
    % Setup dynamic params
    % =====================================================================
    dparam.argsin{1}.load.frm = frm;
    dparam.argsin{1}.load.chunk_idx = chunk_idx;
    dparam.argsin{1}.load.num_chunks = num_chunks;
    if chunk_idx == num_chunks
      dparam.argsin{1}.load.recs = [start_rec + num_rlines_per_chunk*(chunk_idx-1), stop_rec];
    else
      dparam.argsin{1}.load.recs = start_rec + num_rlines_per_chunk*(chunk_idx-1) + [0, num_rlines_per_chunk-1];
    end
    
    for imgs_idx = 1:length(imgs_list)
      dparam.argsin{1}.load.imgs = imgs_list{imgs_idx};
      sub_band_idx = 1;
      dparam.argsin{1}.load.sub_band_idx = sub_band_idx;
    
      % Create a list of waveform/adc pairs for displaying to the string.
      wf_adc_str = '';
      for img = 1:length(dparam.argsin{1}.load.imgs)
        for wf_adc_idx = 1:size(dparam.argsin{1}.load.imgs{img},1)
          wf = abs(dparam.argsin{1}.load.imgs{img}(wf_adc_idx,1));
          adc = abs(dparam.argsin{1}.load.imgs{img}(wf_adc_idx,2));
          if isempty(wf_adc_str)
            wf_adc_str = [wf_adc_str sprintf('wf,adc: %d,%d',wf,adc)];
          else
            wf_adc_str = [wf_adc_str sprintf(' %d,%d',wf,adc)];
          end
        end
      end
      
      % Create success condition
      % =================================================================
      dparam.success = '';
      success_error = 64;
      for subap = 1:length(param.csarp.sub_aperture_steering)
        out_fn_dir = fullfile(csarp_out_dir, ...
          sprintf('%s_data_%03d_%02d_%02d',param.csarp.sar_type,frm, ...
          subap, sub_band_idx));
        for img = 1:length(dparam.argsin{1}.load.imgs)
          if param.csarp.combine_rx
            out_fn = fullfile(out_fn_dir,sprintf('img_%02d_chk_%03d.mat',img,chunk_idx));
            dparam.success = cat(2,dparam.success, ...
              sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
            if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
              delete(out_fn);
            end
          else
            for wf_adc = 1:size(dparam.argsin{1}.load.imgs{img},1)
              wf  = abs(dparam.argsin{1}.load.imgs{img}(wf_adc,1));
              adc = abs(dparam.argsin{1}.load.imgs{img}(wf_adc,2));
              out_fn = fullfile(out_fn_dir,sprintf('wf_%02d_adc_%02d_chk_%03d.mat',wf,adc,chunk_idx));
              dparam.success = cat(2,dparam.success, ...
                sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
              if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
                delete(out_fn);
              end
            end
          end
        end
      end
      
      % Rerun only mode: Test to see if we need to run this task
      % =================================================================
      dparam.notes = sprintf('%s:%s:%s %s_%03d (%d of %d)/%d of %d %s %.0f to %.0f recs', ...
        mfilename, param.radar_name, param.season_name, param.day_seg, frm, frm_idx, length(param.cmd.frms), ...
        chunk_idx, num_chunks, wf_adc_str, (dparam.argsin{1}.load.recs(1)-1)*param.csarp.presums+1, ...
        dparam.argsin{1}.load.recs(2)*param.csarp.presums);
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
      Nx = diff(dparam.argsin{1}.load.recs);
      dparam.cpu_time = 0;
      dparam.mem = 250e6;
      for img = 1:length(dparam.argsin{1}.load.imgs)
        if strcmpi(param.csarp.sar_type,'fk')
          dparam.cpu_time = dparam.cpu_time + 10 + Nx*log2(Nx)*total_num_sam{imgs_idx}{img}*log2(total_num_sam{imgs_idx}{img})*size(dparam.argsin{1}.load.imgs{img},1)*cpu_time_mult;
          dparam.mem = dparam.mem + Nx*total_num_sam{imgs_idx}{img}*size(dparam.argsin{1}.load.imgs{img},1)*mem_mult;
        elseif strcmpi(param.csarp.sar_type,'tdbp')
        end
      end
      
      ctrl = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);
    end
  end
end

ctrl = cluster_save_dparam(ctrl);

ctrl_chain{end+1} = ctrl;

fprintf('Done %s\n', datestr(now));

return;

