function ctrl_chain = sar(param,param_override)
% ctrl_chain = sar(param,param_override)
%
% SAR processor function which breaks up the frame into chunks
% which are processed by sar_task.m
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
%  See run_sar.m for how to run this function directly.
%  Normally this function is called from master.m using the param spreadsheet.
%
% See also: run_master.m, master.m, run_sar.m, sar.m, sar_task.m,
%   sar_coord_task.m

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

if ~isfield(param.sar,'combine_rx') || isempty(param.sar.combine_rx)
  param.sar.combine_rx = false;
end

if ~isfield(param.sar,'force_one_wf_adc_pair_per_job') || isempty(param.sar.force_one_wf_adc_pair_per_job)
  param.sar.force_one_wf_adc_pair_per_job = false;
end

if ~isfield(param.sar,'frm_types') || isempty(param.sar.frm_types)
  param.sar.frm_types = {-1,-1,-1,-1,-1};
end

if ~isfield(param.sar,'Lsar')
  param.sar.Lsar = [];
end
if ~isfield(param.sar.Lsar,'agl')
  param.sar.Lsar.agl = 500;
end
if ~isfield(param.sar.Lsar,'thick')
  param.sar.Lsar.thick = 1000;
end

if ~isfield(param.sar.mocomp,'filter')
  param.sar.mocomp.filter = '';
end

if ~isfield(param.sar,'out_path') || isempty(param.sar.out_path)
  param.sar.out_path = 'out';
end
if ~isfield(param.sar,'coord_path') || isempty(param.sar.coord_path)
  param.sar.coord_path = param.sar.out_path;
end

if ~isfield(param.sar,'presums')
  param.sar.presums = 1;
end

if ~isfield(param.sar,'sar_type') || isempty(param.sar.sar_type)
  param.sar.sar_type = 'fk';
end
if strcmpi(param.sar.sar_type,'f-k')
  error('Deprecated sar_type name. Change param.sar.sar_type from ''f-k'' to ''fk'' in  your parameters (or remove parameter since ''fk'' is the default mode).');
end

if ~isfield(param.sar,'start_eps') || isempty(param.sar.start_eps)
  param.sar.start_eps = er_ice;
end

if ~isfield(param.sar,'surf_filt_dist') || isempty(param.sar.surf_filt_dist)
  param.sar.surf_filt_dist = 3e3;
  warning('Surface filtering is not set. Using default value: %.0f m.',param.sar.surf_filt_dist);
end

if ~isfield(param.sar,'surf_layer') || isempty(param.sar.surf_layer)
  param.sar.surf_layer.name = 'surface';
  param.sar.surf_layer.source = 'layerData';
end
% Never check for the existence of layers
param.array.surf_layer.existence_check = false;

if ~isfield(param.sar,'time_of_full_support') || isempty(param.sar.time_of_full_support)
  param.sar.time_of_full_support = inf;
end

%% Setup processing
% =====================================================================

% Do not apply channel equalization during sar processing unless receivers
% are being combined at this stage (combine_rx is true)
if ~param.sar.combine_rx
  for wf = 1:length(param.radar.wfs)
    param.radar.wfs(wf).chan_equal_dB(:) = 0;
    param.radar.wfs(wf).chan_equal_deg(:) = 0;
  end
end

% Get the standard radar name
[~,~,radar_name] = ct_output_dir(param.radar_name);

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

% SAR output directory
sar_out_dir = ct_filename_out(param, param.sar.out_path);
sar_coord_dir = ct_filename_out(param, param.sar.coord_path);

ctrl_chain = {};

%% Collect waveform information into one structure
%  - This is used to break the frame up into chunks
% =====================================================================
[wfs,~] = data_load_wfs(setfield(param,'load',struct('imgs',{param.sar.imgs})),records);
param.radar.wfs = merge_structs(param.radar.wfs,wfs);

%% Create the SAR coordinate system (used for structuring processing)
% =====================================================================

sar_fn = fullfile(sar_coord_dir,'sar_coord.mat');
fprintf('sar_coord_task %s (%s):\n  %s\n', datestr(now), param.day_seg, sar_fn);
if exist(sar_fn,'file')
  sar_fields = {'Lsar','gps_source','gps_time_offset','type','sigma_x','presums','version','along_track'};
  fcs = load(sar_fn,sar_fields{:});
  for sar_field = sar_fields
    if ~isfield(fcs,sar_field{1})
      error('SAR coordinates file missing %s field. Either delete file or correct file and then rerun.', sar_field{1});
    end
  end
end
Lsar = c/wfs(1).fc*(param.sar.Lsar.agl+param.sar.Lsar.thick/sqrt(er_ice))/(2*param.sar.sigma_x);
if ~exist(sar_fn,'file') ...
    || fcs.Lsar ~= Lsar ...
    || ~strcmpi(fcs.gps_source,records.gps_source) ...
    || any(fcs.gps_time_offset ~= records.param_records.records.gps.time_offset) ...
    || fcs.type ~= param.sar.mocomp.type ...
    || fcs.sigma_x ~= param.sar.sigma_x ...
    || fcs.presums ~= param.sar.presums ...
    || fcs.version ~= 1.0
  %% SAR coordinates file does not exist or needs to be recreated
  
  ctrl = cluster_new_batch(param);
  
  if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','mcrds','seaice','accum2'}))
    cpu_time_mult = 6e-3;
    mem_mult = 64;
    
  elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
    cpu_time_mult = 100e-8;
    mem_mult = 64;
  end
  
  sparam = [];
  sparam.argsin{1} = param; % Static parameters
  sparam.task_function = 'sar_coord_task';
  sparam.num_args_out = 1;
  Nx = numel(records.gps_time);
  sparam.cpu_time = 60 + Nx*cpu_time_mult;
  sparam.mem = 250e6 + Nx*mem_mult;
  sparam.notes = sprintf('%s:%s:%s %s', ...
    sparam.task_function, param.radar_name, param.season_name, param.day_seg);
    
  % Create success condition
  success_error = 64;
  sparam.success = ...
    sprintf('error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, sar_fn);
  if ~ctrl.cluster.rerun_only && exist(sar_fn,'file')
    delete(sar_fn);
  end
  
  ctrl = cluster_new_task(ctrl,sparam,[]);
  
  ctrl_chain{end+1} = ctrl;

else
  %% SAR coordinates file exists and only needs surface updated
  surf_idxs = get_equal_alongtrack_spacing_idxs(fcs.along_track,param.sar.sigma_x);
  if surf_idxs(end) ~= length(fcs.along_track)
    surf_idxs(end+1) = length(fcs.along_track);
  end
  
  % Load surface from layerdata
  tmp_param = param;
  tmp_param.cmd.frms = [];
  surf_layer = opsLoadLayers(tmp_param,param.sar.surf_layer);
  if isempty(surf_layer.gps_time)
    records.surface(:) = 0;
  elseif length(surf_layer.gps_time) == 1
    records.surface(:) = surf_layer.twtt;
  else
    records.surface = interp_finite(interp1(surf_layer.gps_time,surf_layer.twtt,records.gps_time),0);
  end

  surf = sgolayfilt(records.surface(surf_idxs),3,round(param.sar.surf_filt_dist / median(diff(fcs.along_track(surf_idxs)))/2)*2+1);
  fcs.surf_pp = spline(fcs.along_track(surf_idxs),surf);
  
  fprintf('  Exists. Updating with new surface values from param.sar.surf_layer.\n');
  save(sar_fn,'-append','-struct','fcs','surf_pp');
end

%% Create imgs_list
% =====================================================================
% Create imgs_list which divides up imgs into several imgs in a data IO
% efficient way so that all SAR images from a particular data stream are
% done together. The assumption is that each board_idx represents a single
% data stream.
%
% imgs_list: cell array of new broken apart imgs
%  imgs_list{imgs_idx}: cell array of wf-adc pair lists PER task
%   imgs_list{imgs_idx}{task}: array of wf-adc pairs to be processed by a task
imgs_list = {};
if param.sar.combine_rx
  % One SAR image with all wf-adc pairs
  for img = 1:length(param.sar.imgs)
    imgs_list{1}{img} = param.sar.imgs{img};
  end
  
else
  % One SAR image per wf-adc pair

  if ~param.sar.force_one_wf_adc_pair_per_job
    % All SAR images from the same data files per task
    for img = 1:length(param.sar.imgs)
      for wf_adc = 1:size(param.sar.imgs{img},1)
        wf = param.sar.imgs{img}(wf_adc,1);
        adc = param.sar.imgs{img}(wf_adc,2);
        [board,board_idx,profile] = wf_adc_to_board(param,[wf adc]);

        if length(imgs_list) < board_idx
          imgs_list{board_idx} = [];
        end
        if length(imgs_list{board_idx}) < img
          imgs_list{board_idx}{img} = [];
        end
        imgs_list{board_idx}{img}(end+1,:) = param.sar.imgs{img}(wf_adc,:);
      end
    end
    
  else
    % One SAR image per task
    for img = 1:length(param.sar.imgs)
      for wf_adc = 1:length(param.sar.imgs{img})
        imgs_list{end+1}{1}(1,:) = param.sar.imgs{img}(wf_adc,:);
      end
    end
  end
end

%% Create and setup the cluster batch
% =====================================================================
ctrl = cluster_new_batch(param);
cluster_compile({'sar_task.m','sar_coord_task'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);

total_num_sam = {};
if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','mcrds','seaice','accum2'}))
  for imgs_idx = 1:length(imgs_list)
    for img = 1:length(imgs_list{imgs_idx})
      wf = abs(imgs_list{imgs_idx}{img}(1,1));
      total_num_sam{imgs_idx}{img} = wfs(wf).Nt_raw;
    end
  end
  cpu_time_mult = 150e-8/100;
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
sparam.task_function = 'sar_task';
sparam.num_args_out = 1;
for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  if ct_proc_frame(frames.proc_mode(frm),param.sar.frm_types)
    fprintf('%s %s_%03i (%i of %i) (%s)\n', sparam.task_function, param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now));
  else
    fprintf('Skipping frame %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end

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
  num_chunks = round(frm_dist / param.sar.chunk_len);
  
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
      if isempty(imgs_list{imgs_idx})
        continue;
      end
      dparam.argsin{1}.load.imgs = imgs_list{imgs_idx};
      sub_band_idx = 1;
      dparam.argsin{1}.load.sub_band_idx = sub_band_idx;
    
      % Create a list of waveform/adc pairs for displaying to the string.
      wf_adc_str = '';
      for img = 1:length(dparam.argsin{1}.load.imgs)
        for wf_adc = 1:size(dparam.argsin{1}.load.imgs{img},1)
          wf = abs(dparam.argsin{1}.load.imgs{img}(wf_adc,1));
          adc = abs(dparam.argsin{1}.load.imgs{img}(wf_adc,2));
          if isempty(wf_adc_str)
            wf_adc_str = [wf_adc_str sprintf('wf,adc: %d,%d',wf,adc)];
          else
            wf_adc_str = [wf_adc_str sprintf(' %d,%d',wf,adc)];
          end
        end
      end
      
      % Create success condition
      % =================================================================
      dparam.file_success = {};
      for subap = 1:length(param.sar.sub_aperture_steering)
        out_fn_dir = fullfile(sar_out_dir, ...
          sprintf('%s_data_%03d_%02d_%02d',param.sar.sar_type,frm, ...
          subap, sub_band_idx));
        for img = 1:length(dparam.argsin{1}.load.imgs)
          if param.sar.combine_rx
            out_fn = fullfile(out_fn_dir,sprintf('img_%02d_chk_%03d.mat',img,chunk_idx));
            
            dparam.file_success{end+1} = out_fn;
            if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
              delete(out_fn);
            end
          else
            for wf_adc = 1:size(dparam.argsin{1}.load.imgs{img},1)
              wf  = abs(dparam.argsin{1}.load.imgs{img}(wf_adc,1));
              adc = abs(dparam.argsin{1}.load.imgs{img}(wf_adc,2));
              out_fn = fullfile(out_fn_dir,sprintf('wf_%02d_adc_%02d_chk_%03d.mat',wf,adc,chunk_idx));
              
              dparam.file_success{end+1} = out_fn;
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
        sparam.task_function, param.radar_name, param.season_name, param.day_seg, frm, frm_idx, length(param.cmd.frms), ...
        chunk_idx, num_chunks, wf_adc_str, (dparam.argsin{1}.load.recs(1)-1)*param.sar.presums+1, ...
        dparam.argsin{1}.load.recs(2)*param.sar.presums);
      if ctrl.cluster.rerun_only
        % If we are in rerun only mode AND the sar task file success
        % condition passes without error, then we do not run the task.
        if ~cluster_file_success(dparam.file_success)
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
        if strcmpi(param.sar.sar_type,'fk')
          dparam.cpu_time = dparam.cpu_time + 10 + Nx*log2(Nx)*total_num_sam{imgs_idx}{img}*log2(total_num_sam{imgs_idx}{img})*size(dparam.argsin{1}.load.imgs{img},1)*cpu_time_mult;
          dparam.mem = dparam.mem + Nx*total_num_sam{imgs_idx}{img}*size(dparam.argsin{1}.load.imgs{img},1)*mem_mult;
        elseif strcmpi(param.sar.sar_type,'tdbp')
          dparam.cpu_time = ctrl.cluster.max_time_per_job - 40;
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

