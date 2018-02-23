function ctrl_chain = get_heights(param,param_override)
% ctrl_chain = get_heights(param,param_override)
%
% This function generates quick look outputs (CSARP_qlook), tracks the
% surface, and (optionally) stores the surface to records.
%
% param = struct with processing parameters
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

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks
% =====================================================================

if ~isfield(param.get_heights,'frm_types') || isempty(param.get_heights.frm_types)
  param.get_heights.frm_types = {-1,-1,-1,-1,-1};
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

if ~isfield(param.get_heights,'combine_only') || isempty(param.get_heights.combine_only)
  param.get_heights.combine_only = false;
end

if ~isfield(param.get_heights,'img_comb') || isempty(param.get_heights.img_comb)
  param.get_heights.img_comb = [];
end

% Convert inputs to new format that does not use qlook substruct
if isfield(param.get_heights,'qlook')
  warning('The get_heights.qlook field is deprecated. Please remove the qlook portion of each field (i.e. qlook.out_path should be just out_path).');
  qlook_fieldnames = fieldnames(param.get_heights.qlook);
  for name_idx = 1:length(qlook_fieldnames)
    param.get_heights.(qlook_fieldnames{name_idx}) = param.get_heights.qlook.(qlook_fieldnames{name_idx});
  end
  param.get_heights = rmfield(param.get_heights,'qlook');
end

if ~isfield(param.get_heights,'pulse_comp') || isempty(param.get_heights.pulse_comp)
  param.get_heights.pulse_comp = 1;
end

if ~isfield(param.get_heights,'out_path') || isempty(param.get_heights.out_path)
  param.get_heights.out_path = 'qlook';
end

if ~isfield(param.get_heights,'ground_based') || isempty(param.get_heights.ground_based)
  param.get_heights.ground_based = [];
end
if ~isfield(param.get_heights,'imgs') || isempty(param.get_heights.imgs)
  error('No images specified in param.get_heights.imgs. Nothing to do.');
end
if ~isfield(param.get_heights,'block_size') || isempty(param.get_heights.block_size)
  % [Block_Size Overlap]
  param.get_heights.block_size = [10000 0];
end
if numel(param.get_heights.block_size) == 1
  % Overlap
  param.get_heights.block_size(2) = 0;
end

% Check img_comb
if numel(param.get_heights.imgs) == 1 || isempty(param.get_heights.img_comb)
  num_imgs = 1;
else
  num_imgs = length(param.get_heights.imgs);
  if length(param.get_heights.img_comb) ~= 3*(num_imgs-1)
    error('param.get_heights.img_comb not the right length. Since it is not empty, there should be 3 entries for each image combination interface ([Tpd second image for surface saturation, -inf for second image blank, Tpd first image to avoid roll off] is typical). Set correctly here and update param spreadsheet before dbcont.');
  end
end

%% Setup Processing
% =====================================================================

% Get the standard radar name
[~,~,radar_name] = ct_output_dir(param.radar_name);

% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);

% Quick look radar echogram output directory
qlook_out_dir = ct_filename_out(param, param.get_heights.out_path);

% Get version information out of the deconvolution file
if isfield(param.get_heights,'deconvolution') ...
    && ~isempty(param.get_heights.deconvolution) ...
    && param.get_heights.deconvolution == 3
  out_fn_dir = ct_filename_out(param,'', 'CSARP_noise');
  out_segment_fn_dir = fileparts(out_fn_dir);
  out_segment_fn = fullfile(out_segment_fn_dir,sprintf('deconv_%s.mat', param.day_seg));
  spec = load(out_segment_fn,'param_collate');
  
  param.get_heights.deconvolution_sw_version = spec.param_collate.sw_version;
  param.get_heights.deconvolution_params = spec.param_collate.analysis.specular;
end

% Get version information out of the coherent noise file
if isfield(param.get_heights,'coh_noise_method') ...
    && ~isempty(param.get_heights.coh_noise_method) ...
    && any(param.get_heights.coh_noise_method == [17 19])
  
  cdf_fn_dir = fileparts(ct_filename_out(param,param.get_heights.coh_noise_arg{4}, ''));
  cdf_fn = fullfile(cdf_fn_dir,sprintf('coh_noise_simp_%s.nc', param.day_seg));
  
  tmp = netcdf_to_mat(cdf_fn,[],'^sw_version.*');
  param.get_heights.coh_noise_version = tmp.sw_version;
  tmp = netcdf_to_mat(cdf_fn,[],'^param_collate.*');
  param.get_heights.coh_noise_params = tmp.param_collate;
end

%% Create and setup the cluster batch
% =====================================================================
ctrl = cluster_new_batch(param);
cluster_compile({'get_heights_task.m','get_heights_combine_task.m'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);

if any(strcmpi(radar_name,{'acords','hfrds','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  [wfs,~] = load_mcords_wfs(records.settings, param, ...
    1:max(records.param_records.records.file.adcs), param.get_heights);
  for img = 1:length(param.get_heights.imgs)
    wf = abs(param.get_heights.imgs{img}(1,1));
    total_num_sam(img) = wfs(wf).Nt_raw;
  end
  cpu_time_mult = 66e-8;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  total_num_sam = 32000 * ones(size(param.get_heights.imgs));
  cpu_time_mult = 8e-8;
  mem_mult = 64;
  
else
  error('radar_name %s not supported yet.', radar_name);
  
end

%% Load data and create get_heights cluster tasks
% =====================================================================
%
% For each frame load REC_BLOCK_SIZE records at a time (code groups
% by file index, but has to watch negative offset values which imply
% the record starts in a previous file and carries over into the next)
%    --> The last block can range from 0.5 to 1.5 * REC_BLOCK_SIZE
% =====================================================================
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'get_heights_task';
sparam.num_args_out = 1;
sparam.argsin{1}.load.imgs = param.get_heights.imgs;
for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  
  % Check proc_mode from frames file that contains this frames type and
  % make sure the user has specified to process this frame type
  if ct_proc_frame(frames.proc_mode(frm),param.get_heights.frm_types)
    fprintf('%s %s_%03i (%i of %i) %s\n', mfilename, param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now,'HH:MM:SS'));
  else
    fprintf('Skipping %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end
  
  % Clean the array_proc temporary output directory for this frame
  sub_apt_shift_idx = 1;
  sub_band_idx = 1;
  out_fn_dir = fullfile(qlook_out_dir, ...
    sprintf('ql_data_%03d_%02d_%02d',frm,sub_apt_shift_idx,sub_band_idx));
  if ~ctrl.cluster.rerun_only && exist(out_fn_dir,'dir') 
    rmdir(out_fn_dir,'s');
  end

  % recs: Determine the records for this frame
  if frm < length(frames.frame_idxs)
    recs = frames.frame_idxs(frm):frames.frame_idxs(frm+1)-1;
  else
    recs = frames.frame_idxs(frm):length(records.gps_time);
  end
  
  % Determine where breaks in processing blocks are going to occur
  %   Rename variables for readability
  block_size = param.get_heights.block_size(1);
  block_overlap = param.get_heights.block_size(2);
  breaks = 1:block_size:length(recs)-0.5*block_size;
  
  % Create a cluster task for each block
  for break_idx = 1:length(breaks)
    
    % Determine the current records being processed
    % =================================================================
    if break_idx < length(breaks)
      cur_recs_keep = [recs(breaks(break_idx)) recs(breaks(break_idx+1)-1)];
      cur_recs = [max(1,recs(breaks(break_idx))-block_overlap) ...
        recs(breaks(break_idx+1)-1)+block_overlap];
    else
      cur_recs_keep = [recs(breaks(break_idx)) recs(end)];
      cur_recs = [max(1,recs(breaks(break_idx))-block_overlap) min(length(records.gps_time),recs(end)+block_overlap)];
    end
    
    % Fields required for manual submission to Slurm on Ollie
    if strcmp(param.cluster.type,'ollie')
      n_breaks(frm_idx) = length(breaks);
      dynamic_param.frms.(['frm',num2str(frm)]).frm_id = frm;
      dynamic_param.frms.(['frm',num2str(frm)]).breaks.(['break',num2str(break_idx)]).break_id = break_idx;
      dynamic_param.frms.(['frm',num2str(frm)]).breaks.(['break',num2str(break_idx)]).recs = cur_recs;
      dynamic_param.frms.(['frm',num2str(frm)]).breaks.(['break',num2str(break_idx)]).recs_keep = cur_recs_keep;
      continue;
    end
    
    % Prepare task inputs
    % =================================================================
    dparam = [];
    dparam.argsin{1}.load.frm = frm;
    dparam.argsin{1}.load.recs = cur_recs;
    dparam.argsin{1}.load.recs_keep = cur_recs_keep;
    % Set the Nyquist zone field (FMCW radars)
    wf = 1;
    if isfield(frames,'nyquist_zone') && ~isnan(frames.nyquist_zone(frm))
      dparam.argsin{1}.radar.wfs(wf).nyquist_zone = frames.nyquist_zone(frm);
    end
    
    % Create success condition
    % =================================================================
    dparam.success = '';
    for img = 1:length(param.get_heights.imgs)
      out_fn_name = sprintf('qlook_img_%02d_%d_%d.mat',img,cur_recs_keep(1),cur_recs_keep(end));
      out_fn{img} = fullfile(out_fn_dir,out_fn_name);
      if img == 1
        dparam.success = cat(2,dparam.success, ...
          sprintf('if ~exist(''%s'',''file'')', out_fn{img}));
      else
        dparam.success = cat(2,dparam.success, ...
          sprintf(' || ~exist(''%s'',''file'')', out_fn{img}));
      end
    end
    dparam.success = cat(2,dparam.success,sprintf('\n'));
    if 0
      % Enable this check if you want to open each output file to make
      % sure it is not corrupt.
      for img = 1:length(param.get_heights.imgs)
        out_fn_name = sprintf('qlook_img_%02d_%d_%d.mat',img,cur_recs_keep(1),cur_recs_keep(end));
        out_fn{img} = fullfile(out_fn_dir,out_fn_name);
        dparam.success = cat(2,dparam.success, ...
          sprintf('  load(''%s'');\n', out_fn{img}));
      end
    end
    success_error = 64;
    dparam.success = cat(2,dparam.success, ...
      sprintf('  error_mask = error_mask + %d;\n', success_error));
    dparam.success = cat(2,dparam.success,sprintf('end;\n'));
    
    % Rerun only mode: Test to see if we need to run this task
    % =================================================================
    if ctrl.cluster.rerun_only
      % If we are in rerun only mode AND the get heights task success
      % condition passes without error, then we do not run the task.
      error_mask = 0;
      eval(dparam.success);
      if ~error_mask
        fprintf('  %d: Already exists records %d to %d [rerun_only skipping] (%s)\n', ...
          break_idx, cur_recs(1), cur_recs(end), datestr(now));
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
    for img = 1:length(param.get_heights.imgs)
      dparam.cpu_time = dparam.cpu_time + 10 + Nx*total_num_sam(img)*log2(total_num_sam(img))*cpu_time_mult;
      dparam.mem = max(dparam.mem,250e6 + Nx*total_num_sam(img)*mem_mult);
    end
    dparam.notes = sprintf('%s:%s:%s %s_%03d (%d of %d)/%d of %d recs %d-%d', ...
      mfilename, param.radar_name, param.season_name, param.day_seg, frm, frm_idx, length(param.cmd.frms), ...
      break_idx, length(breaks), cur_recs_keep(1), cur_recs_keep(end));
    
    ctrl = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);
  end
end

% Export jobs files for manual submission to Slurm on Ollie
if strcmp(param.cluster.type,'ollie')
  dynamic_param.day_seg = param.day_seg;
  static_param = sparam.argsin{1};
  dynamic_param_file_name = sprintf('%s/qlook_%s_dynamic_param.mat', param.slurm_jobs_path, param.day_seg);
  save(dynamic_param_file_name,'dynamic_param');
  fprintf('Writing %s\n',dynamic_param_file_name);
  
  static_param_file_name = sprintf('%s/qlook_%s_static_param.mat', param.slurm_jobs_path, param.day_seg);
  save(static_param_file_name,'static_param');
  fprintf('Writing %s\n',static_param_file_name);
  
  txt_file_name = sprintf('%s/qlook_%s_parameters.txt', param.slurm_jobs_path, dynamic_param.day_seg);
  fid = fopen(txt_file_name,'w');
  fprintf(fid,'%3s\t %5s\n','frm','break');
  frms = fieldnames(dynamic_param.frms);
  for frm_idx = 1:length(param.cmd.frms)
    breaks = fieldnames(dynamic_param.frms.(frms{frm_idx}).breaks);
    for break_idx = 1:n_breaks(frm_idx)
      params = [dynamic_param.frms.(frms{frm_idx}).frm_id, dynamic_param.frms.(frms{frm_idx}).breaks.(breaks{break_idx}).break_id];
      formatSpec = '%03d\t %03d\n';
      fprintf(fid,formatSpec,params);
    end
  end
  fclose(fid);
  fprintf('Writing %s\n',txt_file_name);
  fprintf('Run batch_qlook.sh and batch_qlook_2.sh\n');
  
  ctrl_chain = {};
  return;
end

ctrl = cluster_save_dparam(ctrl);

ctrl_chain = {ctrl};


%% Create and setup the combine batch
% =====================================================================
ctrl = cluster_new_batch(param);
if param.get_heights.surf.en
  % If surface is enabled, the records file will be updated and this should
  % not be done on the cluster.
  ctrl.cluster.type = 'debug';
end

if any(strcmpi(radar_name,{'acords','hfrds','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  cpu_time_mult = 6e-8;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  cpu_time_mult = 100e-8;
  mem_mult = 24;
end

sparam = [];
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'get_heights_combine_task';
sparam.num_args_out = 1;
sparam.cpu_time = 10;
sparam.mem = 0;
% Add up all records being processed and find the most records in a frame
Nx = 0;
Nx_max = 0;
for frm = param.cmd.frms
  % recs: Determine the records for this frame
  if frm < length(frames.frame_idxs)
    Nx_frm = frames.frame_idxs(frm+1) - frames.frame_idxs(frm);
  else
    Nx_frm = length(records.gps_time) - frames.frame_idxs(frm) + 1;
  end
  if Nx_frm > Nx_max
    Nx_max = Nx_frm;
  end
  Nx = Nx + Nx_frm;
end
% Account for averaging
Nx_max = Nx_max / param.get_heights.decimate_factor / max(1,param.get_heights.inc_ave);
Nx = Nx / param.get_heights.decimate_factor / max(1,param.get_heights.inc_ave);
for img = 1:length(param.get_heights.imgs)
  sparam.cpu_time = sparam.cpu_time + (Nx*total_num_sam(img)*cpu_time_mult);
  if isempty(param.get_heights.img_comb)
    % Individual images, so need enough memory to hold the largest image
    sparam.mem = max(sparam.mem,250e6 + Nx_max*total_num_sam(img)*mem_mult);
  else
    % Images combined into one so need enough memory to hold all images
    sparam.mem = 250e6 + Nx*sum(total_num_sam)*mem_mult;
  end
end
if param.get_heights.surf.en
  sparam.cpu_time = sparam.cpu_time + numel(records.gps_time)/5e6*120;
end
sparam.notes = sprintf('%s:%s:%s %s combine frames', ...
  mfilename, param.radar_name, param.season_name, param.day_seg);

% Create success condition
success_error = 64;
sparam.success = '';
for frm = param.cmd.frms
  out_fn_name = sprintf('Data_%s_%03d.mat',param.day_seg,frm);
  out_fn = fullfile(qlook_out_dir,out_fn_name);
  sparam.success = cat(2,sparam.success, ...
    sprintf('  error_mask = bitor(error_mask,%d*(~exist(''%s'',''file'')));\n', success_error, out_fn));
end

ctrl = cluster_new_task(ctrl,sparam,[]);

ctrl_chain{end+1} = ctrl;
    
fprintf('Done %s\n', datestr(now));

return;

