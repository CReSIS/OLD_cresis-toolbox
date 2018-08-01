function ctrl_chain = qlook(param,param_override)
% ctrl_chain = qlook(param,param_override)
%
% This function generates quick look outputs (CSARP_qlook), tracks the
% surface, and (optionally) stores the surface to a layer data destination
% (default is layerData).
%
% param = struct with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Example:
%  See run_qlook.m for how to run this function directly.
%  Normally this function is called from master.m using the param spreadsheet.
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_qlook.m, qlook.m,
%   qlook_task.m

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

if ~isfield(param.qlook,'frm_types') || isempty(param.qlook.frm_types)
  param.qlook.frm_types = {-1,-1,-1,-1,-1};
end

if ~isfield(param.qlook,'out_path') || isempty(param.qlook.out_path)
  param.qlook.out_path = 'qlook';
end

if ~isfield(param.qlook,'block_size') || isempty(param.qlook.block_size)
  error('param.qlook.block_size must be specified. This is the number of range lines or records to process at a time.');
end

if ~isfield(param.qlook,'presums') || isempty(param.qlook.presums)
  param.qlook.presums = 1;
end

if ~isfield(param.qlook,'imgs') || isempty(param.qlook.imgs)
  param.qlook.imgs = {[1 1]};
end
% block_size: no default

if ~isfield(param.qlook,'imgs') || isempty(param.qlook.imgs)
  param.qlook.imgs = {[1 1]};
end

if ~isfield(param.qlook,'img_comb') || isempty(param.qlook.img_comb)
  param.qlook.img_comb = [];
end

% Check img_comb length
if ~isempty(param.qlook.img_comb) && length(param.qlook.img_comb) ~= 3*(length(param.qlook.imgs)-1)
  error('param.qlook.img_comb not the right length. Since it is not empty, there should be 3 entries for each image combination interface ([Tpd second image for surface saturation, -inf for second image blank, Tpd first image to avoid roll off] is typical).');
end

if ~isfield(param.qlook,'motion_comp') || isempty(param.qlook.motion_comp)
  param.qlook.motion_comp = false;
end

if ~isfield(param.qlook,'dec') || isempty(param.qlook.dec)
  param.qlook.dec = 1;
end
if ~isfield(param.qlook,'B_filter') || isempty(param.qlook.B_filter)
  if param.qlook.dec == 1
    param.qlook.B_filter = 1;
  else
    param.qlook.B_filter = hanning(2*param.qlook.dec+1);
  end
end
if ~mod(length(param.qlook.B_filter),2)
  error('param.qlook.B_filter must be odd length.');
end
param.qlook.B_filter = param.qlook.B_filter(:).'; % Must be row vector
if abs(sum(param.qlook.B_filter)-1) > 1e4*eps
  param.qlook.B_filter = param.qlook.B_filter / sum(param.qlook.B_filter);
end

if ~isfield(param.qlook,'inc_dec') || isempty(param.qlook.inc_dec)
  param.qlook.inc_dec = 1;
end
if ~isfield(param.qlook,'inc_B_filter') || isempty(param.qlook.inc_B_filter)
  if param.qlook.inc_dec == 0 || param.qlook.inc_dec == 1
    param.qlook.inc_B_filter = 1;
  else
    param.qlook.inc_B_filter = hanning(2*param.qlook.inc_dec+1);
  end
end
if ~mod(length(param.qlook.inc_B_filter),2)
  error('param.qlook.inc_B_filter must be odd length.');
end
param.qlook.inc_B_filter = param.qlook.inc_B_filter(:).'; % Must be row vector
if abs(sum(param.qlook.inc_B_filter)-1) > 1e4*eps
  param.qlook.inc_B_filter = param.qlook.inc_B_filter / sum(param.qlook.inc_B_filter);
end

if ~isfield(param.qlook,'resample') || isempty(param.qlook.resample)
  param.qlook.resample = [1 1; 1 1];
end

if ~isfield(param.qlook,'surf') || isempty(param.qlook.surf)
  param.qlook.surf.en = false;
end

if ~isfield(param.qlook,'top_layer') || isempty(param.qlook.top_layer)
  param.qlook.top_layer = [];
end

%% Setup Processing
% =====================================================================

% Get the standard radar name
[~,~,radar_name] = ct_output_dir(param.radar_name);

% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);

% Quick look radar echogram output directory
qlook_out_dir = ct_filename_out(param, param.qlook.out_path);

% % Get version information out of the deconvolution file
% if isfield(param.qlook,'deconvolution') ...
%     && ~isempty(param.qlook.deconvolution) ...
%     && param.qlook.deconvolution == 3
%   out_fn_dir = ct_filename_out(param,'analysis');
%   out_segment_fn_dir = fileparts(out_fn_dir);
%   out_segment_fn = fullfile(out_segment_fn_dir,sprintf('deconv_%s.mat', param.day_seg));
%   spec = load(out_segment_fn,'param_collate');
%   
%   param.qlook.deconvolution_sw_version = spec.param_collate.sw_version;
%   param.qlook.deconvolution_params = spec.param_collate.analysis.specular;
% end
% 
% % Get version information out of the coherent noise file
% if any(param.qlook.coh_noise_method == [17 19])
%   
%   cdf_fn_dir = fileparts(ct_filename_out(param,param.qlook.coh_noise_arg{4}, ''));
%   cdf_fn = fullfile(cdf_fn_dir,sprintf('coh_noise_simp_%s.nc', param.day_seg));
%   
%   tmp = netcdf_to_mat(cdf_fn,[],'^sw_version.*');
%   param.qlook.coh_noise_version = tmp.sw_version;
%   tmp = netcdf_to_mat(cdf_fn,[],'^param_collate.*');
%   param.qlook.coh_noise_params = tmp.param_collate;
% end

%% Create and setup the cluster batch
% =====================================================================
ctrl = cluster_new_batch(param);
cluster_compile({'qlook_task.m','qlook_combine_task.m'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);

total_num_sam = [];
if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','mcrds','seaice','accum2'}))
  [wfs,~] = load_mcords_wfs(records.settings, param, ...
    1:max(records.param_records.records.file.adcs), param.qlook);
  for img = 1:length(param.qlook.imgs)
    wf = abs(param.qlook.imgs{img}(1,1));
    total_num_sam(img) = wfs(wf).Nt_raw;
  end
  cpu_time_mult = 66e-8;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'mcrds'}))
  [wfs,~] = load_mcrds_wfs(records.settings, param, ...
    1:max(records.param_records.records.file.adcs), param.qlook);
  for img = 1:length(param.qlook.imgs)
    wf = abs(param.qlook.imgs{img}(1,1));
    total_num_sam(img) = wfs(wf).Nt_raw;
  end
  cpu_time_mult = 66e-8;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  total_num_sam = 32000 * ones(size(param.qlook.imgs));
  cpu_time_mult = 8e-8;
  mem_mult = 64;
  
else
  error('radar_name %s not supported yet.', radar_name);
  
end

%% Load data and create qlook cluster tasks
% =====================================================================
%
% For each frame load REC_BLOCK_SIZE records at a time (code groups
% by file index, but has to watch negative offset values which imply
% the record starts in a previous file and carries over into the next)
%    --> The last block can range from 0.5 to 1.5 * REC_BLOCK_SIZE
% =====================================================================
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'qlook_task';
sparam.num_args_out = 1;
sparam.argsin{1}.load.imgs = param.qlook.imgs;
for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  
  % Check proc_mode from frames file that contains this frames type and
  % make sure the user has specified to process this frame type
  if ct_proc_frame(frames.proc_mode(frm),param.qlook.frm_types)
    fprintf('%s %s_%03i (%i of %i) (%s)\n', sparam.task_function, param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now));
  else
    fprintf('Skipping %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end
  
  % Create output directory name
  sub_apt_shift_idx = 1;
  sub_band_idx = 1;
  out_fn_dir = fullfile(qlook_out_dir, ...
    sprintf('ql_data_%03d_%02d_%02d',frm,sub_apt_shift_idx,sub_band_idx));

  % recs: Determine the records for this frame
  if frm < length(frames.frame_idxs)
    recs = frames.frame_idxs(frm):frames.frame_idxs(frm+1)-1;
  else
    recs = frames.frame_idxs(frm):length(records.gps_time);
  end
  
  % Determine where breaks in processing blocks are going to occur
  %   Rename variables for readability
  block_size = param.qlook.block_size(1);
  blocks = 1:block_size:length(recs)-0.5*block_size;
  
  % Create a cluster task for each block
  for block_idx = 1:length(blocks)
    
    % Determine the current records being processed
    % =================================================================
    if block_idx < length(blocks)
      cur_recs = [recs(blocks(block_idx)) recs(blocks(block_idx+1)-1)];
    else
      cur_recs = [recs(blocks(block_idx)) recs(end)];
    end
    
    % Fields required for manual submission to Slurm on Ollie
    if strcmp(param.cluster.type,'ollie')
      n_blocks(frm_idx) = length(blocks);
      dynamic_param.frms.(['frm',num2str(frm)]).frm_id = frm;
      dynamic_param.frms.(['frm',num2str(frm)]).blocks.(['block',num2str(block_idx)]).block_id = block_idx;
      dynamic_param.frms.(['frm',num2str(frm)]).blocks.(['block',num2str(block_idx)]).recs = cur_recs;
      dynamic_param.frms.(['frm',num2str(frm)]).blocks.(['block',num2str(block_idx)]).recs_keep = cur_recs;
      continue;
    end
    
    % Prepare task inputs
    % =================================================================
    dparam = [];
    dparam.argsin{1}.load.frm = frm;
    dparam.argsin{1}.load.recs = cur_recs;
    % Set the Nyquist zone field (FMCW radars)
    wf = 1;
    if isfield(frames,'nyquist_zone') && ~isnan(frames.nyquist_zone(frm))
      dparam.argsin{1}.radar.wfs(wf).nyquist_zone = frames.nyquist_zone(frm);
    end
    
    % Create success condition
    % =================================================================
    dparam.success = '';
    for img = 1:length(param.qlook.imgs)
      out_fn_name = sprintf('qlook_img_%02d_%d_%d.mat',img,cur_recs(1),cur_recs(end));
      out_fn{img} = fullfile(out_fn_dir,out_fn_name);
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
      for img = 1:length(param.qlook.imgs)
        out_fn_name = sprintf('qlook_img_%02d_%d_%d.mat',img,cur_recs(1),cur_recs(end));
        out_fn{img} = fullfile(out_fn_dir,out_fn_name);
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
    dparam.notes = sprintf('%s:%s:%s %s_%03d (%d of %d)/%d of %d recs %d-%d', ...
      sparam.task_function, param.radar_name, param.season_name, param.day_seg, frm, frm_idx, length(param.cmd.frms), ...
      block_idx, length(blocks), cur_recs(1), cur_recs(end));
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
    for img = 1:length(param.qlook.imgs)
      dparam.cpu_time = dparam.cpu_time + 10 + Nx*size(param.qlook.imgs{img},1)*total_num_sam(img)*log2(total_num_sam(img))*cpu_time_mult;
      dparam.mem = max(dparam.mem,250e6 + Nx*size(param.qlook.imgs{img},1)*total_num_sam(img)*mem_mult);
    end
    
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
    blocks = fieldnames(dynamic_param.frms.(frms{frm_idx}).blocks);
    for block_idx = 1:n_blocks(frm_idx)
      params = [dynamic_param.frms.(frms{frm_idx}).frm_id, dynamic_param.frms.(frms{frm_idx}).blocks.(blocks{block_idx}).block_id];
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
if param.qlook.surf.en
  % If surface is enabled, the records file will be updated and this should
  % not be done on the cluster.
  ctrl.cluster.type = 'debug';
end

if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','mcrds','seaice','accum2'}))
  cpu_time_mult = 6e-8;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  cpu_time_mult = 100e-8;
  mem_mult = 24;
end

sparam = [];
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'qlook_combine_task';
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
Nx_max = Nx_max / param.qlook.dec / max(1,param.qlook.inc_dec);
Nx = Nx / param.qlook.dec / max(1,param.qlook.inc_dec);
for img = 1:length(param.qlook.imgs)
  sparam.cpu_time = sparam.cpu_time + (Nx*total_num_sam(img)*cpu_time_mult);
  if isempty(param.qlook.img_comb)
    % Individual images, so need enough memory to hold the largest image
    sparam.mem = max(sparam.mem,250e6 + Nx_max*total_num_sam(img)*mem_mult);
  else
    % Images combined into one so need enough memory to hold all images
    sparam.mem = 250e6 + Nx*sum(total_num_sam)*mem_mult;
  end
end
if param.qlook.surf.en
  sparam.cpu_time = sparam.cpu_time + numel(records.gps_time)/5e6*120;
end
sparam.notes = sprintf('%s:%s:%s %s', ...
  sparam.task_function, param.radar_name, param.season_name, param.day_seg);

% Create success condition
success_error = 64;
sparam.success = '';
for frm = param.cmd.frms
  out_fn_name = sprintf('Data_%s_%03d.mat',param.day_seg,frm);
  out_fn = fullfile(qlook_out_dir,out_fn_name);
  sparam.success = cat(2,sparam.success, ...
    sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
  if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
    ct_file_lock_check(out_fn);
    delete(out_fn);
  end
end

ctrl = cluster_new_task(ctrl,sparam,[]);

ctrl_chain{end+1} = ctrl;
    
fprintf('Done %s\n', datestr(now));

return;

