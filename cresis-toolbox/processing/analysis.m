function ctrl_chain = analysis(param,param_override)
% ctrl_chain = analysis(param,param_override)
%
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Example:
%  See run_analysis.m for how to run this function directly.
%  Normally this function is called from master.m using the param spreadsheet.
%
% Authors: John Paden
%
% See also: master.m, run_analysis.m analysis.m,
%   analysis_task.m

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks
% =====================================================================

if ~isempty(param.cmd.frms)
  warning('All frames are always processed with analysis, setting param.cmd.frms to do all frames.');
  param.cmd.frms = []; % All frames
end

if ~isfield(param.analysis,'out_path') || isempty(param.analysis.out_path)
  param.analysis.out_path = 'analysis';
end

if ~isfield(param.analysis,'block_size') || isempty(param.analysis.block_size)
  error('param.analysis.block_size must be specified');
end

if ~isfield(param.analysis,'imgs') || isempty(param.analysis.imgs)
  error('param.analysis.imgs must be specified');
end

if ~isfield(param.analysis,'presums') || isempty(param.analysis.presums)
  param.analysis.presums = 1;
end

if ~isfield(param.analysis,'lever_arm_fh') || isempty(param.analysis.lever_arm_fh)
  param.analysis.lever_arm_fh = [];
end

if ~isfield(param.analysis,'trim_vals') || isempty(param.analysis.trim_vals)
  param.analysis.trim_vals = [0 0];
end

% For each command in the list, set its default settings
for cmd_idx = 1:length(param.analysis.cmd)
  cmd = param.analysis.cmd{cmd_idx};
  
  if ~isfield(cmd,'en') || isempty(cmd.en)
    cmd.en = true;
  end
  
  if ~isfield(cmd,'out_path') || isempty(cmd.out_path)
    cmd.out_path = param.analysis.out_path;
  end
  
  if ~isfield(cmd,'wf_adc_idxs') || isempty(cmd.wf_adc_idxs)
    for img = 1:length(param.analysis.imgs)
      cmd.wf_adc_idxs{img} = [];
    end
  end
  
  for img = 1:length(param.analysis.imgs)
    if length(cmd.wf_adc_idxs) < length(param.analysis.imgs) || isempty(cmd.wf_adc_idxs{img})
      % By default do all wf-adc pairs in the image
      cmd.wf_adc_idxs{img} = 1:size(param.analysis.imgs{img},1);
    end
  end
  
  if ~isfield(cmd,'layer') || isempty(cmd.layer)
    % Set the analysis start time to the beginning of the record
    cmd.layer = -inf;
  end
  
  if ~isfield(cmd,'Nt') || isempty(cmd.Nt)
    % Set the analysis length to the whole record
    cmd.Nt = inf;
  end
  
  if ~isfield(cmd,'B_filter') || isempty(cmd.B_filter)
    % Set the default filter to no filtering (B_filter = 1)
    cmd.B_filter = 1;
  end
  
  if abs(sum(cmd.B_filter)-1) > 1e4*eps
    %warning('B_filter weights are not normalized. They must be normalized so normalizing to one now.')
    cmd.B_filter = cmd.B_filter / sum(cmd.B_filter);
  end
  
  if ~isfield(cmd,'decimate_factor') || isempty(cmd.decimate_factor)
    % Set the default decimation to none (decimate_factor = 1)
    cmd.decimate_factor = 1;
  end

  switch lower(cmd.name)
    case {'qlook'}
      %
    case {'waveform'}
      %
    case {'statistics'}
      %
    case {'saturation'}
      %
    case {'specular'}
      %
    case {'coherent_noise'}
      % Set defaults fro coherent noise analysis method
      
      if ~isfield(cmd,'block_ave') || isempty(cmd.block_ave)
        % Set the default block_ave to 1000
        cmd.block_ave = 1000;
      end
      
      if mod(param.analysis.block_size,cmd.block_ave)
        error('The param.analysis.block_size (%s) must be a multiple of cmd.block_ave (%d).', ...
          param.analysis.block_size, cmd.block_ave);
      end
      
      if ~isfield(cmd,'power_threshold') || isempty(cmd.power_threshold)
        % Set the default power_threshold to inf (i.e. no thresholding)
        cmd.power_threshold = inf;
      end
      
      cmd.coh_noise_method = 0;
      cmd.coh_noise_arg = 0;
      
    case {'burst_noise'}
      %
  end
  
  % Update the command structure
  param.analysis.cmd{cmd_idx} = cmd;
end

%% Setup processing
% =====================================================================

% Get the standard radar name
[~,~,radar_name] = ct_output_dir(param.radar_name);

% Load frames file
load(ct_filename_support(param,'','frames'));

% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);
% Apply presumming
if param.analysis.presums > 1
  records.lat = fir_dec(records.lat,param.analysis.presums);
  records.lon = fir_dec(records.lon,param.analysis.presums);
  records.elev = fir_dec(records.elev,param.analysis.presums);
  records.roll = fir_dec(records.roll,param.analysis.presums);
  records.pitch = fir_dec(records.pitch,param.analysis.presums);
  records.heading = fir_dec(records.heading,param.analysis.presums);
  records.gps_time = fir_dec(records.gps_time,param.analysis.presums);
  records.surface = fir_dec(records.surface,param.analysis.presums);
end

% Compute all estimates with pulse compressed numbers even though raw data
% are loaded
param.analysis.pulse_comp = true;
param.analysis.ft_wind = [];

ctrl_chain = {};

%% Create and setup the cluster batch
% =====================================================================
ctrl = cluster_new_batch(param);
cluster_compile({'analysis_task.m','analysis_combine_task.m'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);

if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  [wfs,~] = load_mcords_wfs(records.settings, param, ...
    1:max(records.param_records.records.file.adcs), param.analysis);
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
out_fn_dir = ct_filename_out(param,param.analysis.out_path);

sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'analysis_task';
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
  actual_cur_recs = [(cur_recs(1)-1)*param.csarp.presums+1, ...
    cur_recs(end)*param.csarp.presums];
  
  % Prepare task inputs
  % =================================================================
  dparam = [];
  dparam.argsin{1}.load.recs = cur_recs;
  
  % Create success condition and set cpu_time, mem requirements
  % =================================================================
  Nx = cur_recs(end)-cur_recs(1)+1;
  dparam.cpu_time = 0;
  dparam.mem = 0;
  dparam.success = '';
  success_error = 64;
  for img = 1:length(param.analysis.imgs)
    for cmd_idx = 1:length(param.analysis.cmd)
      cmd = param.analysis.cmd{cmd_idx}; % cmd: current command
      if ~cmd.en
        continue;
      end
      
      % Create output directory string
      out_fn_dir = ct_filename_out(param,cmd.out_path);
      
      % Load data
      dparam.cpu_time = dparam.cpu_time + 10 + param.csarp.presums*size(param.analysis.imgs{img},1)*Nx*total_num_sam(img)*log2(total_num_sam(img))*cpu_time_mult;

      % Process commands
      switch lower(cmd.name)
        case {'qlook'}
          %
        case {'waveform'}
          out_fn = fullfile(out_fn_dir,sprintf('waveform_img_%02d_%d_%d.mat',img,actual_cur_recs));
          dparam.success = cat(2,dparam.success, ...
            sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
          if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
            delete(out_fn);
          end
          dparam.cpu_time = dparam.cpu_time + 10 + size(param.analysis.imgs{img},1)*Nx*total_num_sam(img)*log2(total_num_sam(img))*cpu_time_mult;
          dparam.mem = max(dparam.mem,250e6 + Nx*total_num_sam(img)*mem_mult);
          
        case {'statistics'}
          %
        case {'saturation'}
          %
        case {'specular'}
          out_fn = fullfile(out_fn_dir,sprintf('specular_img_%02d_%d_%d.mat',img,actual_cur_recs));
          dparam.success = cat(2,dparam.success, ...
            sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
          if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
            delete(out_fn);
          end
          dparam.cpu_time = dparam.cpu_time + 10 + size(param.analysis.imgs{img},1)*Nx*total_num_sam(img)*log2(total_num_sam(img))*cpu_time_mult;
          dparam.mem = max(dparam.mem,250e6 + Nx*total_num_sam(img)*mem_mult);
          
        case {'coherent_noise'}
          out_fn = fullfile(out_fn_dir,sprintf('coh_noise_img_%02d_%d_%d.mat',img,actual_cur_recs));
          dparam.success = cat(2,dparam.success, ...
            sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
          if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
            delete(out_fn);
          end
          dparam.cpu_time = dparam.cpu_time + 10 + size(param.analysis.imgs{img},1)*Nx*total_num_sam(img)*log2(Nx)*cpu_time_mult;
          dparam.mem = max(dparam.mem,250e6 + Nx*total_num_sam(img)*mem_mult);
          
        case {'burst_noise'}
          %
      end
    end
  end
  
  % Rerun only mode: Test to see if we need to run this task
  % =================================================================
  dparam.notes = sprintf('%s:%s:%s %s %d of %d recs %d-%d', ...
    mfilename, param.radar_name, param.season_name, param.day_seg, ...
    break_idx, length(breaks), actual_cur_recs);
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
  ctrl = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);
  
end

ctrl = cluster_save_dparam(ctrl);

ctrl_chain = {ctrl};

%% Create and setup the combine batch
% =====================================================================
ctrl = cluster_new_batch(param);

if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  cpu_time_mult = 6e-6;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  cpu_time_mult = 1000e-8;
  mem_mult = 24;
end

% Create success condition
success_error = 64;
out_fn_dir_dir = fileparts(out_fn_dir);
sparam = [];
sparam.success = '';
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'analysis_combine_task';
sparam.num_args_out = 1;
sparam.cpu_time = 60;
sparam.mem = 0;
% Add up all records being processed and find the most records in a block
Nx = length(records.gps_time);
for img = 1:length(param.analysis.imgs)
  Nt = total_num_sam(img);
  
  for cmd_idx = 1:length(param.analysis.cmd)
    cmd = param.analysis.cmd{cmd_idx};
    if ~cmd.en
      continue;
    end
    
    switch lower(cmd.name)
      case {'qlook'}
        %
      case {'waveform'}
        Nx_cmd = Nx / param.get_heights.decimate_factor;
        if isfinite(param.analysis.surf.Nt)
          Nt = param.analysis.surf.Nt;
        end
        sparam.cpu_time = sparam.cpu_time + size(param.analysis.imgs{img},1)*Nx_cmd*Nt*cpu_time_mult;
        sparam.mem = max(sparam.mem,250e6 + size(param.analysis.imgs{img},1)*Nx_cmd*Nt*mem_mult);
        out_fn = fullfile(out_fn_dir_dir,sprintf('surf_img_%02d.mat',img));
        dparam.success = cat(2,dparam.success, ...
          sprintf('  error_mask = bitor(error_mask,%d*~exist(%s,''file''));\n', success_error, out_fn));
        if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
          delete(out_fn);
        end
        
      case {'statistics'}
        %
      case {'saturation'}
        %
      case {'specular'}
        Nx_cmd = Nx / param.analysis.block_size * param.analysis.specular.threshold_max;
        sparam.cpu_time = sparam.cpu_time + size(param.analysis.imgs{img},1)*Nx_cmd*Nt*cpu_time_mult;
        sparam.mem = max(sparam.mem,250e6 + size(param.analysis.imgs{img},1)*Nx_cmd*Nt*mem_mult);
        out_fn = fullfile(out_fn_dir_dir,sprintf('specular_img_%02d.mat',img));
        dparam.success = cat(2,dparam.success, ...
          sprintf('  error_mask = bitor(error_mask,%d*~exist(%s,''file''));\n', success_error, out_fn));
        if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
          delete(out_fn);
        end
        
      case {'coherent_noise'}
        Nx_cmd = Nx / cmd.block_ave;
        sparam.cpu_time = sparam.cpu_time + size(param.analysis.imgs{img},1)*Nx_cmd*Nt*cpu_time_mult;
        sparam.mem = max(sparam.mem,250e6 + size(param.analysis.imgs{img},1)*Nx_cmd*Nt*mem_mult);
        out_fn = fullfile(out_fn_dir_dir,sprintf('coh_noise_img_%02d.mat',img));
        dparam.success = cat(2,dparam.success, ...
          sprintf('  error_mask = bitor(error_mask,%d*~exist(%s,''file''));\n', success_error, out_fn));
        if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
          delete(out_fn);
        end
        
      case {'burst_noise'}
        %
    end
  end
end
sparam.notes = sprintf('%s:%s:%s %s combine', ...
  mfilename, param.radar_name, param.season_name, param.day_seg);

ctrl = cluster_new_task(ctrl,sparam,[]);

ctrl_chain{end+1} = ctrl;
    
fprintf('Done %s\n', datestr(now));

return
