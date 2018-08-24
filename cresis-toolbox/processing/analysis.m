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

if ~isfield(param.analysis,'surf_layer') || isempty(param.analysis.surf_layer)
  param.analysis.surf_layer.name = 'surface';
  param.analysis.surf_layer.source = 'layerData';
end
% Never check for the existence of files
param.analysis.surf_layer.existence_check = false;

% For each command in the list, set its default settings
for cmd_idx = 1:length(param.analysis.cmd)
  cmd = param.analysis.cmd{cmd_idx};
  
  if ~isfield(cmd,'en') || isempty(cmd.en)
    cmd.en = true;
  end
  
  if ~isfield(cmd,'out_path') || isempty(cmd.out_path)
    cmd.out_path = param.analysis.out_path;
  end
  
  if ~isfield(cmd,'wf_adcs') || isempty(cmd.wf_adcs)
    for img = 1:length(param.analysis.imgs)
      % By default do all wf-adc pairs in the image
      cmd.wf_adcs{img} = 1:size(param.analysis.imgs{img},1);
    end
  end
  
  for img = 1:length(param.analysis.imgs)
    if length(cmd.wf_adcs) < length(param.analysis.imgs)
      % By default do all wf-adc pairs in the image
      cmd.wf_adcs{img} = 1:size(param.analysis.imgs{img},1);
    end
  end
  
  if ~isfield(cmd,'start_time') || isempty(cmd.start_time)
    % Set the analysis start time to the beginning of the record
    cmd.start_time = -inf;
  end
  if isstruct(cmd.start_time)
    if ~isfield(cmd.start_time,'name') || isempty(cmd.start_time.name)
      cmd.start_time.name = 'surface';
    end
    if ~isfield(cmd.start_time,'source') || isempty(cmd.start_time.source)
      cmd.start_time.source = 'layerData';
    end
  end
  
  if ~isfield(cmd,'stop_time') || isempty(cmd.stop_time)
    % Set the analysis stop time to the end of the record
    cmd.stop_time = inf;
  end
  if isstruct(cmd.stop_time)
    if ~isfield(cmd.stop_time,'name') || isempty(cmd.stop_time.name)
      cmd.stop_time.name = 'surface';
    end
    if ~isfield(cmd.stop_time,'source') || isempty(cmd.stop_time.source)
      cmd.stop_time.source = 'layerData';
    end
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

  if ~isfield(cmd,'method') || isempty(cmd.method)
    error('cmd.method must be defined in param.analysis.cmd cell array');
  end
  
  switch lower(cmd.method)
    case {'coh_noise'}
      % Set defaults for coherent noise analysis method
      
      if ~isfield(cmd,'block_ave') || isempty(cmd.block_ave)
        cmd.block_ave = 2000;
      end
      
      if mod(param.analysis.block_size,cmd.block_ave)
        error('The param.analysis.block_size (%s) must be a multiple of cmd.block_ave (%d).', ...
          param.analysis.block_size, cmd.block_ave);
      end
      
      if ~isfield(cmd,'power_threshold') || isempty(cmd.power_threshold)
        % Set the default power_threshold to inf (i.e. no thresholding)
        cmd.power_threshold = inf;
      end
      
    case {'burst_noise'}
      %
    case {'qlook'}
      %
    case {'saturation'}
      %
    case {'specular'}
      % Set defaults for specular analysis method
      
      if ~isfield(cmd,'rlines') || isempty(cmd.rlines)
        cmd.rlines = 128;
      end
      
      if ~isfield(cmd,'max_rlines') || isempty(cmd.max_rlines)
        cmd.max_rlines = 10;
      end
      
      if ~isfield(cmd,'threshold') || isempty(cmd.threshold)
        cmd.threshold = 40;
      end
      
      if ~isfield(cmd,'signal_doppler_bins') || isempty(cmd.signal_doppler_bins)
        cmd.signal_doppler_bins = [1:4 cmd.rlines+(-3:0)];
      end
      
      if ~isfield(cmd,'noise_doppler_bins') || isempty(cmd.noise_doppler_bins)
        cmd.noise_doppler_bins = [12:cmd.rlines-11];
      end
      
      if ~isfield(cmd,'gps_times') || isempty(cmd.gps_times)
        cmd.gps_times = [];
      end
      
    case {'statistics'}
      % Set defaults for statistical analysis method
      
      if ~isfield(cmd,'block_ave') || isempty(cmd.block_ave)
        cmd.block_ave = 2000;
      end
      
      if ~isfield(cmd,'pulse_compress') || isempty(cmd.pulse_compress)
        cmd.pulse_compress = false;
      end
      
      if ~isfield(cmd,'motion_comp') || isempty(cmd.motion_comp)
        cmd.motion_comp = false;
      end
      
      if ~isfield(cmd,'stats') || isempty(cmd.stats)
        error('The statistical command requires that the stats field be set.');
      end
      
    case {'waveform'}
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

[wfs,~] = data_load_wfs(setfield(param,'load',struct('imgs',{param.analysis.imgs})),records);
if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','mcrds','seaice','accum2','accum3'}))
  for img = 1:length(param.analysis.imgs)
    wf = abs(param.analysis.imgs{img}(1,1));
    total_num_sam(img) = wfs(wf).Nt_raw;
  end
  cpu_time_mult = 140e-9;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  total_num_sam = 32000 * ones(size(param.analysis.imgs));
  cpu_time_mult = 4e-8;
  mem_mult = 17;
  
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
  actual_cur_recs = [(cur_recs(1)-1)*param.analysis.presums+1, ...
    cur_recs(end)*param.analysis.presums];
  
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
  % Loading in the data: cpu_time and mem
  dparam.mem = 250e6;
  for img = 1:length(param.analysis.imgs)
    dparam.cpu_time = dparam.cpu_time + 10 + length(param.analysis.imgs{img})*Nx*total_num_sam(img)*log2(Nx)*cpu_time_mult;
    dparam.mem = dparam.mem + length(param.analysis.imgs{img})*Nx*total_num_sam(img)*mem_mult;
  end
  data_load_memory = dparam.mem;
  % Processing the data
  for img = 1:length(param.analysis.imgs)
    for cmd_idx = 1:length(param.analysis.cmd)
      cmd = param.analysis.cmd{cmd_idx}; % cmd: current command
      if ~cmd.en
        continue;
      end
      
      % Create output directory string
      out_fn_dir = ct_filename_out(param,cmd.out_path);
      
      % Load data
      dparam.cpu_time = dparam.cpu_time + 10 + param.analysis.presums*size(param.analysis.imgs{img},1)*Nx*total_num_sam(img)*log2(total_num_sam(img))*cpu_time_mult;

      % Process commands
      switch lower(cmd.method)
        case {'burst_noise'}
          %

        case {'coh_noise'}
          for wf_adc = param.analysis.cmd{cmd_idx}.wf_adcs{img}(:).'
            wf = param.analysis.imgs{img}(wf_adc,1);
            adc = param.analysis.imgs{img}(wf_adc,2);
            out_fn = fullfile(out_fn_dir,sprintf('coh_noise_wf_%d_adc_%d_%d_%d.mat',wf,adc,actual_cur_recs));
            dparam.success = cat(2,dparam.success, ...
              sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
            if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
              delete(out_fn);
            end
            dparam.cpu_time = dparam.cpu_time + 10 + Nx*total_num_sam(img)*log2(Nx)*cpu_time_mult;
            dparam.mem = max(dparam.mem,data_load_memory + Nx*total_num_sam(img)*mem_mult);
          end
          
        case {'qlook'}
          %
          
        case {'saturation'}
          %
          
        case {'specular'}
          for wf_adc = param.analysis.cmd{cmd_idx}.wf_adcs{img}(:).'
            wf = param.analysis.imgs{img}(wf_adc,1);
            adc = param.analysis.imgs{img}(wf_adc,2);
            out_fn = fullfile(out_fn_dir,sprintf('specular_wf_%d_adc_%d_%d_%d.mat',wf,adc,actual_cur_recs));
            dparam.success = cat(2,dparam.success, ...
              sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
            if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
              delete(out_fn);
            end
            dparam.cpu_time = dparam.cpu_time + 10 + Nx*total_num_sam(img)*log2(total_num_sam(img))*cpu_time_mult;
            dparam.mem = max(dparam.mem,data_load_memory + Nx*total_num_sam(img)*mem_mult);
          end
          
        case {'statistics'}
          for wf_adc = param.analysis.cmd{cmd_idx}.wf_adcs{img}(:).'
            wf = param.analysis.imgs{img}(wf_adc,1);
            adc = param.analysis.imgs{img}(wf_adc,2);
            out_fn = fullfile(out_fn_dir,sprintf('stats_wf_%d_adc_%d_%d_%d.mat',wf,adc,actual_cur_recs));
            dparam.success = cat(2,dparam.success, ...
              sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
            if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
              delete(out_fn);
            end
            dparam.cpu_time = dparam.cpu_time + 10 + Nx*total_num_sam(img)*log2(Nx)*cpu_time_mult;
            dparam.mem = max(dparam.mem,data_load_memory + Nx*total_num_sam(img)*mem_mult);
          end
          
        case {'waveform'}
          for wf_adc = param.analysis.cmd{cmd_idx}.wf_adcs{img}(:).'
            wf = param.analysis.imgs{img}(wf_adc,1);
            adc = param.analysis.imgs{img}(wf_adc,2);
            out_fn = fullfile(out_fn_dir,sprintf('waveform_wf_%d_adc_%d_%d_%d.mat',wf,adc,actual_cur_recs));
            dparam.success = cat(2,dparam.success, ...
              sprintf('  error_mask = bitor(error_mask,%d*~exist(''%s'',''file''));\n', success_error, out_fn));
            if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
              delete(out_fn);
            end
            dparam.cpu_time = dparam.cpu_time + 10 + Nx*total_num_sam(img)*log2(total_num_sam(img))*cpu_time_mult;
            dparam.mem = max(dparam.mem,data_load_memory + Nx*total_num_sam(img)*mem_mult);
          end
          
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

if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','mcrds','seaice','accum2','accum3'}))
  cpu_time_mult = 6e-6;
  mem_mult = 8;
  
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
  cpu_time_mult = 2.6e-7;
  mem_mult = 32;
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
    
    switch lower(cmd.method)
      case {'burst_noise'}
        %
        
      case {'coh_noise'}
        Nx_cmd = Nx / cmd.block_ave;
        for wf_adc = param.analysis.cmd{cmd_idx}.wf_adcs{img}(:).'
          wf = param.analysis.imgs{img}(wf_adc,1);
          adc = param.analysis.imgs{img}(wf_adc,2);
          out_fn = fullfile(out_fn_dir_dir,sprintf('coh_noise_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc));
          sparam.success = cat(2,sparam.success, ...
            sprintf('  error_mask = bitor(error_mask,%d*~ct_file_lock_check(''%s'',4));\n', success_error, out_fn));
          if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
            ct_file_lock_check(out_fn,3);
          end
          sparam.cpu_time = sparam.cpu_time + Nx_cmd*total_num_sam(img)*log2(Nx_cmd)*cpu_time_mult;
          sparam.mem = max(sparam.mem,250e6 + Nx_cmd*total_num_sam(img)*mem_mult);
        end
        
      case {'qlook'}
        %
        
      case {'saturation'}
        %
        
      case {'specular'}
        Nx_cmd = Nx / param.analysis.block_size * cmd.max_rlines;
        for wf_adc = param.analysis.cmd{cmd_idx}.wf_adcs{img}(:).'
          wf = param.analysis.imgs{img}(wf_adc,1);
          adc = param.analysis.imgs{img}(wf_adc,2);
          out_fn = fullfile(out_fn_dir_dir,sprintf('specular_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc));
          sparam.success = cat(2,sparam.success, ...
            sprintf('  error_mask = bitor(error_mask,%d*~ct_file_lock_check(''%s'',4));\n', success_error, out_fn));
          if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
            ct_file_lock_check(out_fn,3);
          end
          sparam.cpu_time = sparam.cpu_time + Nx_cmd*total_num_sam(img)*log2(Nx_cmd)*cpu_time_mult;
          sparam.mem = max(sparam.mem,250e6 + Nx_cmd*total_num_sam(img)*mem_mult);
        end
        
      case {'statistics'}
        Nx_cmd = Nx / cmd.block_ave;
        for wf_adc = param.analysis.cmd{cmd_idx}.wf_adcs{img}(:).'
          wf = param.analysis.imgs{img}(wf_adc,1);
          adc = param.analysis.imgs{img}(wf_adc,2);
          out_fn = fullfile(out_fn_dir_dir,sprintf('stats_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc));
          sparam.success = cat(2,sparam.success, ...
            sprintf('  error_mask = bitor(error_mask,%d*~ct_file_lock_check(''%s'',4));\n', success_error, out_fn));
          if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
            ct_file_lock_check(out_fn,3);
          end
          if cmd.block_ave < 64
            % HACK: Assume that if no block averaging is done, that the data size
            % is small. Really need to get a user "hint" here.
            sparam.cpu_time = sparam.cpu_time + Nx_cmd*total_num_sam(img)/128*log2(Nx_cmd)*cpu_time_mult;
            sparam.mem = max(sparam.mem,250e6 + Nx_cmd*total_num_sam(img)*mem_mult);
          else
            sparam.cpu_time = sparam.cpu_time + Nx_cmd*total_num_sam(img)*log2(Nx_cmd)*cpu_time_mult;
            sparam.mem = max(sparam.mem,250e6 + Nx_cmd*total_num_sam(img)*mem_mult);
          end
        end
        
      case {'waveform'}
        Nx_cmd = Nx / param.analysis.dec;
        if isfinite(param.analysis.surf.Nt)
          Nt = param.analysis.surf.Nt;
        end
        for wf_adc = param.analysis.cmd{cmd_idx}.wf_adcs{img}(:).'
          wf = param.analysis.imgs{img}(wf_adc,1);
          adc = param.analysis.imgs{img}(wf_adc,2);
          out_fn = fullfile(out_fn_dir_dir,sprintf('waveform_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc));
          sparam.success = cat(2,sparam.success, ...
            sprintf('  error_mask = bitor(error_mask,%d*~ct_file_lock_check(''%s'',4));\n', success_error, out_fn));
          if ~ctrl.cluster.rerun_only && exist(out_fn,'file')
            ct_file_lock_check(out_fn,3);
          end
          sparam.cpu_time = sparam.cpu_time + Nx_cmd*Nt*log2(Nx_cmd)*cpu_time_mult;
          sparam.mem = max(sparam.mem,250e6 + Nx_cmd*Nt*mem_mult);
        end
        
    end
  end
end
sparam.notes = sprintf('%s:%s:%s %s combine', ...
  mfilename, param.radar_name, param.season_name, param.day_seg);

ctrl = cluster_new_task(ctrl,sparam,[]);

ctrl_chain{end+1} = ctrl;
    
fprintf('Done %s\n', datestr(now));

return
