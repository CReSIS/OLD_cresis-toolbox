function combine_wf_chan_ollie(steady_param_file_name,frm,chunk)
% combine_wf_chan_ollie(steady_param_file_name,frm,chunk)
%
% This script combines the receive channels and outputs the result
% for each waveform. It also combines the waveforms. It takes in
% f-k files one directory at a time and:
%  1. combines the receive channels
%  2. concatenates the results
%  3. square-law detects the data, abs()^2
%  4. takes incoherent averages (multilooks data)
%  5. saves the result in a new directory
%
% The assumption is that the directories in the input_path are named
% using the following convention:
%   PROC-TYPE-STRING_data_#{_SUBAPERTURE-STRING}
% where
%   PROC-TYPE-STRING can be 'fk','tdbp', or 'pc' for f-k migrated,time domain
%   back projected,and pulse compressed respectively ('fk' and tdbp supported)
%   _data_ is always present
%   #, \d+: one or more numbers
%   _SUBAPERTURE-STRING, {_[mp]\d\.\d}: optional subaperture string
% Examples:
%   fk_data_01_01: f-k migrated, frame 1, subaperture 1
%   fk_data_04_02: f-k migrated, frame 4, subaperture 2
%   fk_data_01_03: f-k migrated, frame 1, subaperture 3
%   pc_data_01: pulse compressed only, frame 1
%
% param = struct with processing parameters loaded from file
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_combine_wf_chan.m, combine_wf_chan.m,
%   combine_wf_chan_task.m

% =====================================================================
% General Setup
% =====================================================================

dbstack_info = dbstack;

% Test
if (ischar(frm))
  frm=str2num(frm);
end
if (ischar(chunk))
  chunk=str2num(chunk);
end
  
load(steady_param_file_name,'steady_param');
param=steady_param;

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

% =====================================================================
% Setup processing
% =====================================================================

% Get WGS84 ellipsoid parameters
physical_constants;

if ~isfield(param,'debug_level')
  param.debug_level = 1;
end

if ~isfield(param.sched,'rerun_only') || isempty(param.sched.rerun_only)
  param.sched.rerun_only = false;
end

% Handles multilooking syntax:
%  {{[1 1],[1 2],[1 3],[1 4],[1 5]},{[2 1],[2 2],[2 3],[2 4],[2 5]}}
%  If the image is a cell array it describes multilooking across apertures
if ~iscell(param.combine.imgs{1})
  % No special multilooking, reformat old syntax to new multilooking syntax
  for img = 1:length(param.combine.imgs)
    param.combine.imgs{img} = {param.combine.imgs{img}};
  end
end

for img = 1:length(param.combine.imgs)
  for ml_idx = 1:length(param.combine.imgs{img})
    % Imaginary image indices is for IQ combining during raw data load
    % which we do not need here.
    param.combine.imgs{img}{ml_idx} = abs(param.combine.imgs{img}{ml_idx});
  end
end

img_list = param.combine.imgs;
  
in_path = ct_filename_out(param, ...
  param.combine.in_path, 'CSARP_out');

array_path = ct_filename_out(param, ...
  param.combine.array_path, 'CSARP_out');

out_path = ct_filename_out(param, ...
  param.combine.out_path, sprintf('CSARP_%s', ...
  param.combine.method));

% Create the output directory
if ~exist(out_path,'dir')
  mkdir(out_path);
end

param.surf.manual = 0; % Turn manual pick off

%% Loop through all the frame directories and process the fk
% chunks in those directories
% =====================================================================
retry_fields = {};

%% Input directory for this frame (only look at the first subaperture
% "_01_01" since combine_wf_chan_task will know which subapertures to load)
if strcmpi(param.csarp.sar_type,'f-k')
  sar_type = 'fk';
elseif strcmpi(param.csarp.sar_type,'tdbp')
  sar_type = 'tdbp';
elseif strcmpi(param.csarp.sar_type,'mltdp')
  sar_type = 'mltdp';
end
param.combine.in_path = fullfile(in_path,sprintf('%s_data_%03d_01_01', sar_type, frm));

%% Output directory
param.combine.out_path = fullfile(array_path,sprintf('array_%03d', frm));

%% Get the filenames for each chunk of data processed by csarp
% DEBUG:
%  - Get all the chunk files
%  - Normal operation is 1 to +inf (i.e. all files)
%  - The start/stop funtionality is not used except for debugging
%  - Set start_chunk and stop_chunk to restrict which chunks get processed
start_chunk = chunk;
stop_chunk = chunk;
% Get all time stamps in the directory: the assumption is that csarp.m
% has created all the necessary files for each wf/adc pair required.  So
% to get the time stamps, we just search for all the files for a particular
% pair (in this case the first one in the list).
img = 1;
wf_adc_idx = 1;
filenames = get_filenames(param.combine.in_path,'', ...
  sprintf('wf_%02d_adc_%02d',img_list{img}{1}(wf_adc_idx,1), ...
  img_list{img}{1}(wf_adc_idx,2)),'.mat');
if isempty(filenames)
  error('No filenames found in %s', param.combine.in_path);
end
chunk_ids = {};
for idx = 1:length(filenames)
  [path,name,ext] = fileparts(filenames{idx});
  chunk_idx = str2double(name(end-2:end));
  if chunk_idx >= start_chunk && chunk_idx <=stop_chunk
    chunk_ids{end+1} = name(end-2:end);
  end
end

%% Create and clean the array_proc output directories
if exist(param.combine.out_path,'dir') && ~param.sched.rerun_only
  % If folders do exist, clear them out
    fprintf('  Directory %s already exists, files will be overwritten \n', param.combine.out_path);
    %rmdir(param.combine.out_path,'s');
end
mkdir(param.combine.out_path);

%% Combine Channels: Standard beam-forming, MUSIC, MVDR, etc)
% - This is setup so that multiple fk chunks can be processed in the
%   same task/job.

load(filenames{1},'param_csarp');
if strcmpi(param_csarp.csarp.sar_type,'f-k')
  if max(param.combine.rline_rng) - min(param.combine.rline_rng) > param_csarp.csarp.chunk_overlap
    error('SAR processing chunks will not align properly, chunk_overlap too small');
  end
end
param.combine.sar_type = param_csarp.csarp.sar_type;
num_chunks_per_task = 1;
  
for chunk_idx = 1:num_chunks_per_task:length(chunk_ids)
  %% To make the SAR processed chunks fit together seamlessly without
  % having to resample, we determine the start range line output for
  % each chunk.
  % chunk_idxs: SAR chunks that will be processed by this task,
  %   plus one additional one for calculating rlines(2)
  chunk_idxs = chunk_idx + (0:num_chunks_per_task);
  % chunk_Nx: the number of non-overlapping SAR chunk outputs
  chunk_Nx = floor(param_csarp.csarp.chunk_len / param_csarp.csarp.sigma_x);
  % min_offset: the minimum offset into the SAR chunk which array_proc can
  %   output a full support estimate (since the output uses a neighborhood
  %   of points around the pixel in question, the first output line generally
  %   be from the first input line)
  min_offset = -min(param.combine.rline_rng);
  % rlines(1,:): this will be the first range line output by array_proc
  %   for each SAR chunk this task is array processing
  rlines = [];
  rlines(1,:) = 1+mod(min_offset+param.combine.dline-(chunk_idxs-1)*chunk_Nx, param.combine.dline);
  rlines(rlines<min_offset) = rlines(rlines<min_offset) + ceil(param.combine.dline/(1+min_offset)) * param.combine.dline;
  rlines(2,1:end-1) = chunk_Nx + rlines(1,2:end) - param.combine.dline;
  rlines = rlines(:,1:end-1);
  
  % Check if this is the last chunk. This last chunk could have variable
  % length and we want to return all of the data from this chunk. To tell
  % combine_task to do this, we set rlines(2) to infinity for
  % this chunk
  if chunk_idx+num_chunks_per_task-1 >= length(chunk_ids)
    rlines(2,end) = inf;
  end
  
  %% Get the chunk ids that this task will process
  chunk_idx_last = min(chunk_idx+num_chunks_per_task-1, length(chunk_ids));
  param.combine.chunk_ids = chunk_ids(chunk_idx:chunk_idx_last);
  param.combine.rlines = rlines;
  
  %% Pass in the frame
  param.combine.frm = frm;
  
  %% Rerun only mode: Test to see if we need to run this task
  if param.sched.rerun_only
    % If we are in rerun only mode AND all the combine_wf_chan task output files
    % already exists, then we do not run the task
    file_exists = true;
    for img = 1:length(param.combine.imgs)
%         array_fn = fullfile(param.combine.out_path, ...
%           sprintf('chk_%03d_img_%02d.mat', chunk_ids{chunk_idx}, img));
      array_fn = fullfile(param.combine.out_path, ...
        sprintf('chk_%s_img_%02d.mat', chunk_ids{chunk_idx}, img));
      if ~exist(array_fn,'file')
        file_exists = false;
      end
    end
    if file_exists
      fprintf('  %s already exists [rerun_only skipping] (%s)\n', ...
        param.combine.chunk_ids{1}, datestr(now));
      continue;
    end
  end

  % param.combine = structure controlling the processing
  % Fields used in this function (and potentially array_proc.m)
  % .imgs = images to process
  % .method = 'csarp-combined', 'standard', 'mvdr', 'music', etc.
  % .debug_level = scalar integer specifying debug level
  % .three_dim.en = boolean, enable 3-D mode
  % .Nsv = used to compute sv
  % .sv = steering vector function handle (empty for default)
  % .in_path = input path string
  % .out_path = output path string
  % .chunk_ids = specific csarp output chunk IDs
  %
  % Fields used by array_proc.m only
  % .window
  % .bin_rng
  % .rline_rng
  % .dbin
  % .dline
  % .first_line = required to stitch multiple chunks together
  % .chan_equal
  % .freq_rng
  
  % =====================================================================
  % Process data
  % =====================================================================
  for chunk_idx = 1:length(param.combine.chunk_ids)
    
    tx = 1;
    for img_idx = 1:length(param.combine.imgs)
      ml_list = param.combine.imgs{img_idx};
      
      %% Convert old syntax to new multilooking across SLC images syntax
      % SLC = single look complex SAR echogram
      if ~iscell(ml_list)
        ml_list = {ml_list};
      end
      wf_base = ml_list{1}(1,1);
      
      % -------------------------------------------------------------------
      %% Load data
      %  wf_adc_list = {[1 1],[1 2],[1 3],[1 4],[1 5]}
      %    This will multilook these 5 SLC images (e.g. add incoherently in the case
      %    of periodogram)
      %  wf_adc_list = {[1 1; 1 2; 1 3; 1 4; 1 5]}
      %    This will coherently add these 5 SLC images
      if strcmpi(param.combine.method,'csarp-combined')
        % CSARP-combined has already coherently combined every image in each
        % image list, so we only grab the first image.  Multilooking across
        % SLCs IS NOT supported for this method!
        ml_list = {ml_list{1}(1,:)};
      end
      if ~isfield(param.combine,'subaps') || isempty(param.combine.subaps)
        % If SAR sub-apertures not set, we assume that there is just one
        % subaperture to be passed in for each multilook input
        for ml_idx = 1:length(ml_list)
          param.combine.subaps{ml_idx} = [1];
        end
      end
      if ~isfield(param.combine,'subbnds') || isempty(param.combine.subbnds)
        % If subbands not set, we assume that there is just one
        % subaperture to be passed in for each multilook input
        for ml_idx = 1:length(ml_list)
          param.combine.subbnds{ml_idx} = [1];
        end
      end
      clear fcs chan_equal data;
      for ml_idx = 1:length(ml_list)
        wf_adc_list = ml_list{ml_idx};
        for wf_adc_idx = 1:size(wf_adc_list,1)
          wf = wf_adc_list(wf_adc_idx,1);
          adc = wf_adc_list(wf_adc_idx,2);
          for subap = param.combine.subaps{ml_idx}
            for subbnd = param.combine.subbnds{ml_idx}
              in_path = param.combine.in_path;
              in_path(end-4:end-3) = sprintf('%02d',subap);
              in_path(end-1:end) = sprintf('%02d',subbnd);
              [sar_type_fn,status] = get_filenames(in_path, ...
                sprintf('wf_%02.0f_adc_%02.0f_', wf, adc),['chk_',param.combine.chunk_ids{chunk_idx}],'.mat');
              if param.debug_level >= 1
                fprintf('  Loading SAR processed data %s (%s)\n', sar_type_fn{1}, datestr(now));
              end
              sar_type_file = load(sar_type_fn{1});
              fcs{ml_idx}{wf_adc_idx} = sar_type_file.fcs;
              
              chan_equal{ml_idx}(wf_adc_idx) = 10.^(param.radar.wfs(wf).chan_equal_dB(param.radar.wfs(wf).rx_paths(adc))/20) ...
                .* exp(j*param.radar.wfs(wf).chan_equal_deg(param.radar.wfs(wf).rx_paths(adc))/180*pi);
              
              if strcmpi(param.combine.sar_type,'f-k')
                data{ml_idx}(:,:,subap,subbnd,wf_adc_idx) = sar_type_file.fk_data;
              elseif strcmpi(param.combine.sar_type,'tdbp')
                data{ml_idx}(:,:,subap,subbnd,wf_adc_idx) = sar_type_file.tdbp_data;
              elseif strcmpi(param.combine.sar_type,'mltdp')
                data{ml_idx}(:,:,subap,subbnd,wf_adc_idx) = sar_type_file.mltdp_data;
              end
            end
          end
        end
      end
      
      if isfield(param.combine,'ft_over_sample') && ~isempty(param.combine.ft_over_sample)
        % param.combine.ft_over_sample should be a positive integer
        for ml_idx = 1:length(data)
          data{ml_idx} = interpft(data{ml_idx},size(data{ml_idx},1) * param.combine.ft_over_sample);
        end
        for wf = 1:length(sar_type_file.wfs)
          sar_type_file.wfs(wf).fs = sar_type_file.wfs(wf).fs * param.combine.ft_over_sample;
          sar_type_file.wfs(wf).dt = 1/sar_type_file.wfs(wf).fs;
          sar_type_file.wfs(wf).Nt = sar_type_file.wfs(wf).Nt * param.combine.ft_over_sample;
          sar_type_file.wfs(wf).df = sar_type_file.wfs(wf).fs / sar_type_file.wfs(wf).Nt;
          sar_type_file.wfs(wf).time = sar_type_file.wfs(wf).time(1) + sar_type_file.wfs(wf).dt*(0:sar_type_file.wfs(wf).Nt-1).';
          sar_type_file.wfs(wf).freq = sar_type_file.wfs(wf).fc ...
            + sar_type_file.wfs(wf).df * ifftshift( -floor(sar_type_file.wfs(wf).Nt/2) : floor((sar_type_file.wfs(wf).Nt-1)/2) ).';
        end
      end
      
      % Setup fields for wideband space-time doa estimator passed to
      % array_proc.  This section does the following:
      % 1). Computes the impulse response used to model the amplitude term of
      % the space time correlation matrix,
      if strcmpi(param.combine.method,'wideband')
        % -------------------------------------------------------------------
        % Compute maximum propagation delay across the array (this is used to
        % compute impulse response needed for covariance model (see below))
        
        lever_arm_param.season_name = param.season_name;
        lever_arm_param.radar_name = ct_output_dir(param.radar_name);
        lever_arm_param.gps_source = sar_type_file.param_records.gps_source;
        % Iterate through all wf-adc pairs and create a list of phase centers
        % that are being used to create this image.
        phase_centers = [];
        for ml_idx = 1:length(ml_list)
          wf_adc_list = ml_list{ml_idx};
          for wf_adc_idx = 1:size(wf_adc_list,1)
            wf = wf_adc_list(wf_adc_idx,1);
            adc = wf_adc_list(wf_adc_idx,2);
            % Add phase center for the wf-adc pair to the list
            phase_centers(:,end+1) = lever_arm(lever_arm_param, ...
              param.radar.wfs(wf).tx_weights, param.radar.wfs(wf).rx_paths(adc));
          end
        end
        % Find the two phase centers separated by the maximum distance
        pc_dist = zeros(size(phase_centers,2));
        for pc_idx = 1:size(phase_centers,2)
          for pc_comp_idx = pc_idx+1:size(phase_centers,2)
            pc_dist(pc_idx,pc_comp_idx) ...
              = sqrt(sum(abs(phase_centers(:,pc_idx) - phase_centers(:,pc_comp_idx)).^2));
          end
        end
        max_array_dim = max(max(pc_dist));
        % Convert maximum distance to propagation time in free space
        tau_max = 2*max_array_dim/c;
        
        % Check the value for W from param.combine.W (widening factor)
        W_ideal = 1 + ceil(tau_max / sar_type_file.wfs(wf).dt);
        if param.combine.W ~= W_ideal
          warning('param.combine.W is %d, but should be %d', param.combine.W, W_ideal);
        end
        
        % Compute impulse response
        % -------------------------------------------------------------------
        % NOTE:  The impulse response is computed based on the fast time
        % window used in CSARP.  This could be generalized for a measured
        % impulse response vector but this is not yet supported.
        
        Mt = 10; % Over-sampling factor
        % Create frequency-domain fast-time window identical to what was
        % used for pulse compression
        Hwin = sar_type_file.param_csarp.csarp.ft_wind(length(sar_type_file.wfs(wf).time));
        % Convert to time-domain and over-sample by Mt
        %  - Take real to remove rounding errors that result in imag part
        Hwin = interpft(real(ifft(ifftshift(Hwin).^2)), Mt*length(Hwin));
        % Store the impulse response and corresponding time axis for
        % passing to array_proc
        %  - Ensure we grab enough samples of the impulse response so that
        %    array_proc is always happy.
        Hwin_num_samp = 2 * Mt * (W_ideal + param.combine.W);
        param.combine.imp_resp.vals ...
          = fftshift(Hwin([1:1+Hwin_num_samp, end-Hwin_num_samp+1:end]));
        param.combine.imp_resp.time_vec ...
          = sar_type_file.wfs(wf).dt/Mt * (-Hwin_num_samp:Hwin_num_samp);
      end
      
      if strcmpi(param.combine.method,'csarp-combined')
        %% csarp.m already combined channels, just incoherently average
        
        % Number of fast-time samples in the data
        Nt = size(data{1},1);
        
        param.combine.bins = numel(param.combine.bin_rng)/2+0.5 : param.combine.dbin ...
          : Nt-(numel(param.combine.bin_rng)/2-0.5);
        
        param.combine.lines = param.combine.rlines(1,chunk_idx): param.combine.dline ...
          : min(param.combine.rlines(2,chunk_idx),size(data{1},2)-max(param.combine.rline_rng));
        
        % Perform incoherent averaging
        Hfilter2 = ones(length(param.combine.bin_rng),length(param.combine.rline_rng));
        Hfilter2 = Hfilter2 / numel(Hfilter2);
        tomo.val = filter2(Hfilter2, abs(data{1}).^2);
        
        tomo.val = tomo.val(param.combine.bins,param.combine.lines);
        
        param.combine.chan_equal = chan_equal;
        param.combine.fcs = fcs;
        array_param = param.combine;
        
      else
        %% Regular array processing operation
        
        % Load bin restriction layers if specified
        if isfield(param.combine,'bin_restriction') && ~isempty(param.combine.bin_restriction)
          ops_param = param;
          % Invalid frames will be ignored by opsLoadLayers so we can ignore edge cases
          ops_param.cmd.frms = param.combine.frm + (-1:1);
          param.combine.bin_restriction = opsLoadLayers(ops_param, param.combine.bin_restriction);
          for idx=1:2
            % Interpolate layers onto GPS time of loaded data
            param.combine.bin_restriction(idx).bin = interp1(param.combine.bin_restriction(idx).gps_time, ...
              param.combine.bin_restriction(idx).twtt, sar_type_file.fcs.gps_time);
            % Convert from two way travel time to range bins
            param.combine.bin_restriction(idx).bin = interp1(sar_type_file.wfs(wf).time, ...
              1:length(sar_type_file.wfs(wf).time), param.combine.bin_restriction(idx).bin);
            % Ensure there is a value everywhere
            param.combine.bin_restriction(idx).bin = interp_finite(param.combine.bin_restriction(idx).bin);
          end
        end
        
        if isfield(param.combine,'doa_constraints') && ~isempty(param.combine.doa_constraints)
          % If DOA constraints are specified, they must be specified for
          % every signal. For no constraints, choose 'fixed' method and
          % [-pi/2 pi/2] for init_src_limits and src_limits.
          for src_idx = 1:param.combine.Nsig
            % Load layers for each DOA constraint that needs it
            doa_res = param.combine.doa_constraints(src_idx);
            switch (doa_res.method)
              case {'layerleft','layerright'}
                ops_param = param;
                % Invalid frames will be ignored by opsLoadLayers so we can ignore edge cases
                ops_param.cmd.frms = param.combine.frm + (-1:1);
                param.combine.doa_constraints(src_idx).layer = opsLoadLayers(ops_param, doa_res.params);
                % Interpolate layers onto GPS time of loaded data
                param.combine.doa_constraints(src_idx).layer.twtt = interp1(param.combine.doa_constraints(src_idx).layer.gps_time, ...
                  param.combine.doa_constraints(src_idx).layer.twtt, sar_type_file.fcs.gps_time);
                % Ensure there is a value everywhere
                param.combine.doa_constraints(src_idx).layer.twtt = interp_finite(param.combine.doa_constraints(src_idx).layer.twtt);
            end
          end
        end
        
        param.combine.chan_equal = chan_equal;
        param.combine.fcs = fcs;
        array_param = param.combine;
        array_param.wfs = sar_type_file.wfs(wf);
        array_param.imgs = param.combine.imgs{img_idx};
        if ~iscell(array_param.imgs)
          array_param.imgs = {array_param.imgs};
        end
        array_param.rlines = param.combine.rlines(:,chunk_idx);
        
        % Array Processing Function Call
        [array_param,tomo] = array_proc(array_param,data);
        param.combine.bins = array_param.bins;
        param.combine.lines = array_param.lines;
      end
      
      % -------------------------------------------------------------------
      % Debugging
      if param.debug_level >= 2
        figure(1); clf;
        %imagesc(incfilt(mean(data{1},3),10,4));
        imagesc(10*log10(local_detrend(abs(mean(data{1},3)).^2,[40 100],[10 4],3)));
        title(sprintf('Waveform %d: Boxcar window', wf));
        colorbar
        grid on;
        figure(2); clf;
        %imagesc(10*log10(tomo.val));
        imagesc(10*log10(local_detrend(tomo.val,[20 50],[5 2],3)));
        title(sprintf('Waveform %d: Tomography', wf));
        colorbar
        grid on;
        keyboard
      end
      
      % -------------------------------------------------------------------
      %% Save results
      [sar_type_fn_path sar_type_fn_name sar_type_fn_ext] = fileparts(sar_type_fn{1});
      array_fn = fullfile(param.combine.out_path, [sar_type_fn_name(end-6:end) sprintf('_img_%02d', img_idx) sar_type_fn_ext]);
      fprintf('  Saving combine_wf_chan data %s (%s)\n', array_fn, datestr(now));
      
      if ~isfield(sar_type_file,'lat')
        warning('DEPRECATED CODE: You are using an old SAR file version without lat field. Computing position from the origin of the SAR coordinate system. Support for old SAR files will eventually be removed.');
        [Latitude,Longitude,Elevation] ...
          = ecef2geodetic(sar_type_file.fcs.origin(1,array_param.lines), ...
          sar_type_file.fcs.origin(2,array_param.lines), ...
          sar_type_file.fcs.origin(3,array_param.lines),WGS84.ellipsoid);
        Latitude = Latitude*180/pi;
        Longitude = Longitude*180/pi;
      else
        Latitude = sar_type_file.lat(1,array_param.lines);
        Longitude = sar_type_file.lon(1,array_param.lines);
        Elevation = sar_type_file.elev(1,array_param.lines);
      end
      GPS_time = sar_type_file.fcs.gps_time(array_param.lines);
      Surface = sar_type_file.fcs.surface(array_param.lines);
      Bottom = sar_type_file.fcs.bottom(array_param.lines);
      Roll = sar_type_file.fcs.roll(array_param.lines);
      Pitch = sar_type_file.fcs.pitch(array_param.lines);
      Heading = sar_type_file.fcs.heading(array_param.lines);
      Data = tomo.val;
      Time = sar_type_file.wfs(wf_base).time(array_param.bins);
      param_records = sar_type_file.param_records;
      param_csarp = sar_type_file.param_csarp;
      array_param.sv = []; % Set this to empty because it is so large
      param_combine = param;
      param_combine.array_param = array_param;
      if ~param.combine.three_dim.en
        % Do not save 3D-image
        save('-v7.3',array_fn,'Data','Latitude','Longitude','Elevation','GPS_time', ...
          'Surface','Bottom','Time','param_combine','param_records', ...
          'param_csarp', 'Roll', 'Pitch', 'Heading');
      else
        % Save 3D-image in file
        Topography = tomo;
        save('-v7.3',array_fn,'Topography','Data','Latitude','Longitude','Elevation','GPS_time', ...
          'Surface','Bottom','Time','param_combine','param_records', ...
          'param_csarp', 'Roll', 'Pitch', 'Heading');
      end
    end
  end
  
end

return;