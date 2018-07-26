function success = array_task(param)
% success = array_task(param)
%
% Task/job running on cluster that is called from combine. This
% function is generally not called directly.
%
% param.combine = structure controlling the processing
%   Fields used in this function (and potentially array_proc.m)
%   .imgs = images to process
%   .method = 'csarp-combined', 'standard', 'mvdr', 'music', etc.
%   .three_dim.en = boolean, enable 3-D mode
%   .Nsv = used to compute sv
%   .sv = steering vector function handle (empty for default)
%   .in_path = input path string
%   .out_path = output path string
%
%   Fields used by array_proc.m only
%   .window
%   .bin_rng
%   .rline_rng
%   .dbin
%   .dline
%   .first_line = required to stitch multiple chunks together
%   .chan_equal
%   .freq_rng
%
% success = boolean which is true when the function executes properly
%   if a task fails before it can return success, then success will be
%   empty
%
% Authors: John Paden
%
% See also: combine, array_proc

% =====================================================================
% Preparation

% speed of light, wgs84 ellipsoid
physical_constants;

% Rename for code readability
sar_type = param.combine.sar_type;
data_field_name = sprintf('%s_data', sar_type);

% Input directory (SAR SLC images)
in_fn_dir = fullfile(ct_filename_out(param, param.combine.in_path), ...
  sprintf('%s_data_%03d_01_01', sar_type, param.load.frm));

% Temporary output directory for uncombined array processed images
combine_tmp_dir = fullfile(ct_filename_out(param, param.combine.array_path), ...
  sprintf('%s_%03d', param.combine.method, param.load.frm));
  
% =====================================================================
% Process data
% =====================================================================
for img = 1:length(param.combine.imgs)
  ml_list = param.combine.imgs{img};
  
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
  clear fcs chan_equal data lat lon elev;
  num_next_rlines = 0;
  prev_chunk_failed_flag = false;
  next_chunk_failed_flag = false;
  % num_prev_chunk_rlines: number of range lines to load from the previous chunk
  num_prev_chunk_rlines = max(-param.combine.rline_rng);
  % num_next_chunk_rlines: number of range lines to load from the previous chunk
  num_next_chunk_rlines = max(param.combine.rline_rng);
  % out_rlines: output range lines from SAR processor for the current chunk
  sar_out_rlines = [];
  for ml_idx = 1:length(ml_list)
    wf_adc_list = ml_list{ml_idx};
    for wf_adc_idx = 1:size(wf_adc_list,1)
      wf = wf_adc_list(wf_adc_idx,1);
      adc = wf_adc_list(wf_adc_idx,2);
      for subap = param.combine.subaps{ml_idx}
        for subbnd = param.combine.subbnds{ml_idx}
          in_fn_dir(end-4:end-3) = sprintf('%02d',subap);
          in_fn_dir(end-1:end) = sprintf('%02d',subbnd);
          
          % Load previous chunk data
          % ===============================================================
          
          % If on the first chunk of the frame, then look at the previous
          % frame. Special cases (like frm == 1 and empty
          % param.load.prev_frm_num_chunks are handled by letting the file
          % search fail).
          if param.load.chunk_idx == 1
            load_frm = param.load.frm-1;
            load_chunk_idx = param.load.prev_frm_num_chunks;
          else
            load_frm = param.load.frm;
            load_chunk_idx = param.load.chunk_idx - 1;
          end

          % Create the filename
          in_fn_dir(end-8:end-6) = sprintf('%03d',load_frm);
          if strcmpi(param.combine.method,'csarp-combined')
            [sar_type_fn,status] = get_filenames(in_fn_dir, ...
              sprintf('img_%02.0f_chk_%03.0f', img, load_chunk_idx),'','.mat');
          else
            [sar_type_fn,status] = get_filenames(in_fn_dir, ...
              sprintf('wf_%02.0f_adc_%02.0f_chk_%03.0f', wf, adc, load_chunk_idx),'','.mat');
          end
          
          if ~prev_chunk_failed_flag && ~isempty(sar_type_fn)
            % If file exists, then load it
            sar_data = load(sar_type_fn{1});
            rlines = 1:num_prev_chunk_rlines;
            
            if size(sar_data.(data_field_name),2) >= num_prev_chunk_rlines
              if subap == param.combine.subaps{ml_idx}(1) && subbnd == param.combine.subbnds{ml_idx}(1)
                fcs{ml_idx}{wf_adc_idx}.origin(:,rlines) = sar_data.fcs.origin(:,end-num_prev_chunk_rlines+1:end);
                fcs{ml_idx}{wf_adc_idx}.x(:,rlines) = sar_data.fcs.x(:,end-num_prev_chunk_rlines+1:end);
                fcs{ml_idx}{wf_adc_idx}.y(:,rlines) = sar_data.fcs.y(:,end-num_prev_chunk_rlines+1:end);
                fcs{ml_idx}{wf_adc_idx}.z(:,rlines) = sar_data.fcs.z(:,end-num_prev_chunk_rlines+1:end);
                fcs{ml_idx}{wf_adc_idx}.roll(rlines) = sar_data.fcs.roll(end-num_prev_chunk_rlines+1:end);
                fcs{ml_idx}{wf_adc_idx}.pitch(rlines) = sar_data.fcs.pitch(end-num_prev_chunk_rlines+1:end);
                fcs{ml_idx}{wf_adc_idx}.heading(rlines) = sar_data.fcs.heading(end-num_prev_chunk_rlines+1:end);
                fcs{ml_idx}{wf_adc_idx}.gps_time(rlines) = sar_data.fcs.gps_time(end-num_prev_chunk_rlines+1:end);
                fcs{ml_idx}{wf_adc_idx}.surface(rlines) = sar_data.fcs.surface(end-num_prev_chunk_rlines+1:end);
                fcs{ml_idx}{wf_adc_idx}.bottom(rlines) = sar_data.fcs.bottom(end-num_prev_chunk_rlines+1:end);
                fcs{ml_idx}{wf_adc_idx}.pos(:,rlines) = sar_data.fcs.pos(:,end-num_prev_chunk_rlines+1:end);
              end
              
              data{ml_idx}(:,rlines,subap,subbnd,wf_adc_idx) = sar_data.(data_field_name)(:,end-num_prev_chunk_rlines+1:end);
              lat(rlines) = sar_data.lat(:,end-num_prev_chunk_rlines+1:end);
              lon(rlines) = sar_data.lon(:,end-num_prev_chunk_rlines+1:end);
              elev(rlines) = sar_data.elev(:,end-num_prev_chunk_rlines+1:end);
            else
              prev_chunk_failed_flag = true;
            end
          else
            % If file does not exist, then give up loading previous chunk
            % data
            prev_chunk_failed_flag = true;
          end
          
          % Load current chunk data
          % ===============================================================
          load_frm = param.load.frm;
          load_chunk_idx = param.load.chunk_idx;
          in_fn_dir(end-8:end-6) = sprintf('%03d',load_frm);
          if strcmpi(param.combine.method,'csarp-combined')
            [sar_type_fn,status] = get_filenames(in_fn_dir, ...
              sprintf('img_%02.0f_chk_%03.0f', img, load_chunk_idx),'','.mat');
          else
            [sar_type_fn,status] = get_filenames(in_fn_dir, ...
              sprintf('wf_%02.0f_adc_%02.0f_chk_%03.0f', wf, adc, load_chunk_idx),'','.mat');
          end
          sar_data = load(sar_type_fn{1});
          num_rlines = size(sar_data.(data_field_name),2);
          rlines = num_prev_chunk_rlines + (1:num_rlines);
          
          if subap == param.combine.subaps{ml_idx}(1) && subbnd == param.combine.subbnds{ml_idx}(1)
            fcs{ml_idx}{wf_adc_idx}.Lsar = sar_data.fcs.Lsar;
            fcs{ml_idx}{wf_adc_idx}.gps_source = sar_data.fcs.gps_source;
            fcs{ml_idx}{wf_adc_idx}.origin(:,rlines) = sar_data.fcs.origin;
            fcs{ml_idx}{wf_adc_idx}.x(:,rlines) = sar_data.fcs.x;
            fcs{ml_idx}{wf_adc_idx}.y(:,rlines) = sar_data.fcs.y;
            fcs{ml_idx}{wf_adc_idx}.z(:,rlines) = sar_data.fcs.z;
            fcs{ml_idx}{wf_adc_idx}.roll(rlines) = sar_data.fcs.roll;
            fcs{ml_idx}{wf_adc_idx}.pitch(rlines) = sar_data.fcs.pitch;
            fcs{ml_idx}{wf_adc_idx}.heading(rlines) = sar_data.fcs.heading;
            fcs{ml_idx}{wf_adc_idx}.gps_time(rlines) = sar_data.fcs.gps_time;
            fcs{ml_idx}{wf_adc_idx}.surface(rlines) = sar_data.fcs.surface;
            fcs{ml_idx}{wf_adc_idx}.bottom(rlines) = sar_data.fcs.bottom;
            fcs{ml_idx}{wf_adc_idx}.pos(:,rlines) = sar_data.fcs.pos;
            fcs{ml_idx}{wf_adc_idx}.squint = sar_data.fcs.squint;
            if isfield(sar_data.fcs,'type')
              fcs{ml_idx}{wf_adc_idx}.type = sar_data.fcs.type;
              fcs{ml_idx}{wf_adc_idx}.filter = sar_data.fcs.filter;
            else
              fcs{ml_idx}{wf_adc_idx}.type = [];
              fcs{ml_idx}{wf_adc_idx}.filter = [];
            end
            sar_out_rlines = sar_data.out_rlines;
          
            chan_equal{ml_idx}(wf_adc_idx) = 10.^(param.radar.wfs(wf).chan_equal_dB(param.radar.wfs(wf).rx_paths(adc))/20) ...
              .* exp(j*param.radar.wfs(wf).chan_equal_deg(param.radar.wfs(wf).rx_paths(adc))/180*pi);
          end
          
          data{ml_idx}(:,rlines,subap,subbnd,wf_adc_idx) = sar_data.(data_field_name);
          lat(rlines) = sar_data.lat;
          lon(rlines) = sar_data.lon;
          elev(rlines) = sar_data.elev;
          
          % Load next chunk data
          % ===============================================================
          
          % If on the first chunk of the frame, then look at the previous
          % frame. Special cases (like frm == 1 and empty
          % param.load.prev_frm_num_chunks are handled by letting the file
          % search fail).
          if param.load.chunk_idx == param.load.num_chunks
            load_frm = param.load.frm+1;
            load_chunk_idx = 1;
          else
            load_frm = param.load.frm;
            load_chunk_idx = param.load.chunk_idx + 1;
          end

          % Create the filename
          in_fn_dir(end-8:end-6) = sprintf('%03d',load_frm);
          if strcmpi(param.combine.method,'csarp-combined')
            [sar_type_fn,status] = get_filenames(in_fn_dir, ...
              sprintf('img_%02.0f_chk_%03.0f', img, load_chunk_idx),'','.mat');
          else
            [sar_type_fn,status] = get_filenames(in_fn_dir, ...
              sprintf('wf_%02.0f_adc_%02.0f_chk_%03.0f', wf, adc, load_chunk_idx),'','.mat');
          end
          
          if ~next_chunk_failed_flag && ~isempty(sar_type_fn)
            % If file exists, then load it
            sar_data = load(sar_type_fn{1});
            rlines = num_prev_chunk_rlines + num_rlines + (1:num_next_chunk_rlines);
            
            if size(sar_data.(data_field_name),2) >= num_next_chunk_rlines
              if subap == param.combine.subaps{ml_idx}(1) && subbnd == param.combine.subbnds{ml_idx}(1)
                fcs{ml_idx}{wf_adc_idx}.origin(:,rlines) = sar_data.fcs.origin(:,1:num_next_chunk_rlines);
                fcs{ml_idx}{wf_adc_idx}.x(:,rlines) = sar_data.fcs.x(:,1:num_next_chunk_rlines);
                fcs{ml_idx}{wf_adc_idx}.y(:,rlines) = sar_data.fcs.y(:,1:num_next_chunk_rlines);
                fcs{ml_idx}{wf_adc_idx}.z(:,rlines) = sar_data.fcs.z(:,1:num_next_chunk_rlines);
                fcs{ml_idx}{wf_adc_idx}.roll(rlines) = sar_data.fcs.roll(1:num_next_chunk_rlines);
                fcs{ml_idx}{wf_adc_idx}.pitch(rlines) = sar_data.fcs.pitch(1:num_next_chunk_rlines);
                fcs{ml_idx}{wf_adc_idx}.heading(rlines) = sar_data.fcs.heading(1:num_next_chunk_rlines);
                fcs{ml_idx}{wf_adc_idx}.gps_time(rlines) = sar_data.fcs.gps_time(1:num_next_chunk_rlines);
                fcs{ml_idx}{wf_adc_idx}.surface(rlines) = sar_data.fcs.surface(1:num_next_chunk_rlines);
                fcs{ml_idx}{wf_adc_idx}.bottom(rlines) = sar_data.fcs.bottom(1:num_next_chunk_rlines);
                fcs{ml_idx}{wf_adc_idx}.pos(:,rlines) = sar_data.fcs.pos(:,1:num_next_chunk_rlines);
              end
              
              data{ml_idx}(:,rlines,subap,subbnd,wf_adc_idx) = sar_data.(data_field_name)(:,1:num_next_chunk_rlines);
              lat(rlines) = sar_data.lat(:,1:num_next_chunk_rlines);
              lon(rlines) = sar_data.lon(:,1:num_next_chunk_rlines);
              elev(rlines) = sar_data.elev(:,1:num_next_chunk_rlines);
            else
              next_chunk_failed_flag = true;
            end
          else
            % If file does not exist, then give up loading previous chunk
            % data
            next_chunk_failed_flag = true;
          end
          
        end
      end
    end
  end
  
  if prev_chunk_failed_flag
    % Remove prev chunk data
    for ml_idx = 1:length(ml_list)
      data{ml_idx} = data{ml_idx}(:,num_prev_chunk_rlines+1:end,:,:,:);
      lat = lat(num_prev_chunk_rlines+1:end);
      lon = lon(num_prev_chunk_rlines+1:end);
      elev = elev(num_prev_chunk_rlines+1:end);
      for wf_adc_idx = 1:size(wf_adc_list,1)
        fcs{ml_idx}{wf_adc_idx}.origin = fcs{ml_idx}{wf_adc_idx}.origin(:,num_prev_chunk_rlines+1:end);
        fcs{ml_idx}{wf_adc_idx}.x = fcs{ml_idx}{wf_adc_idx}.x(:,num_prev_chunk_rlines+1:end);
        fcs{ml_idx}{wf_adc_idx}.y = fcs{ml_idx}{wf_adc_idx}.y(:,num_prev_chunk_rlines+1:end);
        fcs{ml_idx}{wf_adc_idx}.z = fcs{ml_idx}{wf_adc_idx}.z(:,num_prev_chunk_rlines+1:end);
        fcs{ml_idx}{wf_adc_idx}.roll = fcs{ml_idx}{wf_adc_idx}.roll(num_prev_chunk_rlines+1:end);
        fcs{ml_idx}{wf_adc_idx}.pitch = fcs{ml_idx}{wf_adc_idx}.pitch(num_prev_chunk_rlines+1:end);
        fcs{ml_idx}{wf_adc_idx}.heading = fcs{ml_idx}{wf_adc_idx}.heading(num_prev_chunk_rlines+1:end);
        fcs{ml_idx}{wf_adc_idx}.gps_time = fcs{ml_idx}{wf_adc_idx}.gps_time(num_prev_chunk_rlines+1:end);
        fcs{ml_idx}{wf_adc_idx}.surface = fcs{ml_idx}{wf_adc_idx}.surface(num_prev_chunk_rlines+1:end);
        fcs{ml_idx}{wf_adc_idx}.bottom = fcs{ml_idx}{wf_adc_idx}.bottom(num_prev_chunk_rlines+1:end);
        fcs{ml_idx}{wf_adc_idx}.pos = fcs{ml_idx}{wf_adc_idx}.pos(:,num_prev_chunk_rlines+1:end);
      end
    end
    num_prev_chunk_rlines = 0;
  end
  
  if next_chunk_failed_flag
    % Remove next chunk data
    for ml_idx = 1:length(ml_list)
      data{ml_idx} = data{ml_idx}(:,1:end-num_next_chunk_rlines,:,:,:);
      lat = lat(1:end-num_next_chunk_rlines);
      lon = lon(1:end-num_next_chunk_rlines);
      elev = elev(1:end-num_next_chunk_rlines);
      for wf_adc_idx = 1:size(wf_adc_list,1)
        fcs{ml_idx}{wf_adc_idx}.origin = fcs{ml_idx}{wf_adc_idx}.origin(:,1:end-num_next_chunk_rlines);
        fcs{ml_idx}{wf_adc_idx}.x = fcs{ml_idx}{wf_adc_idx}.x(:,1:end-num_next_chunk_rlines);
        fcs{ml_idx}{wf_adc_idx}.y = fcs{ml_idx}{wf_adc_idx}.y(:,1:end-num_next_chunk_rlines);
        fcs{ml_idx}{wf_adc_idx}.z = fcs{ml_idx}{wf_adc_idx}.z(:,1:end-num_next_chunk_rlines);
        fcs{ml_idx}{wf_adc_idx}.roll = fcs{ml_idx}{wf_adc_idx}.roll(1:end-num_next_chunk_rlines);
        fcs{ml_idx}{wf_adc_idx}.pitch = fcs{ml_idx}{wf_adc_idx}.pitch(1:end-num_next_chunk_rlines);
        fcs{ml_idx}{wf_adc_idx}.heading = fcs{ml_idx}{wf_adc_idx}.heading(1:end-num_next_chunk_rlines);
        fcs{ml_idx}{wf_adc_idx}.gps_time = fcs{ml_idx}{wf_adc_idx}.gps_time(1:end-num_next_chunk_rlines);
        fcs{ml_idx}{wf_adc_idx}.surface = fcs{ml_idx}{wf_adc_idx}.surface(1:end-num_next_chunk_rlines);
        fcs{ml_idx}{wf_adc_idx}.bottom = fcs{ml_idx}{wf_adc_idx}.bottom(1:end-num_next_chunk_rlines);
        fcs{ml_idx}{wf_adc_idx}.pos = fcs{ml_idx}{wf_adc_idx}.pos(:,1:end-num_next_chunk_rlines);
      end
    end
    num_next_chunk_rlines = 0;
  end
  
  if isfield(param.combine,'ft_over_sample') && ~isempty(param.combine.ft_over_sample)
    % param.combine.ft_over_sample should be a positive integer
    for ml_idx = 1:length(data)
      data{ml_idx} = interpft(data{ml_idx},size(data{ml_idx},1) * param.combine.ft_over_sample);
    end
    for wf = 1:length(sar_data.wfs)
      sar_data.wfs(wf).fs = sar_data.wfs(wf).fs * param.combine.ft_over_sample;
      sar_data.wfs(wf).dt = 1/sar_data.wfs(wf).fs;
      sar_data.wfs(wf).Nt = sar_data.wfs(wf).Nt * param.combine.ft_over_sample;
      sar_data.wfs(wf).df = sar_data.wfs(wf).fs / sar_data.wfs(wf).Nt;
      sar_data.wfs(wf).time = sar_data.wfs(wf).time(1) + sar_data.wfs(wf).dt*(0:sar_data.wfs(wf).Nt-1).';
      sar_data.wfs(wf).freq = sar_data.wfs(wf).fc ...
        + sar_data.wfs(wf).df * ifftshift( -floor(sar_data.wfs(wf).Nt/2) : floor((sar_data.wfs(wf).Nt-1)/2) ).';
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
    lever_arm_param.gps_source = sar_data.param_records.gps_source;
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
    W_ideal = 1 + ceil(tau_max / sar_data.wfs(wf).dt);
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
    Hwin = sar_data.param_csarp.csarp.ft_wind(length(sar_data.wfs(wf).time));
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
      = sar_data.wfs(wf).dt/Mt * (-Hwin_num_samp:Hwin_num_samp);
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
      ops_param.cmd.frms = param.load.frm + (-1:1);
      param.combine.bin_restriction = opsLoadLayers(ops_param, param.combine.bin_restriction);
      for idx=1:2
        % Interpolate layers onto GPS time of loaded data
        param.combine.bin_restriction(idx).bin = interp1(param.combine.bin_restriction(idx).gps_time, ...
          param.combine.bin_restriction(idx).twtt, sar_data.fcs.gps_time);
        % Convert from two way travel time to range bins
        param.combine.bin_restriction(idx).bin = interp1(sar_data.wfs(wf).time, ...
          1:length(sar_data.wfs(wf).time), param.combine.bin_restriction(idx).bin);
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
            ops_param.cmd.frms = param.load.frm + (-1:1);
            param.combine.doa_constraints(src_idx).layer = opsLoadLayers(ops_param, doa_res.params);
            % Interpolate layers onto GPS time of loaded data
            param.combine.doa_constraints(src_idx).layer.twtt = interp1(param.combine.doa_constraints(src_idx).layer.gps_time, ...
              param.combine.doa_constraints(src_idx).layer.twtt, sar_data.fcs.gps_time);
            % Ensure there is a value everywhere
            param.combine.doa_constraints(src_idx).layer.twtt = interp_finite(param.combine.doa_constraints(src_idx).layer.twtt);
        end
      end
    end
    
    param.combine.chan_equal = chan_equal;
    param.combine.fcs = fcs;
    array_param = param.combine;
    array_param.wfs = sar_data.wfs(wf);
    array_param.imgs = param.combine.imgs{img};
    if ~iscell(array_param.imgs)
      array_param.imgs = {array_param.imgs};
    end
    
    % Determine output range lines so that chunks will be seamless (this is
    % necessary because the output is decimated and the decimation may not
    % align with chunk lengths)
    first_rline = find(~mod(sar_out_rlines-1,param.combine.dline),1);
    array_param.rlines = num_prev_chunk_rlines + (first_rline : param.combine.dline : length(sar_out_rlines));
    array_param.rlines = array_param.rlines([1 end]);
    
    % Array Processing Function Call
    [array_param,tomo] = array_proc(array_param,data);
    param.combine.bins = array_param.bins;
    param.combine.lines = array_param.lines;
  end
  
  %% Debugging
  % -------------------------------------------------------------------
  if 0
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
  
  %% Save results
  % -----------------------------------------------------------------------
  array_fn = fullfile(combine_tmp_dir, sprintf('img_%02d_chk_%03d.mat', img, param.load.chunk_idx));
  fprintf('  Saving combine_wf_chan data %s (%s)\n', array_fn, datestr(now));
  combine_tmp_dir = fileparts(array_fn);
  if ~exist(combine_tmp_dir,'dir')
    mkdir(combine_tmp_dir);
  end
  
  Latitude = lat(1,array_param.lines);
  Longitude = lon(1,array_param.lines);
  Elevation = elev(1,array_param.lines);
  GPS_time = fcs{1}{1}.gps_time(array_param.lines);
  Surface = fcs{1}{1}.surface(array_param.lines);
  Bottom = fcs{1}{1}.bottom(array_param.lines);
  Roll = fcs{1}{1}.roll(array_param.lines);
  Pitch = fcs{1}{1}.pitch(array_param.lines);
  Heading = fcs{1}{1}.heading(array_param.lines);
  Data = tomo.val;
  Time = sar_data.wfs(wf_base).time(array_param.bins);
  param_records = sar_data.param_records;
  param_csarp = sar_data.param_csarp;
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

fprintf('%s done %s\n', mfilename, datestr(now));

success = true;
