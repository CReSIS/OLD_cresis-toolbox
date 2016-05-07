function success = combine_task(param)
% success = combine_task(param)
%
% Task/job running on cluster that is called from combine. This
% function is generally not called directly.
%
% param.combine = structure controlling the processing
%   Fields used in this function (and potentially array_proc.m)
%   .imgs = images to process
%   .method = 'csarp-combined', 'standard', 'mvdr', 'music', etc.
%   .debug_level = scalar integer specifying debug level
%   .three_dim.en = boolean, enable 3-D mode
%   .Nsv = used to compute sv
%   .sv = steering vector function handle (empty for default)
%   .in_path = input path string
%   .out_path = output path string
%   .chunk_ids = specific csarp output chunk IDs
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
physical_constants;

% =====================================================================
% Process data
% =====================================================================
if ~isfield(param.combine,'subaps') || isempty(param.combine.subaps)
  param.combine.subaps = {[1]};
end
if ~isfield(param.combine,'subbnds') || isempty(param.combine.subbnds)
  param.combine.subbnds = {[1]};
end

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
    
    [Latitude,Longitude,Elevation] ...
      = ecef2geodetic(sar_type_file.fcs.origin(1,array_param.lines), ...
      sar_type_file.fcs.origin(2,array_param.lines), ...
      sar_type_file.fcs.origin(3,array_param.lines),WGS84.ellipsoid);
    Latitude = Latitude*180/pi;
    Longitude = Longitude*180/pi;
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

success = true;

return;


