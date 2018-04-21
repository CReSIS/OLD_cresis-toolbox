function [wfs,state] = data_load_wfs(param, records)
% [wfs,state] = data_load_wfs(param, records)
%
% Loads waveform/adc information into wfs
% Loads state information for loading raw data into state
% Both are required for data_load.m
%
% param: structure with parameter information
%  .records.file_version: raw file version
%  .load.imgs: cell vector of imgs (each entry is a separate wf_adc_list)
%    Each wf_adc_list is an Nx2 array where N is the number of channels,
%    the first column is the waveform, and the second column is the adc.
%    Absolute index of wf and adc are used in this array.
%  .load.ft_wind: function handle to fast time window that will be applied
%    in the frequency domain
% wfs: 
% state: structure vector that helps data_load load data (state of loader)
%   Each entry in the structure vector corresponds to a specific board.
%     accum(board)
%   Each entry contains lists of how to accumulate and store data for the
%   board when loading from load_mcords2_data.  load_mcords2_data first
%   accumulates the data (presums):
%     accum(board).data{accumulator-instance}
%   Once the presums are finished, the data is processed and stored in
%   the output variable:
%     g_data{img}(fast-time,slow-time,wf_adc_idx)
%   The three fields in accum(board) are all length Kx1 where K is the number
%   of accumulator instances.  The number of accumulator instances is
%   determined by the length of the fields.
%  .adc = K x 1 vector indicating which adc this instance is pulled from
%  .wf = K x 1 vector indicating which waveform this instance is pulled from
%  .wf_adc_idx = an index in the final output array unless adcs are
%    combined in which case size(g_data{:}, 3) == 1.
%  .img = an index in the final output array
%
% Examples: At the bottom of this file
%
% Author: John Paden
%
% See also: data_load.m

save('/tmp/data_load_wfs_test.mat')
keyboard
load('/tmp/data_load_wfs_test.mat')

%% Build raw data loading "state" structure
% =========================================================================

% Create list of wf-adc pairs to determine which boards to load
wf_list = [];
adc_list = [];
board_list = [];
for img = 1:length(param.load.imgs)
  for wf_adc_idx= 1:size(param.load.imgs{img},1)
    for adc_column = 2:2:size(param.load.imgs{img},2)
      wf_list(end+1) = param.load.imgs{img}(wf_adc_idx,adc_column-1);
      adc_list(end+1) = param.load.imgs{img}(wf_adc_idx,adc_column);
      board_list(end+1) = adc_to_board(param.radar_name,adc_list(end));
    end
  end
end
boards = unique(board_list);

% Populate state structure
state = [];
for board = boards
  state(board).adc = [];
  state(board).wf = [];
  state(board).wf_adc_idx = [];
  state(board).img = [];
  state(board).wf_adc_sum = [];
  state(board).wf_adc_sum_done = [];
  state(board).img_comb_idx = [];
  for img = 1:length(param.load.imgs) % For each image img
    for wf_adc_idx = 1:size(param.load.imgs{img},1) % For ach wf-adc pair
      for adc_column = 2:2:size(param.load.imgs{img},2) % For each combined wf-adc pair
        wf = param.load.imgs{img}(wf_adc_idx,adc_column-1); % wf stored in odd columns
        adc = param.load.imgs{img}(wf_adc_idx,adc_column); % adc stored in even columns
        if ~isfield(param.radar.wfs(wf),'wf_adc_sum') || isempty(param.radar.wfs(wf).wf_adc_sum)
          % if wf_adc_sum not specied, then this wf-adc pair is the only
          % one to load
          wf_adc_sum = [wf adc 1];
        else
          % if wf_adc_sum is specified, then extract the list out, each
          % row specifies a wf-adc pair and a complex weight
          % For zero-pi mod with two wf-adc pairs, the weight should be -0.5 or 0.5
          % For IQ mod with two wf-adc pairs, the weights should be 1, j, -j, or -1
          %   [wf adc weight] 
          wf_adc_sum = param.radar.wfs(wf).wf_adc_sum{adc};
        end
        for wf_adc_sum_idx = 1:size(wf_adc_sum,1)
          wf = wf_adc_sum(wf_adc_sum_idx,1);
          adc = wf_adc_sum(wf_adc_sum_idx,2);
          % Add wf-adc pair to state list
          state(board).adc(end+1) = adc;
          state(board).wf(end+1) = wf;
          state(board).wf_adc_idx(end+1) = wf_adc_idx;
          state(board).img(end+1) = img;
          state(board).wf_adc_sum(end+1) = wf_adc_sum(wf_adc_sum_idx,3);
          state(board).wf_adc_sum_done(end+1) = false;
          state(board).img_comb_idx(end+1) = adc_column/2;
        end
        state(board).wf_adc_sum_done(ends) = true;
      end
    end
  end
end

%% Create wfs structure with waveform information
% =========================================================================
[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);
adcs = 1:max(param.records.file.adcs);
for wf = 1:length(param.radar.wfs)
  
  %% Input checks
  % =======================================================================
  if isfield(param.radar.wfs(wf),'Tpd') && ~isempty(param.radar.wfs(wf).Tpd)
    wfs(wf).Tpd = param.radar.wfs(wf).Tpd;
  elseif any(param.records.file_version == [405 406 410]) % [acords mcords]
    wfs(wf).Tpd = records.settings.wfs(1).wfs(wf).Tpd(1);
  end
  if isfield(param.radar.wfs(wf),'fLO') && ~isempty(param.radar.wfs(wf).fLO)
    wfs(wf).fLO = param.radar.wfs(wf).fLO;
  else
    wfs(wf).fLO = 0;
  end
  if isfield(param.radar.wfs(wf),'fmult') && ~isempty(param.radar.wfs(wf).fmult)
    wfs(wf).fmult = param.radar.wfs(wf).fmult;
  else
    wfs(wf).fmult = 1;
  end
  if isfield(param.radar.wfs(wf),'f0') && ~isempty(param.radar.wfs(wf).f0)
    wfs(wf).f0 = param.radar.wfs(wf).f0;
  elseif any(param.records.file_version == [405 406 410]) % [acords mcords]
    wfs(wf).f0 = records.settings.wfs(1).wfs(wf).f0(1);
  end
  wfs(wf).f0 = wfs(wf).f0*wfs(wf).fmult + wfs(wf).fLO;
  if isfield(param.radar.wfs(wf),'f1') && ~isempty(param.radar.wfs(wf).f1)
    wfs(wf).f1 = param.radar.wfs(wf).f1;
  elseif any(param.records.file_version == [405 406 410]) % [acords mcords]
    wfs(wf).f1 = records.settings.wfs(1).wfs(wf).f1(1);
  end
  wfs(wf).f1 = wfs(wf).f1*wfs(wf).fmult + wfs(wf).fLO;
  if isfield(param.radar.wfs(wf),'Tadc_adjust') && ~isempty(param.radar.wfs(wf).Tadc_adjust)
    wfs(wf).Tadc_adjust = param.radar.wfs(wf).Tadc_adjust;
  else
    wfs(wf).Tadc_adjust = 0;
  end
  if isfield(param.radar.wfs(wf),'Tadc') && ~isempty(param.radar.wfs(wf).Tadc)
    wfs(wf).t0_raw    = param.radar.wfs(wf).Tadc + wfs(wf).Tadc_adjust;
  elseif any(param.records.file_version == [405 406 410]) % [acords mcords]
    wfs(wf).t0_raw    = records.settings.wfs(1).wfs(wf).t0(1) + wfs(wf).Tadc_adjust;
  elseif isfield(records.settings.wfs(wf),'t0')
    wfs(wf).t0_raw    = records.settings.wfs(wf).t0 + wfs(wf).Tadc_adjust;
  else
    wfs(wf).t0_raw    = 0 + wfs(wf).Tadc_adjust;
  end
  if isfield(param.radar.wfs(wf),'blank') && ~isempty(param.radar.wfs(wf).blank)
    wfs(wf).blank   = param.radar.wfs(wf).blank;
  else
    wfs(wf).blank   = -inf;
  end
  if isfield(param.radar.wfs(wf),'DDC_mode') && ~isempty(param.radar.wfs(wf).DDC_mode)
    wfs(wf).DDC_mode   = param.radar.wfs(wf).DDC_mode;
  else
    wfs(wf).DDC_mode   = 0;
  end
  if param.records.file_version == 410 % mcords
    wfs(wf).fs_raw = records_wfs.wfs(1).wfs(1).fs;
  else
    if wfs(wf).DDC_mode == 0
      wfs(wf).fs_raw = param.radar.fs;
    else
      wfs(wf).fs_raw = param.radar.fs / 2^(1+wfs(wf).DDC_mode);
    end
  end
  if isfield(param.radar.wfs(wf),'zero_pad') && ~isempty(param.radar.wfs(wf).zero_pad)
    wfs(wf).zero_pad   = param.radar.wfs(wf).zero_pad;
  else
    wfs(wf).zero_pad   = 0;
  end
  if isfield(param.radar.wfs(wf),'ft_dec') && ~isempty(param.radar.wfs(wf).ft_dec)
    wfs(wf).ft_dec = param.radar.wfs(wf).ft_dec;
  else
    if strcmpi(radar_type,'fmcw')
      wfs(wf).ft_dec = [1 1];
    else
      [numerator denominator] = rat((wfs(wf).f1 - wfs(wf).f0) / wfs(wf).fs_raw);
      wfs(wf).ft_dec = [numerator denominator];
    end
  end
  if isfield(param.radar.wfs(wf),'DDC_freq') && ~isempty(param.radar.wfs(wf).DDC_freq)
    wfs(wf).DDC_freq   = param.radar.wfs(wf).DDC_freq;
  else
    wfs(wf).DDC_freq   = 0;
  end
  if isfield(param.radar.wfs(wf),'presums') && ~isempty(param.radar.wfs(wf).presums)
    wfs(wf).presums = param.radar.wfs(wf).presums;
  elseif any(param.records.file_version == [1:8 405 406 410]) % [acords mcords]
    wfs(wf).presums = records.settings.wfs(1).wfs(wf).presums(1);
  else
    wfs(wf).presums = records.settings.wfs(wf).presums;
  end
  if isfield(param.radar.wfs(wf),'BW') && ~isempty(param.radar.wfs(wf).BW)
    wfs(wf).BW_window = param.radar.wfs(wf).BW(1);
  else
    wfs(wf).BW_window = abs(wfs(wf).f1 - wfs(wf).f0);
  end
  if any(param.records.file_version == [405 406 410]) % [acords mcords]
    wfs(wf).Nt_raw = records.settings.wfs(1).wfs(wf).num_sam(1);
  elseif isfield(records.settings.wfs(wf),'num_sam')
    wfs(wf).Nt_raw = records.settings.wfs(wf).num_sam(1);
  else
    wfs(wf).Nt_raw = 0;
  end
  if isfield(param.radar.wfs(wf),'conjugate') && ~isempty(param.radar.wfs(wf).conjugate)
    wfs(wf).conjugate   = param.radar.wfs(wf).conjugate;
  else
    wfs(wf).conjugate   = 0;
  end
  if isfield(param.radar.wfs(wf),'ft_wind_time') && ~isempty(param.radar.wfs(wf).ft_wind_time)
    wfs(wf).ft_wind_time   = param.radar.wfs(wf).ft_wind_time;
  else
    wfs(wf).ft_wind_time   = [];
  end
  if isfield(param.radar.wfs(wf),'tukey') && ~isempty(param.radar.wfs(wf).tukey)
    wfs(wf).tukey   = param.radar.wfs(wf).tukey;
  else
    wfs(wf).tukey   = 0;
  end
  if isfield(param.radar.wfs(wf),'DC_adjust') && ~isempty(param.radar.wfs(wf).DC_adjust)
    tmp = load(fullfile(ct_filename_out(param,'noise','',1), ...
      param.radar.wfs(wf).DC_adjust),'DC_adjust');
    for adc_idx = 1:length(adcs)
      adc = adcs(adc_idx);
      wfs(wf).DC_adjust(adc_idx) = tmp.DC_adjust(adc);
    end
  else
    wfs(wf).DC_adjust   = zeros(size(adcs));
  end
  if isfield(param.radar.wfs(wf),'gain_fn') && ~isempty(param.radar.wfs(wf).gain_fn)
    for adc = adcs
      gain_fn_name = char(param.radar.wfs(wf).gain_fn);
      gain_fn_name = regexprep(gain_fn_name,'%w',sprintf('%.0f',wf));
      gain_fn_name = regexprep(gain_fn_name,'%a',sprintf('%.0f',adc));
      gain_fn = fullfile(ct_filename_out(param,'noise','',1), [gain_fn_name '.mat']);
      
      wfs(wf).gain(adc) = load(gain_fn);
    end
  end
  if isfield(settings,'nyquist_zone')
    wfs(wf).nyquist_zone    = settings.nyquist_zone;
  elseif ~isempty(param.radar.wfs(wf).nyquist_zone)
    % Override nyquist zone
    wfs(wf).nyquist_zone    = param.radar.wfs(wf).nyquist_zone;
  else
    wfs(wf).nyquist_zone    = [];
  end
  wfs(wf).tx_weights = param.radar.wfs(wf).tx_weights;
  wfs(wf).rx_paths = param.radar.wfs(wf).rx_paths;
  wfs(wf).adc_gains = param.radar.wfs(wf).adc_gains;

  %% Compute supporting variables
  % =======================================================================
  wfs(wf).chirp_rate = (wfs(wf).f1-wfs(wf).f0) / wfs(wf).Tpd;
  wfs(wf).fc = (wfs(wf).f1+wfs(wf).f0)/2;
  
  %% Quantization to Voltage conversion
  % =======================================================================
  if any(param.records.file_version == [405 406]) % acords
    num_bit_shifts = records.settings.wfs(1).wfs(wf).bit_shifts(1);
  elseif param.records.file_version == 410 % mcords
    num_bit_shifts = max(0,ceil(log(wfs(wf).presums)/log(2)) - 4);
  elseif isfield(records.settings.wfs(wf),'bit_shifts')
    num_bit_shifts = records.settings.wfs(wf).bit_shifts;
  else
    num_bit_shifts = 0;
  end
  wfs(wf).quantization_to_V ...
    = param.radar.Vpp_scale * 2^num_bit_shifts ...
    / (2^param.radar.adc_bits*wfs(wf).presums);

  
  if strcmpi(radar_type,'fmcw')
    %% FMCW: Create time and frequency axis information
    % =====================================================================
  
  elseif strcmpi(radar_type,'pulsed')
    %% Pulsed: Create time and frequency axis information
    % =====================================================================
    dt = 1/wfs(wf).fs_raw;
    wfs(wf).time_raw = wfs(wf).t0_raw + dt*(0:wfs(wf).Nt_raw-1).';

    nz0 = floor((wfs(wf).f0-wfs(wf).DDC_freq)/wfs(wf).fs_raw*2)
    nz1 = floor((wfs(wf).f1-wfs(wf).DDC_freq)/wfs(wf).fs_raw*2)
    
    df = wfs(wf).fs_raw/wfs(wf).Nt_raw;
    if nz0 == nz1 && wfs(wf).DDC_mode == 0
      % Assume real sampling since signal does not cross Nyquist boundary
      if mod(nz0,2)
        % Negative frequencies first since this is an odd Nyquist zone
        wfs(wf).freq_raw = floor(nz0/2)*wfs(wf).fs_raw + df*(0:wfs(wf).Nt_raw-1);
        wfs(wf).freq_raw(1:floor(wfs(wf).Nt_raw/2)) ...
          = wfs(wf).freq_raw(1:floor(wfs(wf).Nt_raw/2)) - floor(nz0/2)*wfs(wf).fs_raw - ceil(nz0/2)*wfs(wf).fs_raw;
      else
        % Positive frequencies first since this is an odd Nyquist zone
        wfs(wf).freq_raw = floor(nz0/2)*wfs(wf).fs_raw + df*(0:wfs(wf).Nt_raw-1);
        wfs(wf).freq_raw(end-floor(wfs(wf).Nt_raw/2)+1:end) ...
          = wfs(wf).freq_raw(end-floor(wfs(wf).Nt_raw/2)+1:end) - (nz0/2+1)*wfs(wf).fs_raw - nz0/2*wfs(wf).fs_raw;
      end
    else
      % Assume complex sampling since signal crosses Nyquist boundary
      wfs(wf).freq_raw = wfs(wf).DDC_freq ...
        + ifftshift( -floor(wfs(wf).Nt_raw/2)*df : df : floor((wfs(wf).Nt_raw-1)/2)*df ).';
      % Shift the frequencies so they are centered on the signal bandwidth
      wfs(wf).freq_raw = wfs(wf).freq_raw - floor((wfs(wf).freq_raw - (wfs(wf).fc-wfs(wf).fs_raw/2))/wfs(wf).fs_raw)*wfs(wf).fs_raw;
    end
    
    wfs(wf).fs = wfs(wf).fs_raw * wfs(wf).ft_dec(1)/wfs(wf).ft_dec(2);
    
    Nt_ref  = floor(wfs(wf).Tpd * wfs(wf).fs) + 1;
    wfs(wf).pad_length = Nt_ref + wfs(wf).zero_pad - 1;
    wfs(wf).Nt = ceil(wfs(wf).Nt_pc*wfs(wf).ft_dec(1)/wfs(wf).ft_dec(2)) ...
      + wfs(wf).pad_length;
    
    
    % Starts at fc goes to fc+fs/2, fc-fs/2 to fc
    df = wfs(wf).fs/wfs(wf).Nt;
    wfs(wf).freq = wfs(wf).fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    
    % Starts same as raw data minus the padding which is added at the front
    dt = 1/wfs(wf).fs;
    wfs(wf).time = wfs(wf).t0_raw - dt*wfs(wf).pad_length + dt*(0:wfs(wf).Nt-1).';
    
    % Modify reference function so that time vector elements are multiples
    % of dt.
    wfs(wf).time_correction = dt - mod(wfs(wf).time(1),dt);
    wfs(wf).time = wfs(wf).time + wfs(wf).time_correction;
    
    %% Pulsed: Create reference function
    % =====================================================================
    BW = wfs(wf).f1 - wfs(wf).f0;
    time = (0:dt:(wfs(wf).Nt-1)*dt).';
    ref = tukeywin_cont(time/wfs(wf).Tpd-0.5)*exp(1i*2*pi*-BW/2*time + 1i*pi*wfs(wf).chirp_rate*time.^2);
    if ~isempty(wfs(wf).ft_wind_time)
      ref = wfs(wf).ft_wind_time(wfs(wf).Nt_ref).*ref;
    end
    
    for adc = adcs
      ref_fn_name = char(param.radar.wfs(wf).ref_fn);
      ref_fn_name = regexprep(ref_fn_name,'%w',sprintf('%.0f',wf));
      ref_fn_name = regexprep(ref_fn_name,'%a',sprintf('%.0f',adc));
      ref_fn = fullfile(ct_filename_out(param,'noise','',1), [ref_fn_name '.mat']);
      
      if isempty(ref_fn_name) || ~exist(ref_fn,'file')
        wfs(wf).ref{adc} = conj(fft(ref,Nt) ...
          .* exp(-1i*2*pi*freq*param.radar.wfs(wf).Tsys(param.radar.wfs(wf).rx_paths(adc))) );
        wfs(wf).ref_windowed(adc) = false;
        
      else
        % Load reference function from collate_deconv.m (e.g. for deconvolution)
        load(ref_fn,'ref_nonnegative','ref_negative','ref_windowed','ref_window');
        ref_Nt = length(ref_nonnegative)+length(ref_negative);
        if ref_Nt > Nt
          error('Window in ref_fn %s is longer than time axis, increase zero padding to use this reference function or shorten the reference function', ref_fn);
        end
        ref_from_file = [ref_nonnegative; zeros(Nt-ref_Nt,1); ref_negative];
        wfs(wf).ref_windowed(adc) = ref_windowed;
        
        if ref_windowed && ~isequal(ref_window,proc_param.ft_wind)
          error('Window in ref_fn %s is %s and does not match param.load.ft_wind %s', ...
            ref_fn, func2str(ref_window), func2str(proc_param.ft_wind));
        end
        
        ref_from_file = ref_from_file ./ abs(max(ref_from_file));
        wfs(wf).ref{adc} = conj(fft(ref_from_file,Nt) ...
          .* exp(-1i*2*pi*freq*param.radar.wfs(wf).Tsys(param.radar.wfs(wf).rx_paths(adc))) );
      end
      
      if ~wfs(wf).ref_windowed(adc)
        wfs(wf).ref{adc} = wfs(wf).ref{adc} .* ft_wind;
      end
      wfs(wf).ref{adc} = wfs(wf).ref{adc} .* exp(1i*2*pi*freq*wfs(wf).time_correction);

      % Normalize reference function so that it is an estimator
      %  -- Accounts for pulse duration differences
      time_domain_ref = ifft(wfs(wf).ref{adc});
      wfs(wf).ref{adc} = wfs(wf).ref{adc} ...
        ./ dot(time_domain_ref,time_domain_ref);
    end
  end
end

%% Populate the waveform/adc offsets into each record
% =========================================================================

% offset: bytes of data before this data channel
% rec_data_size: number of bytes of data in the record
if any(param.records.file_version == [405 406]) % [acords]
  wf_num_sam = cell2mat({settings.wfs(1).wfs.num_sam}).';
elseif any(param.records.file_version == [403 407 408]) % [mcords3 mcords5]
  wf_num_sam = cell2mat({settings.wfs.num_sam}).';
end

if any(param.records.file_version == [1]) % [fmcw1]
  wfs(wf).num_sam = settings.wfs.num_sam;
end

if any(param.records.file_version == [403 405 406 407 408]) % [acords]
  wf = 1;
  wf_offsets(wf) = 0;
  if ~wfs(wf).DDC_mode
    % Real data
    if any(strcmpi(radar_name,{'acords'}))
      num_elem = length(settings.wfs(wf).wfs(wf).adc_gains(1,:));
      rec_data_size = wf_num_sam(wf)*num_elem*sample_size;
    else
      rec_data_size = wf_num_sam(wf)*sample_size;
    end
  else
    % Complex data
    rec_data_size = 2*wf_num_sam(wf)*sample_size;
  end
  
  for wf = 2:length(param.radar.wfs)
    if isfield(param.radar.wfs(wf),'DDC_mode') && ~isempty(param.radar.wfs(wf).DDC_mode)
      wfs(wf).DDC_mode   = param.radar.wfs(wf).DDC_mode;
    else
      wfs(wf).DDC_mode   = 0;
    end
    
    if ~wfs(wf).DDC_mode
      % Real data
      if any(strcmpi(radar_name,{'acords'}))
        num_elem = length(settings.wfs(wf).wfs(wf).adc_gains(1,:));
        wf_offsets(wf) = wf_offsets(wf-1) + wf_num_sam(wf-1)*sample_size;
        rec_data_size = rec_data_size + wf_num_sam(wf)*num_elem*sample_size;
      else
        wf_offsets(wf) = wf_offsets(wf-1) + wf_num_sam(wf-1)*sample_size;
        rec_data_size = rec_data_size + wf_num_sam(wf)*sample_size;
      end
    else
      % Complex data
      wf_offsets(wf) = wf_offsets(wf-1) + 2*wf_num_sam(wf-1)*sample_size;
      rec_data_size = rec_data_size + 2*wf_num_sam(wf)*sample_size;
    end
  end
  
elseif any(param.records.file_version == [102]) % [accum2]
  wf_offsets = cumsum([0; wf_num_sam(1:end-1)]) ...
    * sample_size;
  rec_data_size = sum(wf_num_sam) * sample_size;
end
