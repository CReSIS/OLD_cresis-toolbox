function [wfs,rec_data_size] = load_mcords_wfs(settings, param, adcs, proc_param)
% [wfs,rec_data_size] = load_mcords_wfs(settings, param, adcs, proc_param)
%
% Creates a struct with all the waveform information.  This loads all
% waveforms (since it loads based off of param.radar) even when only
% a subset is being processed.  It also loads all the receivers that
% were loaded in the records file (and in the same order).
%
% It creates reference functions for each waveform for frequency domain
% pulse compression. These reference functions incorporate the
% param.radar.rx_path().td time delay corrections and the fast-time
% windowing that helps reduce sidelobes.
%
% INPUTS:
% settings = This is the radar configuration information that comes
%   from the data files themselves.
%  .num_sam = number of samples recorded for this waveform
%  .presums = number of presums or coherent averages
%  .which_bits = how much each data sample was divided by so that it would
%    fit inside 16 bits (important when the number of presums is large
%    enough)
%  .t0 = time of first sample relative to the start of the transmission (sec)
%
% .param
%  .radar_name = Used by ct_filename_out
%  .season_name = Used by ct_filename_out
%  .radar = This is all the radar configuration information that
%     should have been in the radar header, but was not (i.e. records_wfs
%     should be sufficient, but because of poor design it is not)
%   .fs = sampling frequency (Hz)
%   .sample_size = number of bytes per data sample (e.g. 2 for 16 bits)
%   .Vpp_scale = full scale signal into ADC corresponds to this in volts
%      peak to peak (Vpp / full-scale-quantization)
%   .wfs = structure vector containing information about each radar
%     waveform
%    .Tpd = pulse duration of linear FM chirp (sec)
%    .f0 = start frequency of linear FM chirp (Hz)
%    .f1 = stop frequency of linear FM chirp (Hz)
%    .ref_fn = filename of reference function file (passed to ct_filename_out)
%       empty uses ideal chirp
%    .tukey = time domain tukey window parameter, tukeywin(Nt,?), to apply
%    .blank = time before which the waveform will be set to zero, leaving
%      blank is same as setting to -inf, i.e. turns feature off (sec)
%    .tx_weights = transmit antenna amplitude weights
%    .rx_paths = which receive paths were used for each ADC channel
%    .rx_path(rx).td = per (waveform,adc) pair time delay
%     correction
%
% adcs = must be an array of the form 1:max(param.records.file.adcs)
%   The param.records.file.adcs MUST BE THE ONE USED TO GENERATE THE
%   RECORDS.   param.records.file.adcs is param_records.file.adcs after
%   it is stored in the records file
%
% proc_param = structure typically passed in from param.get_heights,
%   param.csarp, or param.load_mcords. It needs the following fields:
%  .ft_wind_time = boolean, to apply window in time domain
%    (usually false, but is possible to be true with chirp signals)
%  .ft_wind = function handle to window in fast time
%    (e.g. @hanning, @boxcar, @blackman, @tukeywin, etc)
%  .pulse_comp = boolean, to apply pulse compression
%  .ft_dec = boolean, to apply fast time decimation
%
% OUTPUTS:
% wfs = structure vector containing information for each waveform
%   transmitted.
%   .Tpd = pulse duration (sec)
%   .f0 = start frequency of chirp (Hz)
%   .f1 = stop frequency of chirp (Hz)
%   .t0 = start time of recording relative to tx pulse (sec)
%   .blank = start and stop of receiver blank. This time gate will be
%     set to zero by load_mcords_data.  Leave empty if you want to
%     disable this operation. (sec)
%   .Nt_ref = length of reference in fast-time
%   .Nt_raw = length of raw data in fast-time
%   .Nt_pc = length of pulse compressed data in fast-time
%   .offset = byte offset in the data part of the record
%   .ref = pulse compression reference waveform (conj. freq domain)
%      cell array, index is the absolute adc number
%   .time_raw = time axis of raw data (sec)
%   .freq_inds = frequency inds that will be kept after fast time decimation
%   .dc_shift = DC shift caused by fast time decimation (Hz)
%   .fc = center frequency (Hz)
%   .dt = sample spacing of output data product (sec)
%   .df = frequency spacing of output data product (Hz)
%   .Nt = number of samples in output data product
%   .fs = sampling frequency of output data product (Hz)
%   .freq = frequency axis of output data product (Hz)
%   .time = time axis of output data product (sec)
%   .tx_weights = vector of transmit array weights (see corresponding
%      lever_arm function), linear scale
%   .rx_paths = vector of length equal to the number of adcs
%      the index into the array is the absolute adc, and the value
%      at that position tells which receive antenna is connected to that adc
%   .adc_gains = transmit array weights, linear scale
%      the index into the array is the absolute adc, and the value
%      at that position tells the amplifier gain into the ADC
%
% rec_data_size = number of bytes in the data portion of each record
%
% Author: John Paden

if ~isfield(proc_param,'wf_adc_comb')
  proc_param.wf_adc_comb.en = 0;
end

if param.records.file_version == 406 % ACORDS ver 2
  sample_size = 4;
else
  sample_size = 2;
end

% wf_offsets = bytes of data  before each waveform in a record. This does
%   not include any header bytes.
% rec_data_size = number of bytes of data (not including any header
%   bytes).
if any(strcmpi(param.radar_name,{'acords'}))
  wf_num_sam = cell2mat({settings.wfs(1).wfs.num_sam}).';
else
  wf_num_sam = cell2mat({settings.wfs.num_sam}).';
end
if any(strcmpi(param.radar_name,{'acords','hfrds','mcords','mcords2','mcords3','mcords4','mcords5','seaice'}))
  wf = 1;
  wf_offsets(wf) = 0;
  if isfield(param.radar.wfs(wf),'DDC_mode') && ~isempty(param.radar.wfs(wf).DDC_mode)
    wfs(wf).DDC_mode   = param.radar.wfs(wf).DDC_mode;
  else
    wfs(wf).DDC_mode   = 0;
  end
  if ~wfs(wf).DDC_mode
    % Real data
    if any(strcmpi(param.radar_name,{'acords'}))
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
      if any(strcmpi(param.radar_name,{'acords'}))
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
  
elseif strcmpi(param.radar_name,'accum2')
  wf_offsets = cumsum([0; wf_num_sam(1:end-1)]) ...
    * sample_size;
  rec_data_size = sum(wf_num_sam) * sample_size;
end

if proc_param.wf_adc_comb.en
  % Determine the waveform with the longest time record. All waveforms will
  % be increased to this length for pulse compression to guarantee alignment
  % of the frequency bins. This zero padding will be placed at the end of the
  % waveform.
  
  for wf = 1:length(param.radar.wfs)
    if isfield(param.radar.wfs(wf),'DDC_mode') && ~isempty(param.radar.wfs(wf).DDC_mode)
      wfs(wf).DDC_mode   = param.radar.wfs(wf).DDC_mode;
    else
      wfs(wf).DDC_mode   = 0;
    end
    if wfs(wf).DDC_mode == 0
      fs = param.radar.fs;
    else
      fs = param.radar.fs / 2^(1+wfs(wf).DDC_mode);
    end
    
    Tpd     = param.radar.wfs(wf).Tpd;
    Nt_ref  = floor(Tpd * fs) + 1;
    Nt_raw  = settings.wfs(wf).num_sam;
    Nt_pc(wf)   = Nt_raw + Nt_ref - 1;
  end
  
  Nt_pc_max = max(Nt_pc);
end

for wf = 1:length(param.radar.wfs)
  adc_idx = 1;
  
  if any(strcmpi(param.radar_name,{'acords'}))
    settings.wfs(wf).bit_shifts = settings.wfs(1).wfs(wf).bit_shifts(1);
    settings.wfs(wf).presums = settings.wfs(1).wfs(wf).presums(1);
    settings.wfs(wf).num_sam = settings.wfs(1).wfs(wf).num_sam(1);
    settings.wfs(wf).t0 = settings.wfs(1).wfs(wf).t0(1);
    settings.wfs(wf).prf = settings.wfs(1).wfs(wf).prf(1);
    settings.wfs(wf).f0 = settings.wfs(1).wfs(wf).f0(1);
    settings.wfs(wf).f1 = settings.wfs(1).wfs(wf).f1(1);
    settings.wfs(wf).wf_gen_clk = settings.wfs(1).wfs(wf).wf_gen_clk(1);
    settings.wfs(wf).daq_clk = settings.wfs(1).wfs(wf).daq_clk(1);
    settings.wfs(wf).Tpd = settings.wfs(1).wfs(wf).Tpd(1);
    settings.wfs(wf).tx_win = settings.wfs(1).wfs(wf).tx_win(1);
    settings.wfs(wf).blank = settings.wfs(1).wfs(wf).blank(1);
    settings.wfs(wf).adc_gains = settings.wfs(1).wfs(wf).adc_gains(1,:);
    if param.records.file_version == 406
      settings.wfs(wf).elem_slots = settings.wfs(1).wfs(wf).elem_slots(1,:);
    end
  end
  
  wfs(wf).Tpd     = param.radar.wfs(wf).Tpd;
  wfs(wf).f0      = param.radar.wfs(wf).f0;
  wfs(wf).f1      = param.radar.wfs(wf).f1;
  if isfield(param.radar.wfs(wf),'Tadc_adjust') && ~isempty(param.radar.wfs(wf).Tadc_adjust)
    Tadc_adjust = param.radar.wfs(wf).Tadc_adjust;
  else
    Tadc_adjust = 0;
  end
  if isfield(param.radar.wfs(wf),'Tadc') && ~isempty(param.radar.wfs(wf).Tadc)
    wfs(wf).t0    = param.radar.wfs(wf).Tadc + Tadc_adjust;
  else
    wfs(wf).t0    = settings.wfs(wf).t0 + Tadc_adjust;
  end
  if isfield(param.radar.wfs(wf),'blank') && ~isempty(param.radar.wfs(wf).blank)
    wfs(wf).blank   = param.radar.wfs(wf).blank;
  else
    wfs(wf).blank   = [];
  end
  if isfield(param.radar.wfs(wf),'DDC_mode') && ~isempty(param.radar.wfs(wf).DDC_mode)
    wfs(wf).DDC_mode   = param.radar.wfs(wf).DDC_mode;
  else
    wfs(wf).DDC_mode   = 0;
  end
  if wfs(wf).DDC_mode == 0
    fs = param.radar.fs;
  else
    fs = param.radar.fs / 2^(1+wfs(wf).DDC_mode);
  end
  if isfield(param.radar.wfs(wf),'DDC_freq') && ~isempty(param.radar.wfs(wf).DDC_freq)
    wfs(wf).DDC_freq   = param.radar.wfs(wf).DDC_freq;
  else
    wfs(wf).DDC_freq   = 0;
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
  wfs(wf).presums = settings.wfs(wf).presums;
  wfs(wf).Nt_ref  = floor(wfs(wf).Tpd * fs) + 1;
  wfs(wf).Nt_raw  = settings.wfs(wf).num_sam;
  wfs(wf).Nt_pc   = wfs(wf).Nt_raw + wfs(wf).Nt_ref - 1;
  wfs(wf).pad_length = wfs(wf).Nt_pc - wfs(wf).Nt_raw;
  wfs(wf).offset  = wf_offsets(wf);
  if proc_param.wf_adc_comb.en
    wfs(wf).Nt_pc = Nt_pc_max;
  end

  if isfield(param.radar.wfs(wf),'BW') && ~isempty(param.radar.wfs(wf).BW)
    BW_window = param.radar.wfs(wf).BW(1);
    BW_decimation = param.radar.wfs(wf).BW(2);
  else
    BW_window = abs(wfs(wf).f1 - wfs(wf).f0);
    BW_decimation = abs(wfs(wf).f1 - wfs(wf).f0);
  end
  
  % ===================================================================
  % Quantization to Voltage conversion
  wfs(wf).quantization_to_V ...
    = param.radar.Vpp_scale * 2^settings.wfs(wf).bit_shifts ...
    / (2^14*wfs(wf).presums);
  
  % Other fields from param.radar files
  wfs(wf).tx_weights = param.radar.wfs(wf).tx_weights;
  wfs(wf).rx_paths = param.radar.wfs(wf).rx_paths;
  wfs(wf).adc_gains = param.radar.wfs(wf).adc_gains;
  
  % ===================================================================
  % Create reference waveform
  Tpd = wfs(wf).Tpd;
  f0  = wfs(wf).f0;
  f1  = wfs(wf).f1;
  Nt  = wfs(wf).Nt_ref;
  t0  = wfs(wf).t0;
  
  dt = 1/fs;
  BW = f1-f0;
  time = (0:dt:(Nt-1)*dt).';
  alpha = BW / Tpd;
  fc = (f0 + f1)/2;
  
  if wfs(wf).DDC_mode == 0
    %% DDC Disabled or No DDC
    ref_function = exp(1i*2*pi*f0*time + 1i*pi*alpha*time.^2);
  else
    %% DDC Enabled
    ref_function = exp(1i*2*pi*(f0 - wfs(wf).DDC_freq)*time + 1i*pi*alpha*time.^2);
  end
  if any(strcmpi(param.radar_name,{'acords'}))
    Htukeywin = hamming(Nt);
  else
    Htukeywin = tukeywin(Nt+2,param.radar.wfs(wf).tukey);
    Htukeywin = Htukeywin(2:end-1);
  end
  if proc_param.ft_wind_time && ~isempty(proc_param.ft_wind)
    ref = Htukeywin.*proc_param.ft_wind(Nt).*ref_function;
  else
    ref = Htukeywin.*ref_function;
  end
  
  % Apply receiver delays to reference function
  Nt = wfs(wf).Nt_pc;
  df = 1/(Nt*dt);
  if wfs(wf).DDC_mode == 0
    freq = fs*floor(fc/fs) + (0:df:(Nt-1)*df).';
  else
    freq = wfs(wf).DDC_freq + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
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
        error('Window in ref_fn %s does not match ft_wind parameter', ref_fn);
      end
      
      ref_from_file = ref_from_file ./ abs(max(ref_from_file));
      wfs(wf).ref{adc} = conj(fft(ref_from_file,Nt) ...
        .* exp(-1i*2*pi*freq*param.radar.wfs(wf).Tsys(param.radar.wfs(wf).rx_paths(adc))) );
    end
    
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
  
  % ===================================================================
  % Create raw data with zero padding freq axes variables
  if ~proc_param.pulse_comp
    Nt = wfs(wf).Nt_raw;
  end
  
  df = 1/(Nt*dt);
  %freq = round(fc/fs)*fs + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
  freq = fs*floor(fc/fs) + (0:df:(Nt-1)*df).';
  wfs(wf).time_raw = t0 + (0:dt:(Nt-1)*dt).';
  
  if proc_param.ft_dec
    if wfs(wf).DDC_mode == 0
      %% DDC Disabled or no DDC
      wfs(wf).freq_inds = ifftshift(find(freq >= fc-BW_decimation/2 & freq <= fc+BW_decimation/2));
      wfs(wf).dc_shift = freq(wfs(wf).freq_inds(1))-fc;
    else
      %% DDC Enabled
      freq = wfs(wf).DDC_freq + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
      wfs(wf).freq_inds = find(freq >= min(f0,f1) & freq <= max(f0,f1));
      wfs(wf).dc_shift = fc-freq(wfs(wf).freq_inds(1));
    end
    wfs(wf).fc = fc;
  else
    wfs(wf).freq_inds = 1:length(freq);
    wfs(wf).dc_shift = 0;
    wfs(wf).fc = fs*floor(max(f0,f1)/fs);
  end
  if ~proc_param.ft_wind_time && ~isempty(proc_param.ft_wind) ...
      && proc_param.pulse_comp
    ft_wind = zeros(size(freq));
    
    if wfs(wf).DDC_mode == 0
      %% DDC Disabled or no DDC
      freq_inds = ifftshift(find(freq >= fc-BW_window/2 & freq <= fc+BW_window/2));
    else
      %% DDC Enabled
      freq_inds = find(freq >= min(f0,f1) & freq <= max(f0,f1));
    end
    [~,sorted_freq_inds] = sort(freq(freq_inds));
    sorted_freq_inds = freq_inds(sorted_freq_inds);
    ft_wind(sorted_freq_inds) = proc_param.ft_wind(length(freq_inds));
    
    for adc = adcs
      if ~wfs(wf).ref_windowed(adc)
        wfs(wf).ref{adc} = wfs(wf).ref{adc} .* ft_wind;
      end
    end
  end
  
  for adc = adcs
    % Normalize reference function so that it is an estimator
    %  -- Accounts for pulse duration differences
    time_domain_ref = ifft(wfs(wf).ref{adc});
    wfs(wf).ref{adc} = wfs(wf).ref{adc} ...
      ./ dot(time_domain_ref,time_domain_ref);
  end
  
  % ===================================================================
  % Create output data time/freq axes variables
  
  Nt = length(wfs(wf).freq_inds);
  wfs(wf).dt = 1/(Nt*df);
  wfs(wf).df = df;
  wfs(wf).Nt = Nt;
  wfs(wf).fs = Nt*df;
  dt = 1/wfs(wf).fs;
  if proc_param.ft_dec
    % Starts at fc goes to fc+BW/2, fc-BW/2 to fc
    wfs(wf).freq = fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    wfs(wf).time = t0 + (0:dt:(Nt-1)*dt).';
  else
    % Let ftnz = fast time nyquist zone
    % Starts at ftnz*fs goes to ftnz*fs+fs/2, ftnz*fs-fs/2 to ftnz*fs
    %     wfs(wf).freq = round(fc/fs)*fs ...
    %       + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    wfs(wf).freq = fs*floor(fc/fs) + (0:df:(Nt-1)*df).';
    wfs(wf).time = t0 + (0:dt:(Nt-1)*dt).';
  end
  if proc_param.pulse_comp
    % Assuming pulse compression zero pads the front of the waveform, the
    % output will start earlier by an ammount proportional to the zero
    % padding.
    wfs(wf).time = wfs(wf).time - wfs(wf).pad_length / fs;
    
    % Modify reference function so that time vector elements are multiples
    % of dt.
    wfs(wf).time_correction = dt - mod(wfs(wf).time(1),dt);
    wfs(wf).time = wfs(wf).time + wfs(wf).time_correction;
    
    for adc = adcs
      wfs(wf).ref{adc} = wfs(wf).ref{adc} .* exp(1i*2*pi*freq*wfs(wf).time_correction);
    end
  end
  
end
