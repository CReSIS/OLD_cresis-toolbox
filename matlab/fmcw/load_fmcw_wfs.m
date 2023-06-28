function wfs = load_fmcw_wfs(settings, param, adcs, proc_param)
% wfs = load_fmcw_wfs(settings, param, adcs, proc_param)
%
%
% It creates reference functions for each waveform for frequency domain
% pulse compression. These reference functions incorporate the
% param.radar.Tsys time delay corrections and the fast-time
% windowing that helps reduce sidelobes.
%
% Nc = Number of channels
% Nw = Number of waveforms
%
% INPUTS:
% settings = This is the radar configuration information that comes
%   from the data files themselves (or is overridden by the param
%   spreadsheet).
%  .wfs
%   .num_sam = number of samples recorded for this waveform
%   .presums = number of presums or coherent averages
%   .bit_shifts = how much each data sample was divided by so that it would
%     fit inside 16 bits (important when the number of presums is large
%     enough)
%   .Tadc = time of first sample relative to the start of the transmission (sec)
%   .rx_gain = receiver gain (linear power)   
%  .nyquist_zone (optional)
%  .loopback_mode (optional)
%
% .param = structure from parameter spreadsheet
%  .radar_name = Some system timing information is determined from this
%    variable
%  .radar = This is all the radar and processing configuration information
%     passed in from the parameter spreadsheet.
%   .fs = sampling frequency (Hz)
%   .Vpp_scale = full scale signal into ADC corresponds to this in volts
%      peak to peak (Vpp / full-scale-quantization)
%   .wfs = Nw length structure vector containing information about each radar
%     waveform
%    .f0 = start frequency of linear FM chirp (Hz)
%    .f1 = stop frequency of linear FM chirp (Hz)
%    .fmult = multiplication factor on f0 and f1
%    .fLO = local oscillator to apply after fmult
%    .Tpd = pulse duration of linear FM chirp (sec)
%    .Tsys = system time delay, RF side (sec)
%    .record_start_idx = start index of sampling
%    .presum_override = number of hardware presums
%    .Tadc = ADC time delay, IF side (sec)
%    .nyquist_zone = which Nyquist zone the data is in (zero-indexed)
%    .loopback_mode = 0: normal, 1: loopback
%    .good_rbins = which range bins to use (one-indexed to what was
%      actually recorded)
%    .tx_weights = transmit antenna amplitude weights
%    .rx_paths = Nc length vector holding receiver path for each ADC
%    .adc_gains = Nc length vector holding receiver gain for each ADC
%
% adcs = not used
%
% proc_param = not used
%
% OUTPUTS:
% wfs = structure vector containing information for each waveform
%   transmitted.
%   .fs_raw = raw data sampling frequency (Hz)
%   .f0 = start frequency of chirp (Hz)
%   .f1 = stop frequency of chirp (Hz)
%   .fmult = multiplication factor on f0 and f1
%   .fLO = local oscillator to apply after fmult
%   .Tpd = pulse duration (sec)
%   .Tsys = system time delay, RF side (sec)
%   .record_start_idx = start index of sampling
%   .presum_override = number of hardware presums
%   .Tadc = ADC time delay, IF side (sec)
%   .nyquist_zone = which Nyquist zone the data is in (zero-indexed)
%   .loopback_mode = 0: normal, 1: loopback
%   .quantization_to_V = multiplication factor to turn ADC quantization
%     into volts
%   .tx_weights = vector of transmit array weights (see corresponding
%      lever_arm function), linear scale
%   .rx_paths = vector of length equal to the number of adcs
%      the index into the array is the absolute adc, and the value
%      at that position tells which receive antenna is connected to that adc
%   .adc_gains = transmit array weights, linear scale
%      the index into the array is the absolute adc, and the value
%      at that position tells the amplifier gain into the ADC
%   .fc = RF center frequency (Hz)
%   .chirp_rate = RF chirp rate (Hz/sec)
%   .good_rbins = which recorded samples to use
%   .Nt_raw = length of raw data in fast-time after good rbins
%   .fs = sampling frequency after pulse compression (Hz)
%   .time = RF time axis after pulse compression (sec)
%
% rec_data_size = number of bytes in the data portion of each record
%
% Author: John Paden

for wf = 1:length(param.radar.wfs)
  wfs(wf).fs_raw            = param.radar.fs;
  wfs(wf).f0                = param.radar.wfs(wf).f0;
  wfs(wf).f1                = param.radar.wfs(wf).f1;
  wfs(wf).Tpd               = param.radar.wfs(wf).Tpd;
  wfs(wf).fmult             = param.radar.wfs(wf).fmult;
  
  wfs(wf).fc = abs((wfs(wf).f1 + wfs(wf).f0)*wfs(wf).fmult/2 + param.radar.wfs(wf).fLO);
  
  if isfield(settings,'nyquist_zone')
    wfs(wf).nyquist_zone    = settings.nyquist_zone;
  elseif ~isempty(param.radar.wfs(wf).nyquist_zone)
    % Override nyquist zone
    wfs(wf).nyquist_zone    = param.radar.wfs(wf).nyquist_zone;
  else
    wfs(wf).nyquist_zone    = [];
  end

  wfs(wf).rx_paths = param.radar.wfs(wf).rx_paths;
  wfs(wf).tx_weights = param.radar.wfs(wf).tx_weights;
    
  if param.records.file_version == 1
    wfs(wf).num_sam = settings.wfs.num_sam;
    wfs(wf).presums = settings.wfs.presums;
    if isfield(param.radar.wfs(wf),'presum_override') ...
        && ~isempty(param.radar.wfs(wf).presum_override)
      wfs(wf).presums = param.radar.wfs(wf).presum_override;
    end
  end
  
  %% Determine Chirp Rate
  % Ideal Chirp Rate:
  %   BW = (param.radar.f1-param.radar.f0)*param.radar.fmult;
  %   chirp_rate = BW / param.radar.Tpd;
  % Actual Chirp Rate:
  %   This takes into account the actual hardware implementation
  wfs(wf).chirp_rate = (wfs(wf).f1-wfs(wf).f0)/wfs(wf).Tpd * wfs(wf).fmult;
  
end

return;


















% sample_size: number of bytes in each radar sample
sample_size = 2;

wf_offsets = cumsum([0; settings.wfs(1).num_sam(1:end-1)]) ...
  * sample_size;
rec_data_size = sum(settings.wfs(1).num_sam) ...
  * sample_size;

for wf = 1:length(param.radar.wfs)
    
  %% Input Checking
  if ~isfield(param.radar.wfs(wf),'record_start_idx')
    param.radar.wfs(wf).record_start_idx = [];
  end
  if isnan(settings.wfs(wf).start_idx) && isempty(param.radar.wfs(wf).record_start_idx)
    error('start_idx field from file is NaN, must specify in param spreadsheet');
  end
  if isnan(settings.wfs(wf).presums) && isempty(param.radar.wfs(wf).presum_override)
    error('presums field from file is NaN, must specify in param spreadsheet');
  end
  if isnan(settings.wfs(wf).Tadc) && isempty(param.radar.wfs(wf).Tadc)
    error('t0 field from file is NaN, must specify in param spreadsheet');
  end
  
  wfs(wf).fs_raw            = param.radar.fs;
  wfs(wf).f0                = param.radar.wfs(wf).f0;
  wfs(wf).f1                = param.radar.wfs(wf).f1;
  wfs(wf).fmult             = param.radar.wfs(wf).fmult;
  wfs(wf).fLO               = param.radar.wfs(wf).fLO;
  wfs(wf).Tpd               = param.radar.wfs(wf).Tpd;
  wfs(wf).Tsys              = param.radar.wfs(wf).Tsys; % System delay (RF side)
  
  if isempty(param.radar.wfs(wf).record_start_idx)
    wfs(wf).record_start_idx  = settings.wfs.start_idx;
  else
    wfs(wf).record_start_idx  = param.radar.wfs(wf).record_start_idx;
  end
  if isempty(param.radar.wfs(wf).presum_override)
    wfs(wf).presums         = settings.wfs.presums;
  else
    wfs(wf).presums         = param.radar.wfs(wf).presum_override;
  end
  if isempty(param.radar.wfs(wf).Tadc)
    wfs(wf).Tadc         = settings.wfs.Tadc; % ADC delay (IF side)
  else
    wfs(wf).Tadc         = param.radar.wfs(wf).Tadc; % ADC delay (IF side)
  end
  
  if isempty(param.radar.wfs(wf).nyquist_zone)
    wfs(wf).nyquist_zone    = settings.nyquist_zone;
  else
    wfs(wf).nyquist_zone    = param.radar.wfs(wf).nyquist_zone;
  end
  if isempty(param.radar.wfs(wf).loopback_mode)
    wfs(wf).loopback_mode   = settings.loopback_mode;
  else
    wfs(wf).loopback_mode   = param.radar.wfs(wf).loopback_mode;
  end
  
  % ===================================================================
  % Quantization to Voltage conversion
  wfs(wf).bit_shifts = settings.wfs.bit_shifts(wf);
  wfs(wf).quantization_to_V ...
    = param.radar.Vpp_scale * 2^settings.wfs.bit_shifts(wf) ...
    / (2^param.radar.adc_bits*wfs(wf).presums);
  
  % Other fields from param.radar files
  wfs(wf).tx_weights = param.radar.wfs(wf).tx_weights;
  wfs(wf).rx_paths = param.radar.wfs(wf).rx_paths;
  wfs(wf).adc_gains = param.radar.wfs(wf).adc_gains;
  
  wfs(wf).fc = (wfs(wf).f1 + wfs(wf).f0)*wfs(wf).fmult/2;
  
  % ===================================================================
  % Create output data time/freq axes variables
  
  %% Determine Chirp Rate
  % Ideal Chirp Rate:
  %   BW = (param.radar.f1-param.radar.f0)*param.radar.fmult;
  %   chirp_rate = BW / param.radar.Tpd;
  % Actual Chirp Rate:
  %   This takes into account how the IDL software rounds the radar control
  %   software GUI entries into a 32 bit value to program the DDS
  wfs(wf).chirp_rate = floor((wfs(wf).f1-wfs(wf).f0)/wfs(wf).Tpd * (4/1e9) / 1e9 * 2^32) * 1e9/2^32/(4/1e9) * wfs(wf).fmult;
  
  %% Determine the samples to use in the FFT
  if isempty(param.radar.wfs(wf).good_rbins)
    if any(strcmpi(param.radar_name,{'kuband','snow'}))
      %% Determine good samples automatically
      max_range_delay = param.radar.fs/2 * (param.radar.wfs(wf).nyquist_zone+1) / abs(wfs(wf).chirp_rate);
      if isempty(param.radar.wfs(wf).record_start_idx)
        % Use start index stored in the data file
        start_idx = settings.wfs.Tadc/param.radar.fs; % Start recording index set in radar control software
      else
        % Override start index in the data file
        start_idx = param.radar.wfs(wf).record_start_idx; % Start recording index set in radar control software
      end
      % Based on Ben Panzer's report 1U_DAQ_simple_timing_analysis.docx
      wfs(wf).good_rbins(1) = max(1,2 + ceil((146.85 + max_range_delay*wfs(wf).fs_raw) - start_idx));
      bad_bins_at_end_of_record = 8;
      wfs(wf).good_rbins(2) = min(settings.wfs.num_sam - bad_bins_at_end_of_record, ...
        floor(wfs(wf).good_rbins(1) + (wfs(wf).Tpd-max_range_delay)*wfs(wf).fs_raw));
    elseif any(strcmpi(param.radar_name,{'kuband2','snow2'}))
      %% Determine good samples automatically
      max_range_delay = wfs(wf).fs_raw/2 * (wfs(wf).nyquist_zone+1) / abs(wfs(wf).chirp_rate);
      stop_idx = settings.wfs.start_idx+settings.wfs.num_sam - 8;
      if isempty(param.radar.wfs(wf).record_start_idx)
        % Use start index stored in the data file
        start_idx = settings.wfs.start_idx; % Start recording index set in radar control software
      else
        % Override start index in the data file
        start_idx = param.radar.wfs(wf).record_start_idx; % Start recording index set in radar control software
      end
      % MADE THESE NUMBERS UP!!! AN ANALYSIS NEEDS TO BE DONE!!!
      wfs(wf).good_rbins(1) = max(1,2 + ceil((146.85 + max_range_delay*wfs(wf).fs_raw) - (start_idx/2+231)));
      wfs(wf).good_rbins(2) = min(stop_idx-start_idx, ...
        floor(wfs(wf).good_rbins(1) + (wfs(wf).Tpd-max_range_delay)*wfs(wf).fs_raw));
    elseif any(strcmpi(param.radar_name,{'kaband3','kuband3','snow3'}))
      %% Determine good samples automatically
      max_range_delay = wfs(wf).fs_raw/2 * (wfs(wf).nyquist_zone+1) / abs(wfs(wf).chirp_rate);
      stop_idx = settings.wfs.start_idx+settings.wfs.num_sam - 8;
      if isempty(param.radar.wfs(wf).record_start_idx)
        % Use start index stored in the data file
        start_idx = settings.wfs.start_idx; % Start recording index set in radar control software
      else
        % Override start index in the data file
        start_idx = param.radar.wfs(wf).record_start_idx; % Start recording index set in radar control software
      end
      % MADE THESE NUMBERS UP!!! AN ANALYSIS NEEDS TO BE DONE!!!
      wfs(wf).good_rbins(1) = max(1,2 + ceil((146.85 + max_range_delay*wfs(wf).fs_raw) - (start_idx/2+231)));
      wfs(wf).good_rbins(2) = min(stop_idx-start_idx, ...
        floor(wfs(wf).good_rbins(1) + (wfs(wf).Tpd-max_range_delay)*wfs(wf).fs_raw));
    end
  else
    wfs(wf).good_rbins = param.radar.wfs(wf).good_rbins;
  end
  wfs(wf).Nt_raw = wfs(wf).good_rbins(2) - wfs(wf).good_rbins(1) + 1;
  wfs(wf).fs = abs(wfs(wf).chirp_rate)*wfs(wf).Nt_raw/wfs(wf).fs_raw;
  wfs(wf).Nt = ceil((wfs(wf).Nt_raw+1)/2) - 2;
  
  dt = 1/wfs(wf).fs;
  wfs(wf).time = (0:dt:(wfs(wf).Nt-1)*dt).' + wfs(wf).nyquist_zone*wfs(wf).Nt*dt - param.radar.wfs(wf).Tsys;
  
end
