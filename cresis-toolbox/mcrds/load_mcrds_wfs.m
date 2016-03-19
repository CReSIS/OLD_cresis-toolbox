function [wfs,rec_data_size] = load_mcrds_wfs(records_wfs,param,adcs,proc_param)
% [wfs,rec_data_size] = load_mcrds_wfs(records_wfs,param,adcs,proc_param)
% 
% Creates a struct with all the waveform information.  This loads all
% waveforms (since it loads based off of param.radar) even when only
% a subset is being processed.  It also loads all the receivers that
% were loaded in the records file (and in the same order).
%
% It creates reference functions for each waveform for frequency domain
% pulse compression. These reference functions incorporate the
% param_radar.rx_path().td time delay corrections and the fast-time
% windowing that helps reduce sidelobes.
%
% INPUTS:
% records_wfs = This is the radar configuration information that comes
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
%    should have been in the radar header, but was not (i.e. records_wfs.wfs(1).wfs
%    should be sufficient, but because of poor design it is not)
%  .fs = sampling frequency (Hz)
%  .sample_size = number of bytes per data sample (e.g. 2 for 16 bits)
%  .Vpp_scale = full scale signal into ADC corresponds to this in volts
%     peak to peak (Vpp / full-scale-quantization)
%  .wfs = structure vector containing information about each radar
%    waveform
%   .Tpd = pulse duration of linear FM chirp (sec)
%   .f0 = start frequency of linear FM chirp (Hz)
%   .f1 = stop frequency of linear FM chirp (Hz)
%   .blank = time before which the waveform will be set to zero, leaving
%     blank is same as setting to -inf, i.e. turns feature off (sec)
%   .tx_weights = transmit antenna amplitude weights
%   .rx_paths = which receive paths were used for each ADC channel
%
% adcs = must be the same array of adcs passed to the create records
%   operation (i.e. records.param.records.file.adcs)
%
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
%
% rec_data_size = number of bytes in the data portion of each record
%
% Author: John Paden

if length(records_wfs.wfs) > 1
  warning('Code has limited support for multiple radar settings in one segment');
end

fs = records_wfs.wfs(1).wfs(1).fs;

num_sam = cell2mat({records_wfs.wfs(1).wfs.num_sam}).';
adc_en = cell2mat({records_wfs(1).wfs.wfs.adc_en});

rec_data_size = sum(adc_en * num_sam);
wf_sizes = cell2mat({records_wfs.wfs(1).wfs.num_sam});
wf_offsets = [];
for adc=1:8
  for wf=1:length(records_wfs.wfs(1).wfs)
    wf_offsets(adc,wf) = (adc-1)*sum(wf_sizes) + sum(wf_sizes(1:wf-1)) + 1;
  end
end

for wf = 1:length(records_wfs.wfs(1).wfs)
  % Assume all adc's match adc on some fields (like t0)
  adc = 1;

  wfs(wf).Tpd     = records_wfs.wfs(1).wfs(wf).Tpd;
  wfs(wf).f0      = records_wfs.wfs(1).wfs(wf).f0;
  wfs(wf).f1      = records_wfs.wfs(1).wfs(wf).f1;
  if isfield(param.radar.wfs(wf),'Tadc_adjust') && ~isempty(param.radar.wfs(wf).Tadc_adjust)
    Tadc_adjust = param.radar.wfs(wf).Tadc_adjust;
  else
    Tadc_adjust = 0;
  end
  if isfield(param.radar.wfs(wf),'Tadc') && ~isempty(param.radar.wfs(wf).Tadc)
    wfs(wf).t0    = param.radar.wfs(wf).Tadc + Tadc_adjust;
  else
    wfs(wf).t0    = records_wfs.wfs(1).wfs(wf).Tadc + Tadc_adjust;
  end
  if isfield(param.radar.wfs(wf),'blank')
    wfs(wf).blank   = param.radar.wfs(wf).blank;
  else
    wfs(wf).blank   = [];
  end
  wfs(wf).presums = records_wfs.wfs(1).wfs(wf).presums;
  wfs(wf).Nt_ref  = round(wfs(wf).Tpd * fs)+1;
  num_samples = records_wfs.wfs(1).wfs(wf).num_sam;
  if any(num_samples ~= num_samples(adc))
    error('Code does not handle different range gates for each ADC');
  end
  wfs(wf).Nt_raw  = num_samples(adc);
  if ~isfield(param.radar.wfs(wf),'zero_pad') || isempty(param.radar.wfs(wf).zero_pad)
    Nt_zero = wfs(wf).Nt_ref;
  else
    Nt_zero = round(param.radar.wfs(wf).zero_pad * fs)+1;
  end
  wfs(wf).Nt_pc   = wfs(wf).Nt_raw + Nt_zero - 1;
  wfs(wf).pad_length = wfs(wf).Nt_pc - wfs(wf).Nt_raw;
  wfs(wf).offset  = wf_offsets;

  if isfield(param.radar.wfs(wf),'BW') && ~isempty(param.radar.wfs(wf).BW)
    BW_window = param.radar.wfs(wf).BW(1);
    BW_decimation = param.radar.wfs(wf).BW(2);
  else
    BW_window = abs(wfs(wf).f1 - wfs(wf).f0);
    BW_decimation = abs(wfs(wf).f1 - wfs(wf).f0);
  end
  
  % ===================================================================
  % Quantization to Voltage conversion
  num_bit_shifts = max(0,ceil(log(wfs(wf).presums)/log(2)) - 4);
  wfs(wf).quantization_to_V ...
      = param.radar.Vpp_scale * 2^num_bit_shifts ...
      / (2^(param.radar.adc_bits)*wfs(wf).presums);

  % Other fields from param.radar files
  wfs(wf).tx_weights = param.radar.wfs(wf).tx_weights;
  wfs(wf).rx_paths = param.radar.wfs(wf).rx_paths;
  
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
  
  if isempty(param.radar.wfs(wf).ref_fn)
    ref_function = exp(1i*2*pi*f0*time + 1i*pi*alpha*time.^2);
  else
    param.day_seg = '';
    ref_fn = fullfile(ct_filename_out(param,'','ref_function',1), ...
      param.radar.wfs(wf).ref_fn);
    load(ref_fn);
    if length(ref_function) > length(time)
      ref_function = ref_function(1:length(time));
    else
      ref_function = [ref_function; zeros(length(time)-length(ref_function),1)];
    end
  end
  Htukeywin = tukeywin(Nt+2,param.radar.wfs(wf).tukey);
  Htukeywin = Htukeywin(2:end-1);
  
  if proc_param.ft_wind_time && ~isempty(proc_param.ft_wind)
    ref = Htukeywin.*proc_param.ft_wind(Nt).*ref_function;
  else
    ref = Htukeywin.*ref_function;
  end
 
  % Apply receiver delays to reference function
  Nt = wfs(wf).Nt_pc;
  df = 1/((Nt-1)*dt);
  freq = round(fc/fs)*fs + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
  for adc_idx = 1:length(adcs)
    adc = adcs(adc_idx);
    wfs(wf).ref{adc_idx} = conj(fft(ref,wfs(wf).Nt_pc) ...
      .* exp(-1i*2*pi*freq*param.radar.wfs(wf).Tsys(wfs(wf).rx_paths(adc))) );

    % Normalize reference function so that it is an estimator
    %  -- Accounts for pulse duration differences
    time_domain_ref = ifft(wfs(wf).ref{adc_idx});
    wfs(wf).ref{adc_idx} = wfs(wf).ref{adc_idx} ...
      ./ dot(time_domain_ref,time_domain_ref);
  end
  
  % ===================================================================
  % Create raw data with zero padding freq axes variables
  Nt = wfs(wf).Nt_pc;
  df = 1/((Nt-1)*dt);
  freq = round(fc/fs)*fs + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
  wfs(wf).time_raw = t0 + (0:dt:(Nt-1)*dt).';
  
  if proc_param.ft_dec
    wfs(wf).freq_inds = ifftshift(find(freq >= fc-BW_decimation/2 & freq <= fc+BW_decimation/2));
    wfs(wf).dc_shift = freq(wfs(wf).freq_inds(1))-fc;
    wfs(wf).fc = fc;
  else
    wfs(wf).freq_inds = 1:length(freq);
    wfs(wf).dc_shift = 0;
    wfs(wf).fc = fs*floor(max(f0,f1)/fs);
  end
  if ~proc_param.ft_wind_time && ~isempty(proc_param.ft_wind) ...
      && proc_param.pulse_comp
    ft_wind = zeros(size(freq));
    freq_inds = ifftshift(find(freq >= fc-BW_window/2 & freq <= fc+BW_window/2));
    ft_wind(freq_inds) = ifftshift(proc_param.ft_wind(length(freq_inds)));
    double_window_flag = 0;
    if double_window_flag 
        ft_wind = ft_wind.*ft_wind;
    end
    for adc_idx = 1:length(adcs)
      wfs(wf).ref{adc_idx} = wfs(wf).ref{adc_idx} .* ft_wind;
    end
  end
  
  % ===================================================================
  % Create output data time/freq axes variables
  
  if proc_param.pulse_comp
    Nt = wfs(wf).Nt_pc;
  else
    Nt = wfs(wf).Nt_raw;
  end
  Nt = length(wfs(wf).freq_inds);
  wfs(wf).dt = 1/((Nt-1)*df);
  wfs(wf).df = df;
  wfs(wf).Nt = Nt;
  wfs(wf).fs = (Nt-1)*df;
  dt = 1/wfs(wf).fs;
  if proc_param.ft_dec
    % Starts at fc goes to fc+BW/2, fc-BW/2 to fc
    wfs(wf).freq = fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    wfs(wf).time = t0 + (0:dt:(Nt-1)*dt).';
  else
    % Let ftnz = fast time nyquist zone
    % Starts at ftnz*fs goes to ftnz*fs+fs/2, ftnz*fs-fs/2 to ftnz*fs
    wfs(wf).freq = round(fc/fs)*fs ...
      + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    wfs(wf).time = t0 + (0:dt:(Nt-1)*dt).';
  end
  if proc_param.pulse_comp
    % Assuming pulse compression zero pads the front of the waveform, the
    % output will start earlier by an ammount proportional to the zero
    % padding.
    wfs(wf).time = wfs(wf).time - wfs(wf).pad_length / fs;
  end

end
