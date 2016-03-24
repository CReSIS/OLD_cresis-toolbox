function [data,time,freq] = pulse_compress(data,param)
% [data,time,freq] = pulse_compress(data,param)
%
% Generic pulse compression routine. 
%
% data = data to be pulse compressed
%   first dimension is fast-time and is pulse compressed
% param = structure describing pulse compression
%   .f0: start frequency, Hz
%   .f1: stop frequency, Hz
%   .Tpd: pulse duration, sec
%   .time: column vector for time axis of data, sec
%   .tukey: time domain Tukey window parameter: tukeywin(N,?),
%     default is 0
%   .td_window_func: function handle to time domain window function
%     default is to use tukeywin with param.tukey (setting param.tukey to
%     zero or leaving it undefined has the same effect as doing no
%     windowing). This window is applied to the reference.
%   .window_func: function handle to frequency domain window function
%     default is not to window. This window is applied to the reference.
%   .Tsys: time delay to add into pulse compression operation
%     A positive value of Tsys will have the affect of moving the pulse
%     compressed output to earlier time. i.e. This much delay will be
%     removed from the input signal.
%   .zero_pad: scalar double
%     Negative value: sets zero padding to this length (zeros added to
%       start of waveform). This is useful when the decimation filter
%       creates artifacts in the data which is usually only a problem when
%       boxcar or near boxcar filters are used like tukeywin. In these
%       cases, setting a value to -param.Tpd*C for C > 1 helps avoid
%       circular convolution.
%     Zero: no zero padding
%     Positive: zero padding set to -param.Tpd (i.e. zero padding equal to
%       pulse duration)
%   .BW: 1x2 vectors specifying the bandwidth to use for
%     windowing and decimation/interpolation respectively. For small time
%     bandwidth products it is helpful to increase the windowing bandwidth
%     beyond abs(param.f1-param.f0). Decimation/interpolation allows the
%     pulse compression to build in under/over sampling. The default for
%     both bandwidths is abs(param.f1-param.f0).
%   .decimate: Same as param.baseband. Parameter is misnamed and
%     should use param.baseband.
%   .baseband: logical scalar to complex baseband the signal or not.
%   .Mt = scalar time-domain over-sampling factor, default it 1
%   .DDC_mode: logical scalar to treat data as DDC (assumes the data
%     has been complex basebanded already)
%   .DDC_freq: required if DDC_mode is true, double scalar
%     representing the center frequency which has been shifted to zero
%     frequency
%   .pulse_compress: true to apply the reference function filter to the
%     data or not (all other steps will be run including baseband of data)
%   .stc: sensitivity timing control fields. Is optional. If it is
%     specified then both of the following fields must be specified:
%     .tdelay: apply fast-time delay correction to data
%     .gain: apply fast-time gain correction to data
%   .deconv = deconvolution parameter struct [NOT DONE YET]
%     .freq_rng = [low high] frequencies to use
%     .time_rng = [low high] time range to use
%     .time_delay = remove this amount of time delay from deconvolution
%       waveform
%     .mode = 0 for apply deconvolution supplied by .ref
%        1: return reference waveform
%     .ref
%
% data = pulse compressed data
% time = new time axis, sec
%
% Example:
%  Test code at bottom of this file
%
% Author: John Paden
%
% See also run_load_data, run_basic_load_mcords, run_basic_load_mcords2

if ~isfield(param,'tukey') || isempty(param.tukey)
  param.tukey = 0;
end
if ~isfield(param,'zero_pad') || isempty(param.zero_pad) || param.zero_pad > 0
  param.zero_pad = -param.Tpd;
end
if ~isfield(param,'BW')
  param.BW = [abs(param.f1 - param.f0) abs(param.f1 - param.f0)];
end
if isfield(param,'decimate')
  param.baseband = param.decimate;
end
if ~isfield(param,'baseband') || isempty(param.baseband)
  param.baseband = true;
end
if ~isfield(param,'td') || isempty(param.td)
  param.td = 0;
end
if ~isfield(param,'Mt') || isempty(param.Mt)
  param.Mt = 1;
end
if ~isfield(param,'DDC_mode') || isempty(param.DDC_mode)
  param.DDC_mode = 0;
end
if ~isfield(param,'pulse_compress') || isempty(param.pulse_compress)
  param.pulse_compress = 1;
end
if isreal(data)
  real_data = true;
else
  real_data = false;
end

dt = param.time(2) - param.time(1);

%% Zero pad and FFT input data and create output time axis
Nt_zero = floor(-param.zero_pad/dt) + 1;
Nt_pc = length(param.time) + Nt_zero - 1;
time = -(Nt_zero-1)*dt + (0:dt:(Nt_pc-1)*dt).';
data = fft(cat(1,zeros(Nt_zero-1,size(data,2),size(data,3)), data),[],1);

%% Create reference function (matched filter)
BW = param.f1 - param.f0;
alpha = BW / param.Tpd;
Nt_ref = floor(param.Tpd/dt) + 1;

% Reference always starts at time = 0
ref = exp(1i*2*pi*param.f0*(time-time(1)) ...
  + 1i*pi*alpha*(time-time(1)).^2);
% Reference time domain window
if isfield(param,'td_window_func') && ~isempty(param.td_window_func)
  Hwin = param.td_window_func(Nt_ref);
else
  Hwin = tukeywin_trim(Nt_ref,param.tukey);
end

% Limit reference to pulse duration and apply time-domain window
ref(1:Nt_ref) = ref(1:Nt_ref) .* Hwin;
ref(Nt_ref+1:end) = 0;

%% Create frequency axis
fc = (param.f1 + param.f0)/2;
fs = 1/dt;
nyquist_zone = floor(fc/fs);
df = fs/Nt_pc;
if param.DDC_mode
  freq = param.DDC_freq + ifftshift( -floor(Nt_pc/2)*df : df : floor((Nt_pc-1)/2)*df ).';
  freq_inds = ifftshift(find(freq >= fc-BW/2 & freq <= fc+BW/2));
  ref = ref .* exp(-j*2*pi*param.DDC_freq.*time);
else
  freq = nyquist_zone*fs + (0:df:(Nt_pc-1)*df).';
  freq_inds = find(freq >= min(param.f0,param.f1) & freq <= max(param.f0,param.f1));
end
[~,sorted_freq_inds] = sort(freq(freq_inds));
sorted_freq_inds = freq_inds(sorted_freq_inds);
freq_inds = sorted_freq_inds;

%% Convert reference to frequency domain, adding in an optional time delay
ref = conj(fft(ref) .* exp(-1i*2*pi*freq*param.td));
ref2 = ref;

if isfield(param,'window_func') && ~isempty(param.window_func)
  if fc-param.BW(1)/2 < min(freq) || fc+param.BW(1)/2 > max(freq)
    warning('Windowing bandwidth is larger than original data');
  end
  ft_wind = zeros(size(freq));
  ft_wind(sorted_freq_inds) = param.window_func(length(freq_inds));
  ref = ref.*ft_wind;
end

% Apply match filter/pulse compression in Fourier domain
if param.baseband
  if fc-param.BW(1)/2 < min(freq) || fc+param.BW(1)/2 > max(freq)
    warning('Decimation bandwidth is larger than original data');
  end
  ref = ref(freq_inds);
  ref2 = ref2(freq_inds);
  data = data(freq_inds,:,:);
  Nt = length(freq_inds);
  freq = freq(freq_inds);
  fs = Nt*df;
  dt = 1/fs;
  time = time(1) + (0:dt:(Nt-1)*dt).';
  freq = ifftshift(freq);
  data = ifftshift(data,1);
  ref = ifftshift(ref);
  ref2 = ifftshift(ref2);
end

%% Normalize reference function
time_domain_ref = ifft(ref); % Reference with freq+time domain windowing
time_domain_ref2 = ifft(ref2); % Represents data with time domain windowing
if real_data
  ref = 2*ref ...
    ./ dot(time_domain_ref2,time_domain_ref);
else
  ref = ref ...
    ./ dot(time_domain_ref2,time_domain_ref);
end

% Adjust time axis for start time of data
time = param.time(1) + time;

if param.pulse_compress
  if isfield(param,'stc')
    % Apply sensitivity timing control (STC) corrections
    tmp_data = data;
    
    stc_tdelay = interp1(param.time, param.stc.tdelay, time);
    stc_tdelay = interp_finite(stc_tdelay);
    stc_gain = interp1(param.time, param.stc.gain, time);
    stc_gain = interp_finite(stc_gain);
    data = tmp_data;
    Mt = 100;
    data = ifft(data);
    Nt = length(time);
    data = interpft(data,Nt*Mt);
    time_Mt = time(1) + (time(2)-time(1))/Mt * (0:Nt*Mt-1).';
    data = interp1(time_Mt, data, time+stc_tdelay, 'linear', 0).*stc_gain;
    
    global data_h;
    figure(1); clf;
    plot(param.time,lp(data_h));
    hold on
    plot(time,lp(data));
    plot(time,lp(ifft(tmp_data)));
    hold off;
    
    figure(1); clf;
    plot(param.time,real(data_h));
    hold on
    plot(time,real(data));
    plot(time,real(ifft(tmp_data)));
    plot(time+4.5e-6,10000000*real(ifft(conj(ref))));
    hold off;
    keyboard
    
    data = fft(data);
  end
  for rline = 1:prod(size(data,2)*size(data,3))
    data(:,rline) = ifft(data(:,rline) .* ref);
  end
else
  for rline = 1:prod(size(data,2)*size(data,3))
    data(:,rline) = ifft(data(:,rline));
  end
end

% Adjust time axis for over-sampling
if param.Mt ~= 1
  Nt_oversample = round(Nt*param.Mt);
  if param.baseband
    data = interpft(data,Nt_oversample);
  else
    % Insert zeros as far from fc as possible
    error('Not supported');
  end
  time = time(1) + dt*Nt/Nt_oversample * (0:Nt_oversample-1);
end

return;

% =======================================================================
% Test code
% =======================================================================

%% Example #1
% Showing effects of start/stop frequency shifting to create delay
clear pc_param;
pc_param.f0 = 180e6;
pc_param.f1 = 210e6;
pc_param.Tpd = 10e-6;
pc_param.Mt = 10;
pc_param.baseband = true;
pc_param.window_func = @hanning;
pc_param.Mt = 100;
pc_param.tukey = 0.1;
fs = 1e9/9;
Nt = 20000;
t0 = 0e-6;
dt = 1/fs;
BW = pc_param.f1-pc_param.f0;
alpha = BW/pc_param.Tpd;
pc_param.time = t0 + (0:dt:(Nt-1)*dt).';
signal = cos(2*pi*pc_param.f0*pc_param.time + pi*alpha*pc_param.time.^2);
% Apply Tukey window
good_mask = pc_param.time>=0 & pc_param.time<=pc_param.Tpd;
signal(~good_mask) = 0;
signal(good_mask) = signal(good_mask) .* tukeywin_trim(sum(good_mask),pc_param.tukey);

[pc_signal,pc_time] = pulse_compress(signal,pc_param);

pc_param.f0 = 180e6 + 5e-9*alpha;
pc_param.f1 = 210e6 + 5e-9*alpha;
[pc_signal_shift,pc_time] = pulse_compress(signal,pc_param);

plot(pc_time*1e6, lp(pc_signal,2));
hold on;
plot(pc_param.time*1e6, lp(signal,2),'r');
plot(pc_time*1e6, lp(pc_signal_shift,2),'g');
hold off;
xlabel('Time (us)');
ylabel('Relative power (dB)');

%% Examples #2
% Accum: Loading raw data and pulse compressing it
fs = 1e9;
fn = '/mnt/backup-iu/array1/20131217/accum_20131217/accum2_00_20131217_025800_0028.bin';
[hdr,data] = basic_load_accum2(fn,struct('clk',fs));
param.f0 = 900000000;
param.f1 = 600000000;
param.Tpd = 2.048e-6;
param.time = 1/fs * (0:size(data{1},1)-1).';
param.tukey = 0.2;
param.window_func = @hanning;
wf = 1;
[pc_data,pc_time] = pulse_compress(data{wf} - repmat(mean(data{wf},2),[1 size(data{wf},2)]),param);
figure(1); clf;
imagesc(lp(pc_data));
figure(3); clf;
plot(lp(pc_data(:,9063)));

%% Example #3
% UWB MCoRDS4: Showing pulse compression with real only versus with complex
clear pc_param;
pc_param.f0 = 200e6;
pc_param.f1 = 450e6;
pc_param.Tpd = 10e-6;
pc_param.window_func = @hanning;
pc_param.tukey = 0.05;
fs = 1e9/2;
Nt = 50000;
t0 = 0e-6;
dt = 1/fs;
BW = pc_param.f1-pc_param.f0;
alpha = BW/pc_param.Tpd;
pc_param.time = t0 + (0:dt:(Nt-1)*dt).';
signal = cos(2*pi*pc_param.f0*pc_param.time + pi*alpha*pc_param.time.^2);
% Apply Tukey window
good_mask = pc_param.time>=0 & pc_param.time<=pc_param.Tpd;
signal(~good_mask) = 0;
signal(good_mask) = signal(good_mask) .* tukeywin_trim(sum(good_mask),pc_param.tukey);

[pc_signal,pc_time] = pulse_compress(signal,pc_param);

signalIQ = exp(1i*(2*pi*pc_param.f0*pc_param.time + pi*alpha*pc_param.time.^2));
% Apply Tukey window
signalIQ(~good_mask) = 0;
signalIQ(good_mask) = signalIQ(good_mask) .* tukeywin_trim(sum(good_mask),pc_param.tukey);

[pc_signalIQ,pc_time] = pulse_compress(signalIQ,pc_param);

plot(pc_time*1e6, lp(pc_signal,2));
hold on;
plot(pc_param.time*1e6, lp(signal,2),'r');
plot(pc_time*1e6, lp(pc_signalIQ,2),'g');
hold off;
xlabel('Time (us)');
ylabel('Relative power (dB)');

%% Example #4
% MCRDS: showing use of zero padding and variable bandwidth fields to deal
% with small time bandwidth products
clear pc_param;
pc_param.f0 = 140e6;
pc_param.f1 = 160e6;
pc_param.Tpd = 2.5e-6;
pc_param.BW = [25e6 25e6];
pc_param.zero_pad = -4e-6;
pc_param.baseband = true;
pc_param.window_func = inline('tukeywin_trim(N,0.1)');
pc_param.Mt = 10;
fs = 120e6;
Nt = 20000;
t0 = 0e-6;
dt = 1/fs;
BW = pc_param.f1-pc_param.f0;
alpha = BW/pc_param.Tpd;
pc_param.time = t0 + (0:dt:(Nt-1)*dt).';
signal = cos(2*pi*pc_param.f0*pc_param.time + pi*alpha*pc_param.time.^2);
signal(pc_param.time<0 | pc_param.time>pc_param.Tpd) = 0;

pc_param.tukey = 0;
[pc_signal,pc_time] = pulse_compress(signal,pc_param);

figure(1); clf;
plot(pc_time*1e6, lp(pc_signal,2));
hold on;
plot(pc_param.time*1e6, lp(signal,2),'r');
hold off;
xlabel('Time (us)');
ylabel('Relative power (dB)');

%% Example #4
% MCRDS: showing use of zero padding and variable bandwidth fields to deal
% with small time bandwidth products
clear pc_param;
pc_param.f0 = 150e6;
pc_param.f1 = 152e6;
pc_param.Tpd = 3e-6;
pc_param.BW = [2e6 2e6];
pc_param.zero_pad = -3e-6;
pc_param.baseband = true;
pc_param.window_func = inline('tukeywin_trim(N,0.1)');
pc_param.Mt = 10;
fs = 120e6;
Nt = 20000;
t0 = 0e-6;
dt = 1/fs;
BW = pc_param.f1-pc_param.f0;
alpha = BW/pc_param.Tpd;
pc_param.time = t0 + (0:dt:(Nt-1)*dt).';
signal = cos(2*pi*pc_param.f0*pc_param.time + pi*alpha*pc_param.time.^2);
signal(pc_param.time<0 | pc_param.time>pc_param.Tpd) = 0;

pc_param.tukey = 0;
[pc_signal,pc_time] = pulse_compress(signal,pc_param);

figure(1); clf;
plot(pc_time*1e6, lp(pc_signal,2));
hold on;
plot(pc_param.time*1e6, lp(signal,2),'r');
hold off;
xlabel('Time (us)');
ylabel('Relative power (dB)');

hold on;
pc_param.BW = [25e6 25e6];
[pc_signal,pc_time] = pulse_compress(signal,pc_param);
plot(pc_time*1e6, lp(pc_signal,2),'g');

