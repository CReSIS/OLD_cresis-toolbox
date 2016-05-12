% script run_basic_file_loader
%
% Example script for running basic_file_loader.m and using its outputs.
%
% Author: John Paden

%% Load data
if 1
  [param,defaults] = default_radar_params_2016_Greenland_Polar6_mcords;
  
  [data,fn,settings,default,gps,hdr,pc_param] = basic_file_loader(param,defaults);
  param.img = pc_param.img;
end

%% Check Saturation
if 1
  wf = abs(param.img(1,1));
  adc = abs(param.img(1,2));
  figure(1); clf;
  Nt = size(data,1);
  dt = pc_param.time(2)-pc_param.time(1);
  time = pc_param.time;
  SATURATION_CORRECTION_FACTOR = 1722.5/4096;
  h_plot = plot([1 size(data,2)],+SATURATION_CORRECTION_FACTOR*(2^12/2)*hdr.wfs(wf).presums/2^hdr.wfs(wf).bit_shifts*[1 1],'LineWidth',4,'LineStyle','--');
  color = get(h_plot,'color');
  hold on;
  plot([1 size(data,2)],-SATURATION_CORRECTION_FACTOR*(2^12/2-1)*hdr.wfs(wf).presums/2^hdr.wfs(wf).bit_shifts*[1 1],'color',color,'LineWidth',4,'LineStyle','--');
  h_plot = plot(max(real(data)),'color',color);
  plot(max(imag(data)),'color',color);
  plot(min(real(data)),'color',color);
  plot(min(imag(data)),'color',color);
  title(sprintf('wf: %d, adc: %d', wf, adc));
end

%% Remove DC
if 1
  for wf_adc = 1:size(param.img,1)
    wf = abs(param.img(wf_adc,1));
    if isfield(default.radar.wfs(wf),'DC_adjust') && ~isempty(default.radar.wfs(wf).DC_adjust)
      adc = abs(param.img(wf_adc,2));
      tmp = load(fullfile(ct_filename_out(param,'noise','',1), ...
        default.radar.wfs(wf).DC_adjust),'DC_adjust');
      data(:,:,wf_adc) = data(:,:,wf_adc) - tmp.DC_adjust(adc);
    end
  end
end

%% Convert from quantization to voltage @ ADC
if 1
  wf = abs(param.img(1,1));
  data = data ...
    * default.radar.adc_full_scale/2^default.radar.adc_bits ...
    * 2^hdr.wfs(abs(wf)).bit_shifts / hdr.wfs(wf).presums;
end

%% Additional software presums
if 1
  param.presums = 10;
  for wf_adc = 1:size(data,3)
    data(:,1:floor(size(data,2)/param.presums),wf_adc) = fir_dec(data(:,:,wf_adc),param.presums);
  end
  data = data(:,1:floor(size(data,2)/param.presums),:);
  hdr.radar_time = fir_dec(hdr.radar_time,param.presums);
  hdr.gps_time = fir_dec(hdr.gps_time,param.presums);
  hdr.lat = fir_dec(hdr.lat,param.presums);
  hdr.lon = fir_dec(hdr.lon,param.presums);
  hdr.elev = fir_dec(hdr.elev,param.presums);
  hdr.roll = fir_dec(hdr.roll,param.presums);
  hdr.pitch = fir_dec(hdr.pitch,param.presums);
  hdr.heading = fir_dec(hdr.heading,param.presums);
end

%% Pulse compression
if 1
  [pc_signal,pc_time] = pulse_compress(data,pc_param);
end

%% Deconvolution
if 1
  wf = param.img(1,1);
  adc = param.img(1,2);
  ref_fn_name = char(default.radar.ref_fn);
  ref_fn_name = regexprep(ref_fn_name,'%w',sprintf('%.0f',wf));
  ref_fn_name = regexprep(ref_fn_name,'%a',sprintf('%.0f',adc));
  ref_fn = fullfile(ct_filename_out(param,'noise','',1), [ref_fn_name '.mat']);
  
  deconv_param = pc_param;
  deconv_param.ref_fn = ref_fn;
  [deconv_signal,deconv_time] = pulse_compress(data,deconv_param);
end

return;

%% Motion compensation
% ======================================================================
param.averaging_fh = @mean;
param.time = pc_time;
dt = pc_time(2) - pc_time(1);
Nt = length(pc_time);
df = 1/(Nt*dt);
param.freq = pc_param.DDC_freq + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';

if all(gps.roll==0)
  param.mocomp_type = 2;
else
  param.mocomp_type = 4;
end
param.tx_weights = double(settings.DDS_Setup.Ram_Amplitude);
param.rx_paths = {}; param.rx_paths{wf} = default.radar.rx_paths;
param.lever_arm_fh = @lever_arm;

param.combine_channels = false;
param.snr_threshold = 10;
param.phase = default.radar.wfs(1).chan_equal_deg;
param.amp = default.radar.wfs(1).chan_equal_dB;
param.td = default.radar.wfs(1).chan_equal_Tsys;

freq = param.freq;
data = pc_signal;

mocomp_param.type = param.mocomp_type;
mocomp_param.tx_weights = param.tx_weights;
mocomp_param.season_name = param.season_name;
mocomp_param.radar_name = param.radar_name;
mocomp_param.gps_source = hdr.gps_source;
for wf_adc_idx = 1:size(data,3)
  wf = abs(param.img(wf_adc_idx,1));
  adc = param.img(wf_adc_idx,2);
  mocomp_param.rx = param.rx_paths{wf}(adc);
  
  % drange = change in range (positive is longer range)
  drange = basic_motion_comp(mocomp_param,param.lever_arm_fh,hdr.roll, ...
    hdr.pitch,hdr.heading,hdr.lat,hdr.lon,hdr.elev);
  % dtime = Convert to time (in air), positive is longer range/time-delay
  dtime = 2*drange/3e8;
  Nt = size(data,1);
  Nx = size(data,2);
  if 0
    figure(1); clf;
    plot(dtime*1e9);
    title('adc %d');
    pause;
  end
  % Time shift data in the frequency domain
  %   A positive dtime implies a larger negative phase delay (since
  %   delay is negative/lagging phase)
  data(:,:,wf_adc_idx) = ifft(fft(data(:,:,wf_adc_idx)) ...
    .*exp(-1i*2*pi*repmat(freq,1,Nx).*repmat(dtime,Nt,1)));
end

%% Apply time correction
%   Time delay is removed (positive moves targets closer in range)
% =======================================================================
w = 2*pi*freq;
param.averaging_fh = @mean;
param.time = pc_time;
dt = pc_time(2) - pc_time(1);
Nt = length(pc_time);
df = 1/(Nt*dt);
param.freq = pc_param.DDC_freq + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';

param.rx_paths = {}; param.rx_paths{wf} = default.radar.rx_paths;
param.phase = default.radar.wfs(1).chan_equal_deg;
param.amp = default.radar.wfs(1).chan_equal_dB;
param.td = default.radar.wfs(1).chan_equal_Tsys;

for wf_adc_idx = 1:size(data,3)
  wf = abs(param.img(wf_adc_idx,1));
  adc = param.img(wf_adc_idx,2);
  data(:,:,wf_adc_idx) = ifft(fft(data(:,:,wf_adc_idx)).*exp(1i*repmat(w,1,size(data,2))*param.td(param.rx_paths{wf}(adc))));
end

%% Apply amplitude and phase correction
%   Amp/phase are DIVIDED out as opposed to being multiplied
% =======================================================================

for wf_adc_idx = 1:size(data,3)
  wf = abs(param.img(wf_adc_idx,1));
  adc = param.img(wf_adc_idx,2);
  data(:,:,wf_adc_idx) = data(:,:,wf_adc_idx) ./ (10^(param.amp(param.rx_paths{wf}(adc))/20).*exp(1i*param.phase(param.rx_paths{wf}(adc))/180*pi));
end

%% Apply elevation correction
physical_constants;

wf_adc = 1;
elev_filt_len = round(length(hdr.elev(wf_adc,:))/20)*2+1;
hdr.elev(wf_adc,:) = sgolayfilt(hdr.elev(wf_adc,:), 2, elev_filt_len, hanning(elev_filt_len));

wf = abs(param.img(wf_adc,1));
adc = abs(param.img(wf_adc,2));

% Create frequency axis
dt = pc_time(2) - pc_time(1);
Nt = length(pc_time);
T = Nt*dt;
df = 1/T;
freq = pc_param.DDC_freq + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';

% Correct all the data to a constant elevation (no zero padding is
% applied so wrap around could be an issue)
for rline = 1:size(data,2)
  elev_dt = (hdr.elev(rline) - hdr.elev(1)) / (c/2);
  data(:,rline,wf_adc) = ifft(fft(data(:,rline,wf_adc)) .* exp(1i*2*pi*freq*elev_dt));
end

%% Track surface
param.ref_wf_adc = 1;
ml_data = lp(fir_dec(abs(data(:,:,param.ref_wf_adc)).^2,ones(1,5)/5,1));
good_time_bins = find(pc_time > pc_param.Tpd*1.1 & pc_time > default.basic_surf_track_min_time);
[max_value,surf_bin] = max(ml_data(good_time_bins,:));
surf_bin = surf_bin + good_time_bins(1)-1;


%% Track time and phase delay of surface
best_idx = 1150;
rlines = 1100:1200;

max_value = zeros(1,size(data,2));
max_idx_unfilt = zeros(1,size(data,2));
Mt = 100;
for idx = 1:size(data,2)
  oversampled_rline = interpft(data(:,idx),size(data,1)*Mt);
  start_bin = surf_bin(idx)*Mt-1000;
  [max_value(idx),max_idx_unfilt(idx)] ...
    = max(oversampled_rline(start_bin:end));
  max_idx_unfilt(idx) = max_idx_unfilt(idx) + start_bin - 1;
end

% Filter the delay (max_idx) and phase (max_value) vectors
max_idx = sgolayfilt((max_idx_unfilt-1)/Mt+1,3,17);
max_idx = max_idx - mean(max_idx(rlines));
phase_corr = sgolayfilt(double(unwrap(angle(max_value))),3,17);

%% Compensate range lines for phase and delay variance
% in the peak value
T = Nt*dt;
df = 1/T;
fc = (pc_param.f0+pc_param.f1)/2;
freq = pc_param.DDC_freq + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';

% Apply true time delay shift to flatten surface
comp_data = ifft(fft(data(:,:,wf_adc)) .* exp(1i*2*pi*freq*max_idx/Mt*dt) );
% Apply phase correction (compensating for phase from time delay shift)
comp_data = comp_data .* repmat(exp(-1i*(phase_corr + 2*pi*fc*max_idx/Mt*dt)), [Nt 1]);

