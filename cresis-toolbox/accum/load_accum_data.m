function [data,time] = load_accum_data(param,records)
% [data,time] = load_accum_data(param,records)
%
% param = struct controlling the loading and processing
% records = records struct with "gps_time" GPS time vector and "elev"
%   WGS-84 elevation vector
%
% data = pulse compressed output
% time = fast-time axis
%
% Author: John Paden

% =======================================================================
% Check Arguments

ENABLE_DECONVOLUTION_WF_CREATION = false;
if ENABLE_DECONVOLUTION_WF_CREATION
  %% Currently, creation of deconvolution waveforms is done by hand
  param.proc.elev_correction = true;
end

if ~isfield(param.proc,'elev_correction') || isempty(param.proc.elev_correction)
  param.proc.elev_correction = false;
end

%% Load data
% =======================================================================

rec = param.load.recs(1);
rec_load_offset = 0;

wf = 1;
total_recs = param.load.recs(2)-param.load.recs(1)+1;
rline = 0;
data = zeros(1,total_recs,16,'single');
while rec <= param.load.recs(end)
  % Determine which file contains this record
  file_idx = find(param.load.file_rec_offset{1} > rec,1) - 1;
  [tmp fn_name fn_ext] = fileparts(param.load.filenames{1}{file_idx});
  fn = fullfile(param.load.filepath, [fn_name fn_ext]);
  fprintf('  Loading %s\n', fn);
  
  % Determine offset to first record to load from this file
  start_rec = rec - param.load.file_rec_offset{1}(file_idx);
  
  % Determine how many records we can load from this file
  num_rec = min(param.load.file_rec_offset{1}(file_idx+1) - rec, ...
    param.load.recs(end) - rec + 1);
  
  if param.load.file_version == 101
    %% File Version 101: John Ledford and Carl Leuschen 1U-DAQ system
    % Open file
    [fid,msg] = fopen(fn,'r','ieee-be');
    if fid < 1
      fprintf('Could not open file %s\n', fn);
      error(msg);
    end
    
    % Read in records
    HEADER_SIZE = 8*4; % in bytes
    WF_HEADER_SIZE = 16*8; % in bytes
    rec_idx = 1;
    fseek(fid,param.load.offset(rec_load_offset+rec_idx) + HEADER_SIZE,-1);
    num_sam = mod(fread(fid,1,'uint32'), 2^14) - 16 - 4;
    tmp = fread(fid,1,'uint32');
    bit_shifts = mod(floor(tmp/2^24), 2^4);
    t0 = mod(floor(tmp/2^10), 2^14) / param.radar.fs - 8.3450e-06;
    presums = 1 + mod(tmp, 2^10);
    
    % Read in records
    for rec_idx = 1:num_rec
      rline = rline + 1;
      fseek(fid,param.load.offset(rec_load_offset+rec_idx) + HEADER_SIZE + WF_HEADER_SIZE,-1);
      for wf = 1:16
        fseek(fid,32,0);
        data(1:num_sam,rline,wf) = fread(fid,num_sam,'uint16=>single');
        fseek(fid,8,0);
      end
    end
    % Close file
    fclose(fid);
    
    % rec: This is the absolute record number in the segment
    rec = rec + num_rec;
    % rec_load_offset: This is the relative record number of this particular call to
    % load_fmcw_data
    rec_load_offset = rec_load_offset + num_rec;
    
  else
    error('Unsupported file version');
    
  end
end

%% Remove digital error bursts
% =======================================================================
%  Probably should make these values programmable from spreadsheet, not
%  sure if they are always 44047 and 3840
test_idxs = find(data(:) == 44047).';
for test_idx = test_idxs
  if test_idx >= 1 && data(test_idx-1) < data(test_idx) ...
      && test_idx < numel(data) && data(test_idx+1) < data(test_idx)
    data(max(1,test_idx-2):min(test_idx+1,numel(data))) ...
      = 2^(param.radar.adc_bits-1)*presums/2^bit_shifts;
  end
end
test_idxs = find(data(:) == 3840).';
for test_idx = test_idxs
  if test_idx >= 1 && data(test_idx-1) < data(test_idx) ...
      && test_idx < numel(data) && data(test_idx+1) < data(test_idx)
    data(max(1,test_idx-2):min(test_idx+1,numel(data))) ...
      = 2^(param.radar.adc_bits-1)*presums/2^bit_shifts;
  end
end

%% Shift data from being 1's complement to 2's complement
% =======================================================================
data = data - 2^(param.radar.adc_bits-1)*presums/2^bit_shifts;

%% Presum
% =======================================================================
if param.proc.presums > 1
  % Elevation correction just corrects for deviations from the mean elevation
  % within the presum window.
  drange = fir_dec(records.elev,param.proc.presums);
  drange = reshape(repmat(drange,[param.proc.presums 1]),[1 param.proc.presums*length(drange)]);
  drange = records.elev(1:length(drange)) - drange;
  if ENABLE_DECONVOLUTION_WF_CREATION
    drange = records.elev(1:length(drange)) - mean(records.elev);
  end
  num_rec = length(drange)/param.proc.presums;
  data = data(:,1:length(drange),:);
  Nt = size(data,1);
  df = param.radar.fs/Nt;
  physical_constants;
  param.proc.elev_correction = 0;
  
  for wf = 1:size(data,3)
    if param.proc.elev_correction
      freq_relative = df*ifftshift(-floor(Nt/2) : floor((Nt-1)/2)).';
      freq = -(param.radar.wfs.step.fLO + (wf-1)*param.radar.wfs.step.f_step) + freq_relative;
      freq(freq_relative<0) = freq(freq_relative<0) - 2*freq(1);

      % This sign of the exp(-j ...) correction seems backwards, yet works...
      % Greater range implies more negative phase, so the correction should be
      % a more positive phase and this is the opposite. Not sure why??? May have
      % something to do with how we do complex baseband with FFT... maybe this
      % is causing a conjugation of the phase.
      data(:,:,wf) = ifft(fft(data(:,:,wf)) .* exp(j*2*pi*freq*drange/(c/2)));
      
    end
    data(:,1:num_rec,wf) = fir_dec(data(:,:,wf),param.proc.presums);
  end
  records.gps_time(:,1:num_rec) = fir_dec(records.gps_time,param.proc.presums);
  records.elev(:,1:num_rec) = fir_dec(records.elev,param.proc.presums);
  records.gps_time = records.gps_time(:,1:num_rec);
  records.elev = records.elev(:,1:num_rec);
  data = data(:,1:num_rec,:);
end

%% Conversion from quantization to voltage
% =======================================================================
adc_bits = param.radar.adc_bits;
Vpp_scale = param.radar.Vpp_scale;
data = data * Vpp_scale/2^adc_bits * 2^bit_shifts / presums;

%% Pulse compress parameters
% =======================================================================
pc_param.f0 = abs(param.radar.wfs.step.f0 - param.radar.wfs.step.fLO);
pc_param.f1 = abs(param.radar.wfs.step.f1 - param.radar.wfs.step.fLO);
pc_param.Tpd = param.radar.wfs.Tpd;
Nt_ref = floor(pc_param.Tpd*param.radar.fs) + 1;
dt = 1/param.radar.fs;
pc_param.tukey = param.radar.wfs.tukey;
pc_param.zero_pad = false;
pc_param.td_window_func = '';
pc_param.window_func = param.radar.wfs.step.ft_window;
pc_param.decimate = false;

if 0
  % Add a time domain window before performing FFT: Does not seem to help.
  if ~isempty(param.radar.wfs.step.start_time_window)
    param.radar.wfs.step.start_time_window = param.radar.wfs.step.start_time_window(1:end/2);
    % Window first part
    data(1:length(param.radar.wfs.step.start_time_window)/2,:,:) ...
      = data(1:length(param.radar.wfs.step.start_time_window)/2,:,:) ...
      .* repmat(param.radar.wfs.step.start_time_window(1:end/2),[1 size(data,2) size(data,3)]);
    % Window last part
    data(end-length(param.radar.wfs.step.start_time_window)/2+1:end,:,:) ...
      = data(end-length(param.radar.wfs.step.start_time_window)/2+1:end,:,:) ...
      .* repmat(param.radar.wfs.step.start_time_window(end/2+1:end),[1 size(data,2) size(data,3)]);
  end
end

% Adjust zero-pad so that frequency bins from each step chirp overlap
Nt_pc = Nt_ref-1 + size(data,1);
step_integer = ceil(Nt_pc*abs(param.radar.wfs.step.f_step)/param.radar.fs);
step_integer_mult = lcm(param.radar.fs,abs(param.radar.wfs.step.f_step))/param.radar.fs;
step_integer = step_integer + mod(step_integer_mult-mod(step_integer,step_integer_mult),step_integer_mult);
Nt_pc_aligned = step_integer*param.radar.fs/abs(param.radar.wfs.step.f_step);
Nt_ref = Nt_ref + (Nt_pc_aligned - Nt_pc);
pc_param.time = t0 - param.radar.wfs.Tsys - (Nt_ref-1)*dt + dt*(0:(size(data,1)+Nt_ref-1)-1).';

% Zero-pad for pulse compression
data = cat(1,zeros(Nt_ref-1,size(data,2),size(data,3)),data);

clear pc_data;
for wf_idx = 1:size(data,3)
  [pc_data(:,:,wf_idx),pc_time] ...
    = pulse_compress(data(:,:,wf_idx),pc_param);
end
clear data;

%% Combine Waveforms
% =======================================================================

alpha = (param.radar.wfs.step.f1 - param.radar.wfs.step.f0)/param.radar.wfs.Tpd;
comb_f0 = param.radar.wfs.step.f0 + alpha*param.radar.wfs.Tpd*0.5 + abs(param.radar.wfs.step.f_step/2) + param.radar.wfs.step.f_offset;
comb_f1 = comb_f0 + param.radar.wfs.step.f_step * size(pc_data,3) + param.radar.wfs.step.f_offset;

Nt = size(pc_data,1);
df = param.radar.fs/Nt;

Nt_step = abs(param.radar.wfs.step.f_step) / df;

min_f0 = min(comb_f0,comb_f1);
comb_Nt = 16*Nt_step;
comb_freq = min_f0 + df*(0:comb_Nt-1);

pc_freq = param.radar.wfs.step.fLO + -df*(0:Nt-1);
[~,pc_data_idx] = min(abs(pc_freq - comb_f0));

pc_data = fft(pc_data);
comb_freq = ifftshift(comb_freq,1);

% Determine pulse compression phase offset
if 1
  pc_phase_offset(1) = 0;
  for wf=2:16
    pc_phase_offset(wf) = angle(mean(mean(conj(pc_data(pc_data_idx+Nt_step+(-10:10),:,wf-1)) .* conj(conj(pc_data(pc_data_idx+(-10:10),:,wf))),2)));
  end
  pc_phase_offset = cumsum(pc_phase_offset);
  pc_phase_offset
end

data = zeros(length(comb_freq),size(pc_data,2));
for wf = 1:16
  % Combine chirps
  % 1. Flip spectrum and conjugation of phase because negative side of frequencies were used
  % 2. Pulse compression phase shift because the phase of each of the pulse compression filters
  %    won't align without this.
  data(size(data,1)-(wf-1)*Nt_step + (0:-1:-Nt_step+1),:) = exp(j*pc_phase_offset(wf))*conj(pc_data(pc_data_idx+(0:Nt_step-1),:,wf));
end
data = ifftshift(data,1);

%% Deconvolution Creation
% =======================================================================
if ENABLE_DECONVOLUTION_WF_CREATION
  %% Debug code for creating deconvolution file: Run with the specific frame and block that you need
  % 20110325_01_015, block 10 in get heights, deconv_rlines = 480:540;
  %deconv_rlines = 480:540;
  %deconv_wind_len = 1000;
  % 20100525_02_022, block 4 in get heights, deconv_rlines = 201:216;
  %deconv_rlines = 201:216;
  %deconv_wind_len = 300;
  % 20100525_02_024, block 4 in get heights, deconv_rlines = 16:27;
  deconv_rlines = [1     2     3     6    10    11    13    14    19    20    21    22    23    24    25    26    31    36    39    40];
  deconv_wind_len = 300;
  % 20100507_01_049, block 11 in get heights, deconv_rlines = 138:149;
%   deconv_rlines = 138:149;
%   deconv_wind_len = 300;
  
  [max_val,max_idx] = max(interpft(ifft(data(:,deconv_rlines)),20*size(data,1)));
  max_idx = max_idx/20
  max_idx = max_idx - round(mean(max_idx));
  Nt = size(data,1);
  freq_norm = ifftshift(-floor(Nt/2) : floor((Nt-1)/2)).';
  
  data(:,deconv_rlines) = data(:,deconv_rlines) .* exp(j*2*pi*freq_norm*max_idx/Nt);

  [max_val,max_idx] = max(interpft(ifft(data(:,deconv_rlines)),20*size(data,1)));
  max_idx = max_idx/20
  figure(1); clf;
  plot(abs(max_idx))
  
  % 
  [max_val,max_idx] = max(ifft(data(:,deconv_rlines)));
  figure(1); clf;
  plot(abs(max_idx))
  
  figure(3); clf;
  plot(abs(max_val))
  angle_fix = unwrap(angle(max_val));
  plot(angle_fix);
  
  deconv_filt = mean(data(:,deconv_rlines) .* repmat(1./max_val,[size(data,1) 1]) ,2);
  
  deconv_filt = fft(circshift(ifft(deconv_filt),round(mean(-max_idx))));
  tukeywin_pad = ifftshift(tukeywin(deconv_wind_len));
  tukeywin_pad = [tukeywin_pad(1:deconv_wind_len/2); zeros(length(deconv_filt)-deconv_wind_len,1); ...
    tukeywin_pad(deconv_wind_len/2+1:deconv_wind_len)];
  
  %deconv_filt = fft(ifft(deconv_filt) .* tukeywin_pad);
  
  fn_deconv = fullfile(ct_filename_out(param,'','',1),param.radar.wfs.step.fn_deconv);
  
  %param.proc.ft_wind = inline('boxcar(N)'); % Override window for testing
  deconv_filt_windowed = ifftshift(param.proc.ft_wind(length(deconv_filt))) ./ deconv_filt;
  deconv_data = data .* repmat(deconv_filt_windowed,[1 size(data,2)]);
  deconv_data = ifft(deconv_data);
  
  colormap(1-gray(256))
  test = lp(fir_dec(abs(deconv_data).^2,ones(1,11)/11));
  test2 = lp(fir_dec(abs(ifft(data)).^2,ones(1,11)/11));
  figure(3); clf;
  clf;
  imagesc(test);
  colormap(1-gray(256));
  %ylim([1440 1540])

  % Create before/after deconvolution results
  deconv_result = mean(deconv_data(:,deconv_rlines),2);
  deconv_orig = mean(abs(ifft(data(:,deconv_rlines))).^2,2);
  
  % Normalize the deconvolution filter
  deconv_filt = deconv_filt * sqrt(max(abs(deconv_result).^2) / max(deconv_orig));
  
  fprintf('Saving deconvolution results %s\n', fn_deconv);
  save(fn_deconv,'deconv_filt','param','deconv_rlines','pc_phase_offset','deconv_result','deconv_orig');
  
end

%% Deconvolution
% =======================================================================
clear pc_data; % No longer needed

if ~isempty(param.radar.wfs.step.fn_deconv)
  fn_deconv = fullfile(ct_filename_out(param,'','',1),param.radar.wfs.step.fn_deconv);
  deconv = load(fn_deconv);
  deconv.deconv_filt = ifftshift(param.proc.ft_wind(length(deconv.deconv_filt))) ./ deconv.deconv_filt;
  if param.radar.wfs.step.f_offset ~= deconv.param.radar.wfs.step.f_offset
    error('Pulse compression frequency offsets must match otherwise deconvolution will fail.')
  end
  if ~isequal(param.radar.wfs.step.ft_window,deconv.param.radar.wfs.step.ft_window)
    error('Stepped frequency pulse compression window must be the same for deconvolution.')
  end

  data = data .* repmat(deconv.deconv_filt,[1 size(data,2)]);
end

data = ifft(data);

% Create time axis
BW = abs(param.radar.wfs.step.f_step)*16;
dt = 1/BW;
Nt = size(data,1);
time = pc_time(1) + dt*(0:Nt-1);

if 0
  %% Debug Test Code
  figure(1); clf;
  colormap(1-gray(256))
  test = lp(fir_dec(abs(data).^2,ones(1,11)/11));
  clf;
  imagesc(test);
  ylim([1180 1260])
end

end
