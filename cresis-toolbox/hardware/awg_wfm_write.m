% script awg_wfm_write
%
% Creates WFM files for the Tektronix 25 GSPS AWG70002
% Full scale is +/- 1 floating point scale (32 bit float)
% Full scale output voltage is +/-0.25 V  (so +/- 1.0 in the file produces 0.5 Vpp)
%
% Edit "User Settings" section to create the signal waveform that you want.
%
% Todo: Ability to set markers individually not determined.
%
% Authors: Mahmood Hameed, John Paden

% =========================================================================
%% User Settings
% =========================================================================

fn_dir = 'd:\tmp\paden_20151023\';

% Create signals (enable/disable creation by setting "if 0" and "if 1")
signal_name = {};
signal = {};
marker = {};

if 1
  fs = 25e9;
  dt = 1/fs;
  Tpd = 240e-6*6/16;
  Nt = (Tpd + 16e-6)*fs;
  time = dt * (0:Nt-1).';
  f0 = 2e9;
  f1 = 8e9;
  BW = (f1-f0);
  K = BW/Tpd;
  td = 0e-6;
  
  signal_name{end+1} = 'chirp2to8';
  
  signal{end+1} = tukeywin_cont((time-Tpd/2-td)/Tpd, 0) .* cos(2*pi*f0*(time-td) + pi*K*(time-td).^2);

  % Marker: creates 1 us pulse at start on both marker channels
  marker{end+1} = zeros(size(signal{end}),'uint8');
  marker{end}(time<1e-6,:) = 255;
end

if 1
  fs = 25e9;
  dt = 1/fs;
  Tpd = 240e-6*6/16;
  Nt = (Tpd + 16e-6)*fs;
  time = dt * (0:Nt-1).';
  f0 = 2e9;
  f1 = 8e9;
  BW = (f1-f0);
  K = BW/Tpd;
  td = 3e-6;
  
  signal_name{end+1} = 'chirp2to8_delayed2us';
  
  signal{end+1} = tukeywin_cont((time-Tpd/2-td)/Tpd, 0) .* cos(2*pi*f0*(time-td) + pi*K*(time-td).^2);

  % Marker: creates 1 us pulse at start on both marker channels
  marker{end+1} = zeros(size(signal{end}),'uint8');
  marker{end}(time<1e-6,:) = 255;
end

if 1
  fs = 25e9;
  dt = 1/fs;
  Tpd = 240e-6*6/16;
  Nt = (Tpd + 16e-6)*fs;
  time = dt * (0:Nt-1).';
  f0 = 2e9;
  f1 = 2e9;
  BW = (f1-f0);
  K = BW/Tpd;
  td = 0e-6;
  
  signal_name{end+1} = 'tone2GHz';
  
  signal{end+1} = tukeywin_cont((time-Tpd/2-td)/Tpd, 0) .* cos(2*pi*f0*(time-td) + pi*K*(time-td).^2);

  % Marker: creates 1 us pulse at start on both marker channels
  marker{end+1} = zeros(size(signal{end}),'uint8');
  marker{end}(time<1e-6,:) = 255;
end

% =========================================================================
%% Automated Section
% =========================================================================

crlf = hex2dec({'0d' '0a'});
for sig_idx = 1:length(signal)
  fn = fullfile(fn_dir, sprintf('%s.WFM', signal_name{sig_idx}));
  samples = size(signal{sig_idx},1);
  
  fprintf('Writing %d samples to %s\n', samples, fn);
  
  fid = fopen(fn, 'w+');
  
  fprintf('  Writing header\n');
  fwrite(fid, 'MAGIC 1000');
  fwrite(fid, crlf);

  sample_bytes = num2str(samples*5);
  length_bytes = num2str(length(sample_bytes));
  header = ['#' length_bytes sample_bytes];
  fwrite(fid, header);
  
  fprintf('  Writing waveforms\n');
  % Convert signal to 32 bit float
  start_pos = ftell(fid);
  fwrite(fid, single(signal{sig_idx}(1)), 'float32');
  fwrite(fid, single(signal{sig_idx}(2:end)), 'float32', 1); % skip marker bytes
  
  fprintf('  Writing markers\n');
  fseek(fid,start_pos,-1);
  fwrite(fid, marker{sig_idx}(1:samples), 'uint8', 4); % skip waveform bytes
  
  fprintf('  Writing footer\n');
  fwrite(fid, 'CLOCK 25e+09');
  fwrite(fid, crlf);
  
  fclose(fid);
  fprintf('  Done\n');
end

return;
