% script create_accum_DDS_amp_correction
%
% Two Options:
% 1. Can create an undistorted codex (RAM) file for the 1U-DAQ DDS
%    Set create_base_codex_file_flag to true.
% 2. Analyze measured results and optionally create a predistorted codex
% (RAM file) for the 1U-DAQ DDS.
%    Set create_new_codex_file_flag to create optional codex file.
%
% Specifically designed for the accum2 digital system (20120219). This
% function does predistortion for the DDS output and uses 900 MHz
% down to 600 MHz sampling alias of the 1 GHz sampled DDS.
%
% Steps:
% 1. Run "Creation of BASE/TESTING codex.dat file for magnitude
%    measurement"
% 2. Take data on the digital storage oscilloscope (DSO) which properly
%    captures the DDS/accumulation-radar output
%   a. 5 GSPS, 100k points, 128 averages, CSV file
%   b. DO NOT CHANGE THESE SETTINGS OR THE TRIGGER HOLD OFF... KEEP THEM
%      FOR ALL YOUR TESTING
% 3. Set dso_fn = 'tek0007.csv'; to the actual DSO filename
% 4. Run script
% 5. Use figure 2 to select truncation_bins of DSO output (avoid
%    including transients at beginning/end that will be hard to match/remove)
%    Adjust beginning of time scale with fudge_factor.
% 6. After figure 3 looks good, run "Creation of Corrected codex.dat file"
% 7. Recapture DSO data with the same settings and trigger holdoff
%    (otherwise you'll be stuck retruncating the data)
% 8. Set dso_fn to this new file, run script, and then check figure 3
%    to see how flat the output was (don't bother looking at the
%    fitted/red-line... its not important that it looks good).
%
% Author: Cameron Lewis, John Paden

% ====================================================================
% User Settings
% ====================================================================

% create_base_codex_file_flag: boolean, set to true if you want the base
%   (undistorted) codex file to be generated. Nothing else runs if this
%   is set to true
create_base_codex_file_flag = true;

% fc = scalar center frequency (Hz)
fc = 750e6;

% Tpd = pulse duration of chirp (sec)
Tpd = 1.024e-6;

CODEX_LENGTH = 256*1; % 256 for 1.024e-6, 512 for 2.048e-6 Tpd
CODEX_ZERO_BINS = 20;
CODEX_MAX_VALUE = 15341;
TUKEY_WIN = 0.1;

% create_new_codex_file_flag: boolean, set to true if you want a
%   predistorted codex file to be generated, requires dso_fn be set to the
%   measurement file
create_new_codex_file_flag = true;

% dso_fn: string containing the Tektronix DSO measurement file (you must
%   save time/amplitude in the CSV file)
dso_fn = 'D:\tmp\accum_NA\codex\accum_DSO_capture_20150215.mat';

out_fn_dir = 'D:\tmp\accum_NA\codex\';

% truncation_bins: which bins in the DSO output (after downsampling) are
%   the bins containing the chirp output. Use figure 2 to figure out
%   which bins are the correct bins to truncate.
truncation_bins = [];

% NOTE: THIS IS PROTOTYPE STAGE CODE, YOU WILL PROBABLY HAVE TO MODIFY
% SOMETHING BELOW THIS POINT...

% ====================================================================
% Automated Section
% ====================================================================

if create_base_codex_file_flag
  %% Creation of BASE/TESTING codex.dat file for magnitude measurement
  
  % Start with flat output (staying away from full scale here to avoid
  % risk of saturation)
  codex = ones(CODEX_LENGTH,1)*CODEX_MAX_VALUE;
  % Set first and last sample to zero
  % Also: Due to DDS programming method, first 7-8 samples must be 0, in this
  % case we do first 8 samples
  codex([1:CODEX_ZERO_BINS-1, end]) = 0;
  %codex(codex~=0) = round(max(codex)*tukeywin_trim(sum(codex~= 0),TUKEY_WIN));
  % DDS reads it in backwards
  codex = flipud(codex);
  
  out_fn = fullfile(out_fn_dir,sprintf('codex_flat%.0f.dat',max(codex)));
  fprintf('Writing %s\n', out_fn)
  dlmwrite(out_fn,codex,'newline','pc','precision','%.0f');
  %     fid = fopen('codex35.dat','w');
  %     fprintf(fid,'%.0f\r\n',codex);
  %     fclose(fid);
  
  return;
end

[tmp tmp dso_fn_ext] = fileparts(dso_fn);
if strcmpi(dso_fn_ext,'.csv')
  [time,data] = textread(dso_fn,'%f%f','delimiter',',','headerlines',15);
elseif strcmpi(dso_fn_ext,'.mat')
  tmp = load(dso_fn);
  figure(1); clf;
  plot(tmp.data(:,4))
  %data = tmp.data(49060:69245,4);
  
  Nt = round((Tpd-CODEX_ZERO_BINS*4e-9)/dt);
  data = tmp.data(49171+(0:Nt-1),4);
  time = tmp.time(49171+(0:Nt-1));
  plot(data);

  time = reshape(time,[length(time) 1]);
  data = reshape(data,[length(data) 1]);
  dt = mean(diff(time));
  fs = 1/dt;
  fprintf('Sampling frequency: %.12f\n', fs);
else
  error('File format not supported %s.', dso_fn);
end

data = data - mean(data);
plot(data)

if 0
  %% Create hand picked envelope
  
  if 0
    envelope_x = [];
    envelope_y = [];
    xlim([0 2400]);
    while 1
      [envelope_x(end+1),envelope_y(end+1)] = ginput(1);
      xlim([0 2400] + envelope_x(end)-1000);
      xx = xlim;
      ylim([0 max(data)]);
    end
    save('Accum_hand_envelope','envelope_x','envelope_y','data','time');
  end
  
  load('Accum_hand_envelope');
  [envelope_x sort_idxs] = sort(envelope_x);
  envelope_y = envelope_y(sort_idxs);
  
  [envelope_x sort_idxs] = unique(envelope_x);
  envelope_y = envelope_y(sort_idxs);
  
  envelope_y = interp1(envelope_x,envelope_y,1:length(data));
  envelope_x = 1:length(data);
  
  envelope_y(isnan(envelope_y)) = 0;
[B,A] = butter(2,1/100);
envelope_y = filtfilt(B,A,envelope_y);

  figure(1); clf;
  plot(data)
  hold on
  plot(envelope_x,envelope_y,'r')
  plot([3890 3890],[-1 1]*max(data),'k');
  plot(10240+[3890 3890],[-1 1]*max(data),'k');
  hold off;
  
  envelope_x = interp1(1:length(envelope_x),time,envelope_x) - time(3890);
  
  if create_new_codex_file_flag
    %% Creation of Corrected codex.dat file
    
    % Start with flat output (staying away from full scale here to avoid
    % risk of saturation)
    codex = ones(512,1)*35000;
    % Set first and last sample to zero
    % Also: Due to DDS programming method, first 7-8 samples must be 0, in this
    % case we do first 8 samples
    codex([1:9 end]) = 0;
    % Create 0.2 (20%) tukey window, accounting for zeroed out bins (-10)
    % and tukeywin zero-padding (+2), then remove the tukeywin trailing
    % and leading zeros
    Hwin = tukeywin(512+2-10,0.2);
    Hwin = Hwin(2:end-1);
    % Create time axis for output
    time_codex = linspace(0,Tpd,512-10).';
    % Apply correction and window
    codex(10:end-1) = codex(10:end-1) ./ interp1(envelope_x,envelope_y,time_codex) .* Hwin;
    % Normalize (65535 is max DDS output)
    codex = codex / max(codex) * 65000;
    % DDS reads it in backwards
    codex = flipud(codex);
    dlmwrite('codex35_handcorr.dat',codex,'newline','pc','precision','%.0f');
    figure(4); clf;
    plot(codex)
    title('Corrected waveform');
    grid on;
    %     fid = fopen('codex35.dat','w');
    %     fprintf(fid,'%.0f\r\n',codex);
    %     fclose(fid);
  end


end

% fLO = Local oscillator signal for down conversion to baseband
fLO = exp(-j*2*pi*fc*time);

% data_BB = baseband signal
data_BB = data .* fLO;

% low pass filter baseband signal
[B,A] = butter(2,200/2500);
data_BB = filtfilt(B,A,data_BB);

% decimate baseband signal (probably not necessary)
Mt = 10;
data_BB = interpft(data_BB,floor(length(data_BB)/Mt));

Nt = size(data,1);
dt = time(2)-time(1);
Nt2 = floor(Nt/Mt);
dt2 = dt*Mt;
time2 = (0:Nt2-1).'*dt2;
T2 = Nt2*dt2;
df2 = 1/T2;
freq2 = (0:Nt2-1).'*df2;

% Plot results
figure(1); clf;
fig1_plot = 10*log10(abs(fft(data_BB).^2));
h = plot(freq2/1e6, fig1_plot);
xlabel('frequency (MHz)');
ylabel('relative power (dB)');
ylim(max(fig1_plot) + [-30 3]);
grid on;

figure(2); clf;
plot(10*log10(abs(data_BB).^2));
hold on;
plot(10*log10(abs(max(reshape(data(1:floor(length(data)/Mt)*Mt),[Mt floor(length(data)/Mt)]))).^2/Mt)+3,'r');
hold off;
xlabel('sample (range bin)');
ylabel('relative power (dB)');
grid on;

if isempty(truncation_bins)
  truncation_bins = 1:length(data_BB);
end

data_trunc = 10*log10(abs(data_BB(truncation_bins)).^2);

time_trunc = time2(truncation_bins);
% fudge_factor because we truncated the beginning of the chirp due to
% a transient... this would normally be zero.
fudge_factor = 0.048e-6;
time_trunc = time_trunc - time_trunc(1) + fudge_factor;
polynomial_order = 3;
p = polyfit(time_trunc,data_trunc,polynomial_order);
data_fit = polyval(p,time_trunc);

figure(3); clf;
plot(data_trunc);
hold on;
plot(data_fit,'r');
hold off;
grid on;
legend('raw data envelope','polynomial fitted envelope','location','best');

if create_new_codex_file_flag
  %% Creation of Corrected codex.dat file
  
  % Start with flat output (staying away from full scale here to avoid
  % risk of saturation)
  codex = ones(CODEX_LENGTH,1);
  % Set first and last sample to zero
  % Also: Due to DDS programming method, first 7-8 samples must be 0, in this
  % case we do first 8 samples
  codex([1:CODEX_ZERO_BINS-1, end]) = 0;
  Hwin = tukeywin_trim(CODEX_LENGTH-CODEX_ZERO_BINS,TUKEY_WIN);

  % Create time axis for output
  time_codex = linspace(0,Tpd,CODEX_LENGTH-CODEX_ZERO_BINS).';
  
  % Apply correction and window
  codex(CODEX_ZERO_BINS:end-1) = codex(CODEX_ZERO_BINS:end-1) ./ 10.^(polyval(p,time_codex)/20) .* Hwin;
  
  % Normalize (65535 is max DDS output)
  codex = codex / max(codex) * CODEX_MAX_VALUE;
  
  % DDS reads it in backwards
  codex = flipud(codex);
  out_fn = fullfile(out_fn_dir,sprintf('codex_corr%.0f.dat',max(codex)));
  fprintf('Writing %s\n', out_fn)
  dlmwrite(out_fn,codex,'newline','pc','precision','%.0f');
  figure(4); clf;
  plot(codex)
  title('Corrected waveform');
  grid on;
  %     fid = fopen('codex35.dat','w');
  %     fprintf(fid,'%.0f\r\n',codex);
  %     fclose(fid);
end

return;
