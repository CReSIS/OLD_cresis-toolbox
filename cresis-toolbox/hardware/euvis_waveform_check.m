% script euvis_waveform_check
%
% Script for checking the performance of the EUVIS 8 GSPS 11 bit DAC (AWG801)
% 1. Reads in a Lecroy DSO .dat file
% 2. Converts to Matlab file .mat
% 3. Pulse compresses and plots
% 4. Plots spectrum
% 5. Optionally creates corrections file used by euvis_waveform.m
%
% For consistency, use these DSO Settings
% 20 us/div, -30 us trigger delay, hold off 110 us
% 50 mV/div
% 400 mV offset
% 10 MS
%
% Authors: Tyler Berry, John Paden
%
% See also: euvis_waveform.m

% =========================================================================
%% User Settings
% =========================================================================

base_dir = 'C:\Users\tberry\Desktop\Lecroy\test1\';
base_fn = 'C3chirp-lpf-test1-iter';
create_corrections_file_en = true;
iteration = 2;
samples_to_use = [1399958:3400080]; % Leave empty to manually select later

plot_en = true;

% Setup pulse compression parameters
Tpd = 100e-6;
f0 = 0.5e9;
f1 = 3.5e9;

% =========================================================================
%% Automated Section
% =========================================================================

%% Convert data file from .dat to .mat for quicker loading


% Load in a batch of files stored by the Lecroy
fn_name = sprintf('%s%05i.mat',base_fn,iteration);
fn = fullfile(base_dir,fn_name);
if exist(fn,'file')
  load(fn);
else
  % MAT file does not exist, create the MAT file and continue
  dat_fn_name = sprintf('%s%05i.dat',base_fn,iteration);
  dat_fn = fullfile(base_dir,dat_fn_name);
  
  fprintf('Loading file %s\n', dat_fn);
  % Load the data from the file
  A = load(dat_fn);
  
  if isempty(samples_to_use)
    plot(A(:,2));
    fprintf('Set samples_to_use and then run dbcont\n');
    keyboard
  end
  time = A(samples_to_use,1);
  dso = A(samples_to_use,2);
  
  fprintf('  Saving output %s\n', fn);
  save(fn,'time','dso');
end

if plot_en
  fig_h = figure(11); clf;
  fig_fn_name = sprintf('fig_%s%05i_%i.jpg',base_fn,iteration,fig_h);
  fig_fn = fullfile(base_dir,fig_fn_name);
  plot(time,dso);
  title(sprintf('DSO captured %s\n', fn));
  xlabel('time (sec)');
  ylabel('voltage (V)');
  % xlim([-2e-9 0]+100e-6)
  grid on;
  saveas(fig_h,fig_fn);
end

%% Setup frequency axis
dt = time(2)-time(1);
fs_measured = 1/dt;
Nt = length(dso);
T = dt*Nt;
df = 1/T;
freq = df * (0:Nt-1).';

fprintf('Sampling freq %.1f MHz\n', fs_measured/1e6);

%% Pulse compress data

% Step 1: EUVIS clock is inaccurate and substantially different from the
% DSO. We can search for the offset, but we have already done this and know
% the clock is off by 1.22 MHz.

% dfs = delta-fs (Hz)
dfs_range = 1220000;
% dfs_range = linspace(1218000,1225000,21);

for dfs = dfs_range
  fs = 20.000e9 + dfs;
  
  % Create time axis and derived pulse compression parameters
  dt = 1/fs;
  alpha = (f1-f0)/Tpd;
  Nt = length(dso);
  Nt_ref = Tpd*fs;
  time = dt*(0:Nt-1).' + 50e-9;
  
  % Create the reference pulse compression waveform
  ref = exp(j*2*pi*f0*time + j*pi*alpha*time.^2);
  ref(time < 0 & time>Tpd) = 0;
  ref = fft(ref,2*Nt);
  
  % Convert DSO data to frequency domain
  dso_freq = fft(dso,2*Nt);
  
  % Plot Frequency Spectrum
  fig_h = figure(12); clf;
  fig_fn_name = sprintf('fig_%s%05i_%i.jpg',base_fn,iteration,fig_h);
  fig_fn = fullfile(base_dir,fig_fn_name);
  plot(lp(ref));
  hold on;
  plot(lp(dso_freq),'r');
  hold off;
  grid on;
  saveas(fig_h,fig_fn);
  
  % Create new time/frequency axis for pulse compressed data
  T = dt*2*Nt;
  df = 1/T;
  freq = df * (0:2*Nt-1).';
  
  % Filter and complex base band the data by selecting and using just
  % the frequencies of interest
  good_idxs = find(freq > f0 & freq < f1);
  ref = ref(good_idxs);
  dso_freq = dso_freq(good_idxs);
  
  % Apply a window to the reference
  ref_wind = conj(ref) .* hanning(length(ref));
  
  % Perform pulse compression in the frequency domain
  dso_pc = ifft(dso_freq .* ref_wind);
  
  % Normalize output
  dso_pc = dso_pc./max(dso_pc);
  
  % Oversample output
  Nt = length(dso_pc)*10;
  fs = df*Nt;
  dt = 1/fs;
  time = dt*(0:Nt-1).';
  physical_constants;
  range = c/2/sqrt(1.53)*time;
  
  if 0
    % Use range bin axis
    fig_h = figure(13); clf;
    fig_fn_name = sprintf('fig_%s%05i_%i.jpg',base_fn,iteration,fig_h);
    fig_fn = fullfile(base_dir,fig_fn_name);
    plot(lp(ifft(fft(dso_pc)*10,10*length(dso_pc))),'r');
    grid on;
    xlim([1000 2000]);
    ylim([-70 0]);
    saveas(fig_h,fig_fn);
    
  else
    % Use range axis
    fig_h = figure(13); hold on;
    fig_fn_name = sprintf('fig_%s%05i_%i.jpg',base_fn,iteration,fig_h);
    fig_fn = fullfile(base_dir,fig_fn_name);
    plot(range, lp(ifft(fft(dso_pc)*10,10*length(dso_pc))),'b');
    xlabel('range (m)');
    ylabel('relative power (dB)');
    grid on;
    xlim([4 9])
    ylim([-70 0]);
    saveas(fig_h,fig_fn);
    
  end
  
  if length(dfs_range) > 1
    % If looping through lots of delta-fs values, pause after each one
    % except for the last one.
    dfs
    pause;
  end
end

if create_corrections_file_en
  %% Create correction (equalization) coefficients to be used by euvis_waveform.m
  
  % Create low pass filter (LPF) coefficients which will be used to average out
  % high frequency noise in data
  [Bequal,Aequal] = butter(2,1/5000);
  % Median filter and then LPF the power spectrum
  Hequal = filtfilt(Bequal,Aequal,medfilt1(abs(dso_freq).^2,101));
  % Replace/overwrite transient at beginning
  Hequal(1:10000) = Hequal(10000);
  %   Hequal(end-5000:end) = Hequal(end-5000);
  % Convert power to amplitude scale
  Hequal = sqrt(Hequal);
  %Get the frequency axis for the equalization coefficients
  freq_equal = freq(good_idxs);
  % Plot results for debugging
  fig_h = figure(14); clf;
  fig_fn_name = sprintf('fig_%s%05i_%i.jpg',base_fn,iteration,fig_h);
  fig_fn = fullfile(base_dir,fig_fn_name);
  plot(freq_equal/1e9,lp(dso_freq,2),'b')
  hold on;
  plot(freq_equal/1e9,lp(Hequal,2),'r')
  hold off;
  saveas(fig_h,fig_fn);
  
  % Divide out the refernce function, then find the angle and median filter
  % the angle to remove outliers, unwrap the phase, detrend the phase,
  % and then LPF the phase.
  %   Phase unwrapping:
  phase_corr = filtfilt(Bequal,Aequal,detrend(unwrap(medfilt1(angle(dso_freq ./ ref),51))));
  % Replace/overwrite transient at beginning
  phase_corr(1:10000) = phase_corr(10000);
  phase_corr = phase_corr - phase_corr(1);
  Hequal = Hequal .* exp(j*phase_corr);
  
  fig_h = figure(15); clf;
  fig_fn_name = sprintf('fig_%s%05i_%i.jpg',base_fn,iteration,fig_h);
  fig_fn = fullfile(base_dir,fig_fn_name);
  plot(freq_equal/1e9,detrend(unwrap(medfilt1(angle(dso_freq ./ ref),51)))*180/pi)
  hold on;
  plot(freq_equal/1e9,angle(Hequal)*180/pi,'r');
  hold off;
  %   Hequal(end-5000:end) = Hequal(end-5000);
  saveas(fig_h,fig_fn);
  
  if iteration >= 1
    % If the second iteration or later, we need to load in the old
    % correction coefficients to determine the new coefficients
    corr_fn_name = sprintf('correction_%s%05i.mat',base_fn,iteration-1);
    corr_fn = fullfile(base_dir,corr_fn_name);
    fprintf('  Loading old corrections file %s\n', corr_fn);
    old = load(corr_fn,'freq_equal','Hequal');
    Hequal = Hequal .* old.Hequal ./ max(abs(old.Hequal));
  end
  
  %% The EUVIS AWG has problems generating signals smaller than 50% of the full scale
  % so we clamp corrections to 50% or above here. Note: the amplitude will be
  % 1/voltage_correction so we have to invert everything
  volt_correction = abs(Hequal);
  % Normalize voltage correction to a maximum of 100%
  volt_correction = volt_correction / min(volt_correction);
  % Find all indices below 50%, these are the ones we want to clamp and
  % set those to slightly higher than 50% to avoid rounding/sampling issues
  clamp_idxs = find(volt_correction > 1/0.53);
  volt_correction(clamp_idxs) = 1/0.53;
  Hequal = volt_correction .* exp(j*angle(Hequal));
  
  corr_fn_name = sprintf('correction_%s%05i.mat',base_fn,iteration);
  corr_fn = fullfile(base_dir,corr_fn_name);
  fprintf('  Saving corrections file %s\n', corr_fn);
  save(corr_fn,'freq_equal','Hequal');
end

return;
