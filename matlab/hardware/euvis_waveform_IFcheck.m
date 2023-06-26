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

base_dir = 'C:\Users\tberry\Desktop\Lecroy\IF_test0\';
base_fn = 'C3IF-chirp-iter2-test';
create_corrections_file_en = false;
iteration = 0;
samples_to_use = [1399958 + 60000:3400080 - 60000]; % Leave empty to manually select later

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
BW = f1 - f0;
alpha = BW/Tpd;
time = freq / alpha;
dt = time(2)-time(1);

fprintf('Sampling freq %.1f MHz\n', fs_measured/1e6);

%% Pulse compress data

dso_pc = fft(dso .* hanning(length(dso)));

% Oversample output
Mt = 10;
Nt = length(dso_pc)*Mt;
dt = dt/Mt;
time = dt*(0:Nt-1).';
physical_constants;
range = c/2*time;

df = df/Mt;
freq = df * (0:Nt-1).';

if 0
  % Use range bin axis
  fig_h = figure(13); clf;
  fig_fn_name = sprintf('fig_%s%05i_%i.jpg',base_fn,iteration,fig_h);
  fig_fn = fullfile(base_dir,fig_fn_name);
  plot(lp(ifft(fft(dso_pc)*Mt,Mt*length(dso_pc))),'r');
  grid on;
  xlim([1000 2000]);
  ylim([-70 0]);
  saveas(fig_h,fig_fn);
  
else
  % Use range axis
  fig_h = figure(13); clf;
  fig_fn_name = sprintf('fig_%s%05i_%i.jpg',base_fn,iteration,fig_h);
  fig_fn = fullfile(base_dir,fig_fn_name);
  db = lp(ifft(fft(dso_pc)*Mt,Mt*length(dso_pc)));
  plot(range, db-max(db),'b');
  min_db = min(db);
  ylim([-60 0]);
  xlabel('range (m)');
%   plot(freq/1e6, lp(ifft(fft(dso_pc)*Mt,Mt*length(dso_pc))),'b');
%   xlabel('IF frequency (MHz)');
  ylabel('relative power (dB)');
  grid on;
xlim([437 445])
  ylim manual;
  saveas(fig_h,fig_fn);
  
end

return;
