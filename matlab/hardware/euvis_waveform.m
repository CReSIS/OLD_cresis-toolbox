% Script euvis_waveform
%
% Subroutine to generate user defined waveform for Euvis AWG801.
%
% Authors: Austin Arnett, John Paden, Tyler Berry

% =========================================================================
%% User Settings
% =========================================================================

base_dir = 'C:\Users\tberry\Desktop\Lecroy\test1';
base_fn = 'C3chirp-lpf-test1-iter';
iteration = 2;

%type specifies whether the GUI or the user file controls the markers
%#type=1 - GUI controls the markers
%#type=5 - text file controls markers
type=5;

% The number of amplitude bits the AWG has, for AWG801 use 11
bits = 11;
data_range = 2^bits;

% Sampling frequency of AWG, for AWG801 use 8e9
fs = 8e9;
dt = 1/fs;

% Chirp parameters
f0 = 500e6;
f1 = 3500e6;
Tpd = 100e-6;


%marker_clock_frequency tells the frequency of the marker clock (50% duty)
%for accurate clock frequency (fs/marker_clock_frequency) must be an integer
marker_clock_frequency = 8e3/1;

%pads the waveform with zeros from the end of the user generated waveform
%to waveform_length. waveform_length must be divisible by the number of
%samples for a clock cycle. (waveform_length*marker_clock_frequency)/fs must be an integer
waveform_length = 1000000;

% =========================================================================
%% Automated Section
% =========================================================================

euvis_waveform_tstart = tic;

if ~exist('base_dir','dir')
  fprintf('Creating output directory %s\n', base_dir);
  mkdir(base_dir);
end

wave_IP = zeros(1,waveform_length);
alpha = (f1-f0)/Tpd;
Nt = Tpd/dt;
time = dt * (0:Nt-1);
freq = f0 + alpha*time;

if iteration > 0
  %% Load corrections file
  corr_fn_name = sprintf('correction_%s%05i.mat',base_fn,iteration-1);
  corr_fn = fullfile(base_dir,corr_fn_name);
  fprintf('  Loading old corrections file %s\n', corr_fn);
  tmp = load(corr_fn);
  
  %% Interpolate correction coefficients to the correct frequencies
  
  % phase_correction: represents the phase of the uncorrected waveform
  % (i.e. we are going to subtract this phase away)
  phase_correction = interp1(tmp.freq_equal,angle(tmp.Hequal),freq,'linear','extrap');
  
  % volt_correction: represents the amplitude waveform that you'll get
  % with no predistortion (i.e. we are going to divide this out)
  volt_correction = interp1(tmp.freq_equal,abs(tmp.Hequal),freq,'linear','extrap');
  
  wave_IP(2:Nt+1) = sin(2*pi*f0*time + pi*alpha*time.^2 - phase_correction) ./ volt_correction;
else
  wave_IP(2:Nt+1) = sin(2*pi*f0*time + pi*alpha*time.^2);
end

% Scale and shift data to go from -1 to +1
wave_IP = 2*wave_IP / (max(wave_IP) - min(wave_IP));

% Create waveform plot
fig_h = figure(1); clf;
fig_fn_name = sprintf('fig_%s%05i_%i.jpg',base_fn,iteration-1,fig_h);
fig_fn = fullfile(base_dir,fig_fn_name);
plot(wave_IP);
saveas(fig_h,fig_fn);

%% Quantize and scale data to 0 to 2047 (11 bit DAC)
wave_quant = round((data_range-1) * (0.5 + 0.5*wave_IP));
wave_quant(wave_quant<0) = 0;
wave_quant(wave_quant>data_range-1) = data_range-1;

% Convert quantized waveform to HEX
wave_hex = dec2hex(wave_quant,3);

%% Creates the "marker clock" used for triggering
marker_clock_length = ceil(1/(marker_clock_frequency*dt));
marker_clock=zeros(1,length(wave_quant));
marker_clock(1:1e4) = 1;

%% write data to file

out_fn_name = sprintf('%s%05i.uda',base_fn,iteration);
out_fn = fullfile(base_dir,out_fn_name);

fprintf('Creating output %s\n', out_fn);
[fid,msg] = fopen(out_fn,'w');
if (fid < 0)
  error(msg);
end

% create header
fprintf(fid,'#type=%i\r\n',type);
fprintf(fid,'#hex=1\r\n',type);

% Fill file with waveform and marker data
for i = 1:length(wave_hex)
  fprintf(fid,'%s %d \t; [%d]\r\n',wave_hex(i,:).',marker_clock(i),i);
end

% Close file so it can be read by GUI program
fclose('all');

fprintf('  Done (%.1f sec)\n', toc(euvis_waveform_tstart));

return;

