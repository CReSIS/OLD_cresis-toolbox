% script plot_rad_patterns
%
% First generate the individual element patterns using
% generate_complex_svLUT.m

if 0
  %% 2014_Antarctica_DC8
  
elseif 0
  %% 2015_Greenland_C130
  % fns = list of filenames to load and concatenate data from
  fn = 'D:\rds\2015_Greenland_C130\CSARP_noise\surf_20150313_14.mat';
  elements_fn = 'D:\rds\2015_Greenland_C130\CSARP_noise\sv_table_2015_Greenland_C130.mat';
  output_fn = 'D:\rds\2015_Greenland_C130\CSARP_noise\combined_pattern_2015_Greenland_C130.mat';
  
  % roll_to_ant_mapping(ANTENNA) = index into surf_vals containing the
  % combined antenna transmit data
  roll_to_ant_mapping = [5];
  
  % rlines = restrict which range lines are used (or use all)
  rlines = [];
  
  ref_idx = 1;     % Receive channel index into surf_vals from fn
  ref_ant = 1;     % Receive channel index into sv_deviation_approx from elements_fn
  surf_bins = 6;     % The relative range bin into surf_vals that we will use for extracting values from
  fc = (180e6 + 450e6)/2; % Center frequency
  
  good_mask_min_samples = 100;
  good_mask_min_angle = -50;
  good_mask_max_angle = 40;
  
  extra_degrees_of_freedom = 6; % spatial filter fitting (# of antennas + this many)
  
elseif 1
  %% 2015_Greenland_Polar6
  % fns = list of filenames to load and concatenate data from
  fn = 'J:\rds\2015_Greenland_Polar6\CSARP_noise\surf_20150911_17.mat';
  elements_fn = 'J:\rds\2015_Greenland_Polar6\CSARP_noise\sv_table_2015_Greenland_Polar6.mat';
  output_fn = 'J:\rds\2015_Greenland_Polar6\CSARP_noise\combined_pattern_2015_Greenland_Polar6.mat';
  
  % roll_to_ant_mapping(ANTENNA) = index into surf_vals containing the
  % combined antenna transmit data
  roll_to_ant_mapping = [2];
  
  % rlines = restrict which range lines are used (or use all)
  rlines = [];
  
  ref_idx = 1;     % Receive channel index into surf_vals from fn (used only for generating complex rad patterns)
  ref_ant = 1;     % Receive channel index into sv_deviation_approx from elements_fn (should correspond to the receiver that is used in fn)
  surf_bins = 11;     % The relative range bin into surf_vals that we will use for extracting values from
  fc = (165e6 + 520e6)/2; % Center frequency  
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -35;
  good_mask_max_angle = 35; 
  
  degrees_of_freedom = 4; % spatial filter fitting (usually around the # of antennas)  
  
end

%% Load data
data = load(fn);

if 0
  close all
  %% DEBUG
  figure(1); clf;
  subplot(3,1,1);
  plot(lp(data.surf_vals(surf_bins,:,roll_to_ant_mapping)),'.')
  a1 = gca;
  subplot(3,1,2);
  plot(angle(data.surf_vals(surf_bins,:,roll_to_ant_mapping) .* conj(data.surf_vals(surf_bins,:,roll_to_ant_mapping)) ),'.')
  a2 = gca;
  subplot(3,1,3);
  plot(data.roll);
  a3 = gca;
  figure(2); clf;
  imagesc(lp(data.surf_vals(:,:,roll_to_ant_mapping)));
  a4 = gca;
  linkaxes([a1 a2 a3 a4],'x');
  return;
end

%% Create Lookup Table
% Bin all the roll measurements into 1 degree large bins
% Determine the roll angle bins that we will put each measurement into
[roll_binned,~,roll_idxs] = unique(round((data.roll)*180/pi));

figure(1); clf;
hist(roll_binned(roll_idxs), -90:90);
xlabel('Angle (deg)');
ylabel('Number of measurements');

% sv_table: 
%  First row: roll angle
%  Rows 2 to Nc+1: complex values of surface return for channels 1 to Nc
%  Each column represents a different roll angle from roll_binned
sv_table = zeros(length(roll_to_ant_mapping), length(roll_binned));
power_table = zeros(size(sv_table));

for ant = 1:length(roll_to_ant_mapping)
  input_idx = roll_to_ant_mapping(ant);
  
  % Eventually need to deal with channels that are not EPRI synchronized
  %[epri rlines] = intersect(data.epri{ant_idx}, epri);
  
  powers = mean(abs(double(data.surf_vals(surf_bins,:,input_idx))).^2,1);
  complex_vals = mean(data.surf_vals(surf_bins,:,input_idx) .* exp(-1i*angle(data.surf_vals(surf_bins,:,ref_idx))),1);
  
  % Average all the data falling within each angle/roll bin
  for roll_idx = 1:length(roll_binned)
    sv_table(ant,roll_idx) = mean(complex_vals(roll_idxs == roll_idx));
    power_table(ant,roll_idx) = mean(powers(roll_idxs == roll_idx));
  end
end

%% Different ways to just keep good measurements
N = hist(roll_binned(roll_idxs), roll_binned);
good_mask = N > good_mask_min_samples;
good_mask = good_mask & roll_binned > good_mask_min_angle & roll_binned < good_mask_max_angle;
plot(roll_binned, good_mask)
roll_binned = roll_binned(good_mask);
sv_table = sv_table(:,good_mask);
power_table = power_table(:,good_mask);

%% Receiver equalization (force sv to be all ones at nadir)
nadir_idx = find(roll_binned==0);
rx_equalization = sv_table(:,nadir_idx);
fprintf('Equalization (deg):\n');
fprintf('%5.1f ', angle(rx_equalization)*180/pi);
fprintf('\n');
fprintf('Equalization (dB):\n');
fprintf('%5.1f ', 10*log10(abs(rx_equalization)));
fprintf('\n');
sv_table = sv_table ./ repmat(rx_equalization,[1 size(sv_table,2)]);

%% Load the individual element pattern
sv_elements = load(elements_fn);

sv_table = sv_table ./ sv_elements.sv_deviation_approx(ref_ant,:);

ky_relative = sind(roll_binned);
[B,A] = invfreqz(sv_table,ky_relative,'complex',degrees_of_freedom,0);
sv_table_approx = freqz(B,A,ky_relative);

sv_table = sv_table/max(sv_table_approx);
sv_table_approx = sv_table_approx/max(sv_table_approx);

figure(2); clf;
plot(roll_binned,10*log10(abs(sv_table_approx).^2))
hold on
plot(roll_binned,10*log10(abs(sv_table).^2),'.')
grid on;
ylabel('Power (dB)');
xlabel('Angle (deg)');
legend('Fit','Measured');

save(output_fn,'roll_binned','sv_table','sv_table_approx');

return;
