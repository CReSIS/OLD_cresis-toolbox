% script generate_complex_svLUT.m
%
% Generates complex steering vector lookup table using roll data collected
% with an array system.

if 0
  %% 2014_Antarctica_DC8
  % fn = output from coh_noise surf tracker
  
elseif 0
  %% 2015_Greenland_C130
  % fn = output from coh_noise surf tracker
  fn = 'D:\rds\2015_Greenland_LC130\CSARP_noise\surf_20150313_14.mat';
  output_fn = 'D:\rds\2015_Greenland_LC130\CSARP_noise\sv_table_2015_Greenland_C130.mat';
  
  % roll_to_ant_mapping(ANTENNA) = index into surf_vals containing this
  % antenna's data (transmit and receive with the same antenna)
  roll_to_ant_mapping = [1 4]; % 2015_Greenland_C130
  
  % rlines = restrict which range lines are used (or use all)
  rlines = [];
  
  ref_ant = 1;     % Antenna phase reference channel (usually a center element)
  surf_bins = 1;     % The relative range bin into surf_vals that we will use for extracting values from
  fc = (180e6 + 450e6)/2; % Center frequency
  
  good_mask_min_samples = 100;
  good_mask_min_angle = -50;
  good_mask_max_angle = 40;
  
  extra_degrees_of_freedom = 6; % spatial filter fitting (# of antennas + this many)
  
elseif 0
  %% 2013_Antarctica_Basler
  % fn = output from coh_noise surf tracker
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2013_Antarctica_Basler/CSARP_noise/surf_20131216_05.mat';
  output_fn = '/cresis/snfs1/dataproducts/ct_data/rds/2013_Antarctica_Basler/CSARP_noise/sv_table_2013_Antarctica_Basler.mat';
  
  % roll_to_ant_mapping(ANTENNA) = index into surf_vals containing this
  % antenna's data (transmit and receive with the same antenna)
  roll_to_ant_mapping = [1:8]; % 2013_Antarctica_Basler
%   roll_to_ant_mapping = [9:16]; % 2013_Antarctica_Basler
%   roll_to_ant_mapping = [16 12 13 9 15 11 14 10]; % 2013_Antarctica_Basler
  
  % rlines = restrict which range lines are used (or use all)
  rlines = [];
  
  ref_ant = 4;     % Antenna phase reference channel (usually a center element)
  surf_bins = 6;     % The relative range bin into surf_vals that we will use for extracting values from
  fc = (200e6 + 450e6)/2; % Center frequency
  
  good_mask_min_samples = 27;
  good_mask_min_angle = -57;
  good_mask_max_angle = 51;
  
  extra_degrees_of_freedom = -3; % spatial filter fitting (# of antennas + this many)
  
elseif 1
  %% 2015_Greenland_Polar6
  % fn = output from coh_noise surf tracker
  fn = 'J:/rds/2015_Greenland_Polar6/CSARP_noise/surf_20150911_15.mat';
  output_fn = 'J:/rds/2015_Greenland_Polar6/CSARP_noise/sv_table_2015_Greenland_Polar6.mat';
  
  % roll_to_ant_mapping(ANTENNA) = index into surf_vals containing this
  % antenna's data (transmit and receive with the same antenna)
  roll_to_ant_mapping = [1:8]; % 2015_Greenland_Polar6
%   roll_to_ant_mapping = [9:16]; % 2013_Antarctica_Basler
%   roll_to_ant_mapping = [16 12 13 9 15 11 14 10]; % 2013_Antarctica_Basler
  
  % rlines = restrict which range lines are used (or use all)
  rlines = [];
  
  ref_ant = 4;     % Antenna phase reference channel (usually a center element)
  surf_bins = 6;     % The relative range bin into surf_vals that we will use for extracting values from
  fc = (165e6 + 510e6)/2; % Center frequency
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -35;
  good_mask_max_angle = 35; 
  
  extra_degrees_of_freedom = -3; % spatial filter fitting (# of antennas + this many)
    
end

%% Load data
data = load(fn);

Hchan = chebwin(8,30).';
for chan = roll_to_ant_mapping
  data.surf_vals(:,:,chan) = data.surf_vals(:,:,chan) ./ Hchan(chan);
end

if 0
  close all
  %% DEBUG
  figure(1); clf;
  subplot(3,1,1);
  plot(lp(data.surf_vals(surf_bins,:,1)),'.')
  a1 = gca;
  subplot(3,1,2);
  plot(angle(data.surf_vals(surf_bins,:,4) .* conj(data.surf_vals(surf_bins,:,1)) ),'.')
  a2 = gca;
  subplot(3,1,3);
  plot(data.roll);
  a3 = gca;
  linkaxes([a1 a2 a3],'x');
  figure(2); clf;
  imagesc(lp(data.surf_vals(:,:,1)));
  return;
end

%% Retrack surface
if 1
  ref_idx = roll_to_ant_mapping(ref_ant);
  
  threshold = lp(mean(abs(data.surf_vals(1,:,ref_idx)).^2)) + 4;
  surf_bin = NaN*zeros(1,size(data.surf_vals,2));
  ml_data = lp(fir_dec(abs(data.surf_vals(:,:,ref_idx)).^2,ones(1,5)/5,1));
  for rline = 1:size(data.surf_vals,2)
    cur_threshold = max([ml_data(1,rline)+7; ml_data(:,rline)-13]);
    tmp = find(ml_data(:,rline) > cur_threshold,1);
    if ~isempty(tmp)
      [~,max_offset] = max(ml_data(tmp+(0:2),rline));
      tmp = tmp-1 + max_offset;
      surf_bin(rline) = tmp;
    end
  end
  
  figure(1); clf;
  imagesc(ml_data);
  hold on;
  plot(surf_bin);
  
  for rline = 1:size(data.surf_vals,2)
    if ~isnan(surf_bin(rline))
      data.surf_vals(:,rline,:) = circshift(data.surf_vals(:,rline,:),[6-surf_bin(rline) 0 0]);
    end
  end
  
  ml_data = lp(fir_dec(abs(data.surf_vals(:,:,ref_idx)).^2,ones(1,5)/5,1));
  figure(2); clf;
  imagesc(ml_data);
  
  keyboard
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
  ref_idx = roll_to_ant_mapping(ref_ant);
  
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

if 0
  close all
  %% DEBUG
  figure(1); clf;
  imagesc(roll_binned,[],lp(power_table));
  colorbar;
  figure(2); clf;
  imagesc(roll_binned,[],angle(sv_table));
  colorbar;
  figure(30); clf;
  imagesc(lp(data.surf_vals(:,:,8) ))
  return;
end

%% Different ways to just keep good measurements
N = hist(roll_binned(roll_idxs), roll_binned);
good_mask = N > good_mask_min_samples;
good_mask = good_mask & roll_binned > good_mask_min_angle & roll_binned < good_mask_max_angle;
if 0
  %% DEBUG
  plot(roll_binned, good_mask)
end
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

%% Load the ideal steering vectors
physical_constants;
clear phase_centers;
for tx_chan = 1:length(roll_to_ant_mapping)
  tx_weights = zeros(1,length(roll_to_ant_mapping));
  tx_weights(tx_chan) = 1;
  phase_centers(:,tx_chan) = lever_arm(data.param_analysis, tx_weights, tx_chan);
end
[theta,sv_ideal] = array_proc_sv({'theta' roll_binned/180*pi}, fc, phase_centers(2,:).', phase_centers(3,:).');
sv_ideal = sv_ideal ./ repmat(sv_ideal(ref_ant,:),[size(sv_ideal,1) 1]);
rx_equalization = sv_ideal(:,nadir_idx);
sv_ideal = sv_ideal ./ repmat(rx_equalization,[1 size(sv_ideal,2)]);

sv_deviation = sv_table ./ sv_ideal;
figure(2); clf;
imagesc(roll_binned,[],angle(sv_deviation)*180/pi)
h = colorbar;
set(get(h,'YLabel'),'String','Phase (deg)')
xlabel('Angle (deg)');
ylabel('Antenna');

figure(3); clf;
plot(roll_binned,angle(sv_deviation)*180/pi)
ylabel('Phase deviation (deg)')
xlabel('Angle (deg)');
legend_str = {};
for tx_chan = 1:length(roll_to_ant_mapping)
  legend_str{tx_chan} = sprintf('Ant %d', tx_chan);
end
legend(legend_str,'location','best');

figure(4); clf;
plot(roll_binned,10*log10(abs(sv_deviation).^2))
ylabel('Power (dB)');
xlabel('Angle (deg)');
for tx_chan = 1:length(roll_to_ant_mapping)
  legend_str{tx_chan} = sprintf('Ant %d', tx_chan);
end
legend(legend_str,'location','best');

ky_relative = sind(roll_binned);
sv_deviation_approx = [];
for ant = 1:size(sv_deviation,1)
  [B,A] = invfreqz(sv_deviation(ant,:),ky_relative,'complex',length(roll_to_ant_mapping)+extra_degrees_of_freedom,0);
  sv_deviation_approx(ant,:) = freqz(B,A,ky_relative);
  
  figure(100+ant); clf;
  plot(roll_binned,10*log10(abs(sv_deviation_approx(ant,:)).^2))
  hold on;
  plot(roll_binned,10*log10(abs(sv_deviation(ant,:)).^2),'.')
  grid on;
  ylabel('Power (dB)');
  xlabel('Angle (deg)');
  legend('Approximated','Measured','Location','best');
  set(100+ant,'WindowStyle','docked')
  title(sprintf('Ant %d', ant));
  
  figure(200+ant); clf;
  plot(roll_binned,180/pi*angle(sv_deviation_approx(ant,:)))
  hold on;
  plot(roll_binned,180/pi*angle(sv_deviation(ant,:)),'.')
  grid on;
  ylabel('Phase deviation (deg)')
  xlabel('Angle (deg)');
  legend('Approximated','Measured','Location','best');
  set(200+ant,'WindowStyle','docked')
  title(sprintf('Ant %d', ant));
end

roll = roll_binned;
sv_deviation_approx_angle = unwrap(angle(sv_deviation_approx).').';
sv_deviation_approx_angle = sv_deviation_approx_angle - repmat(sv_deviation_approx_angle(:,nadir_idx),[1 size(sv_deviation_approx_angle,2)]);
sv_deviation_approx = sqrt(abs(sv_deviation_approx)) .* exp(1i*0.5*sv_deviation_approx_angle);

sv_ideal_angle = unwrap(angle(sv_ideal).').';
sv_ideal_angle = sv_ideal_angle - repmat(sv_ideal_angle(:,nadir_idx),[1 size(sv_ideal_angle,2)]);
sv_ideal = sqrt(abs(sv_ideal)) .* exp(1i*0.5*sv_ideal_angle);

save(output_fn,'roll','sv_deviation_approx','sv_ideal');

figure(5); clf;
plot(roll,10*log10(abs(sv_deviation_approx .* sv_ideal).^2));
hold on;
% plot(roll,10*log10(abs(sqrt(sv_table)).^2),'.');
grid on;
ylabel('Power (dB)');
xlabel('Angle (deg)');
for tx_chan = 1:length(roll_to_ant_mapping)
  legend_str{tx_chan} = sprintf('Ant %d', tx_chan);
end
legend(legend_str,'location','best');
% xlim([-45 45]);

nadir_idx = find(roll==0);
final_phase = unwrap(angle(sv_deviation_approx .* sv_ideal).').';
final_phase = final_phase - repmat(final_phase(:,nadir_idx),[1 size(final_phase,2)]);
measured = unwrap(angle(sv_table).').'/2;
measured = measured - repmat(measured(:,nadir_idx),[1 size(measured,2)]);
ideal = unwrap(angle(sv_ideal).').';
ideal = ideal - repmat(ideal(:,nadir_idx),[1 size(ideal,2)]);
figure(6); clf;
plot(roll,final_phase*180/pi);
hold on;
% plot(roll,measured*180/pi,'.');
plot(roll,ideal*180/pi,'k-');
grid on;
ylabel('Phase (deg)')
xlabel('Angle (deg)');
for tx_chan = 1:length(roll_to_ant_mapping)
  legend_str{tx_chan} = sprintf('Ant %d', tx_chan);
end
legend(legend_str,'location','best');
% xlim([-45 45]);
