% script plot_rad_patterns
%
% First generate the individual element patterns using
% generate_complex_svLUT.m

% fns = list of filenames to load and concatenate data from
fns = {'D:\roll_measurements_2014_Antarctica_DC8.mat','D:\roll_measurements_2014_Antarctica_DC8_2.mat'};
elements_fn = 'D:\sv_table_2014_Antarctica_DC8.mat';
output_fn = 'D:\combined_pattern_2014_Antarctica_DC8.mat';

% roll_to_ant_mapping(ANTENNA) = index into surf_vals containing the
% combined antenna transmit data
roll_to_ant_mapping = [7];

% rlines = restrict which range lines are used (or use all)
rlines = [];

ref_ant = 1;     % Antenna phase reference channel (usually a center element)
surf_bins = 1:5;     % The relative range bin into surf_vals that we will use for extracting values from
fc = (165e6 + 215e6)/2; % Center frequency

%% Load data
data = [];
for fn_idx = 1:length(fns)
  tmp = load(fns{1});
  if fn_idx == 1
    data.roll = [];
    data.epri = cell(1,length(tmp.epri));
    data.surf_vals = cell(1,length(tmp.surf_vals));
  end
  data.roll = cat(2,data.roll,tmp.roll);
  for idx = 1:length(tmp.epri)
    data.epri{idx} = cat(2,data.epri{idx},tmp.epri{idx});
  end
  for idx = 1:length(tmp.epri)
    data.surf_vals{idx} = cat(2,data.surf_vals{idx},tmp.surf_vals{idx});
  end
end

%% Align each of the surf_vals results
% epri = data.epri{1};
% for ant_idx = 2:length(data.epri)
%   epri = intersect(epri,data.epri{ant_idx});
% end
% [epri rlines] = intersect(data.epri{1}, epri);

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
  
  powers = mean(abs(double(data.surf_vals{input_idx}(surf_bins,:))).^2,1);
  complex_vals = mean(data.surf_vals{input_idx}(surf_bins,:) .* exp(-1i*angle(data.surf_vals{ref_idx}(surf_bins,:))),1);
  
  % Average all the data falling within each angle/roll bin
  for roll_idx = 1:length(roll_binned)
    sv_table(ant,roll_idx) = mean(complex_vals(roll_idxs == roll_idx));
    power_table(ant,roll_idx) = mean(powers(roll_idxs == roll_idx));
  end
end

%% Different ways to just keep good measurements
N = hist(roll_binned(roll_idxs), roll_binned);
good_mask = N > 100;
% good_mask = abs(roll_binned) < 45;
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

sv_table = sv_table ./ sv_elements.sv_deviation_approx(2,:);

ky_relative = sind(roll);
[B,A] = invfreqz(sv_table,ky_relative,'complex',7,0);
sv_table_approx = freqz(B,A,ky_relative);

sv_table = sv_table/max(sv_table_approx);
sv_table_approx = sv_table_approx/max(sv_table_approx);

figure(2); clf;
plot(roll,10*log10(abs(sv_table_approx).^2))
hold on
plot(roll,10*log10(abs(sv_table).^2),'.')
grid on;
ylabel('Power (dB)');
xlabel('Angle (deg)');
legend('Fit','Measured');

save(output_fn,'roll','sv_table','sv_table_approx');

return;
