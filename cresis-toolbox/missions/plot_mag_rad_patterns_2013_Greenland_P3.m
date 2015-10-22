%% Plot the magnitude and angle of measured steering vectors from single radiation pattern for rds_2013_Greenland
tmp = load('C:\Users\y006z286\Documents\scripts\Matlab\20140311_sv\rds_radiation_pattern_single.mat');
% tmp = load('C:\Users\y006z286\Documents\scripts\Matlab\20140311_sv\rds_radiation_pattern_pingpong.mat');

% Slow time, range lines.
rlines = 1:size(tmp.surf_vals{1},2);

color_modes = [0 0 0
    0.3 0.3 0.3
  1 0 0
  1 1 0
  0 1 1
  0 1 0
  0 0 1
  1 0 1];

% figure(3); clf; 
legend_txt = {};
ant_vals = [];

% epri = tmp.epri{1};
% for ant_idx = 2:length(tmp.epri)
%   epri = intersect(epri,tmp.epri{ant_idx});
% end
% [epri rlines_master] = intersect(tmp.epri{1}, epri);
rlines_master = rlines;

% fake_shift = zeros(1,length(tmp.surf_vals));

ref_chan = 4;     % Phase reference channel
surf_bin = 1;     % The bin surface fall into

%% Create Lookup Table
% Bin all the roll measurements into 1 degree large bins
% Determine the roll angle bins that we will put each measurement into
[roll_binned,~,roll_idxs] = unique(round((tmp.roll)*180/pi));

% look_up_table: 
%  First row: roll angle
%  Rows 2 to Nc+1: complex values of surface return for channels 1 to Nc
%  Each column represents a different roll angle from roll_binned
look_up_table = zeros(length(tmp.surf_vals)+1, length(roll_binned));
mag_table = zeros(size(look_up_table(2:end,:)));

look_up_table(1,:) = roll_binned;
for ant_idx = 1:length(tmp.surf_vals)

  % Eventually need to deal with channels that are not EPRI synchronized
  %[epri rlines] = intersect(tmp.epri{ant_idx}, epri);
  
  mag_val = double(abs(tmp.surf_vals{ant_idx}(surf_bin,:)));
  surf_vals = tmp.surf_vals{ant_idx}(surf_bin,:) .* conj(tmp.surf_vals{ref_chan}(surf_bin,:));
  
  % Average all the data falling within each angle/roll bin
  for roll_idx = 1:length(roll_binned)
    look_up_table(ant_idx+1,roll_idx) = mean(surf_vals(roll_idxs == roll_idx));
    mag_table(ant_idx,roll_idx) = mean(mag_val(roll_idxs == roll_idx));
  end
end

nadir_idx = find(roll_binned==0);   % 40
look_up_table = look_up_table .* repmat(conj(look_up_table(:,nadir_idx)),[1 size(look_up_table,2)]);

% Moving average smooth for magnitude
for i = 1:size(mag_table,1)
    mag_table(i,:) = smooth(mag_table(i,:),100);
end
angle_unwrap = zeros(length(tmp.surf_vals), length(roll_binned));
% unwrap_midRef function aims to avoid phase jump from the reference.
angle_unwrap(:,1:nadir_idx) = unwrap_midRef(angle(look_up_table(2:end,1:nadir_idx))); 
angle_unwrap(:,nadir_idx:end) = unwrap(angle(look_up_table(2:end,nadir_idx:end)),[],2);
% angle_unwrap = angle_unwrap(:,5:end);
% x_axis = 1:length(angle_unwrap);

% Moving average smooth for angle
for i = 1:size(angle_unwrap,1)
    angle_unwrap(i,:) = smooth(angle_unwrap(i,:),5);
end

%% Construct the sv table
% Extrapolate to -45~45 range
roll_range = [-45:1:45];
roll_basic = double(roll_binned);
angle_fit = zeros(size(angle_unwrap,1),length(roll_range));
mag_fit = zeros(size(angle_unwrap,1),length(roll_range));

for ant_idx = 1:7
    angle_fit(ant_idx,:) = interp1(roll_basic, angle_unwrap(ant_idx,:), roll_range, 'linear', 'extrap');
    mag_fit(ant_idx,:) = interp1(roll_basic, mag_table(ant_idx,:), roll_range, 'linear', 'extrap');
end
angle_fit = angle_fit - repmat(angle_fit(ref_chan,:),[size(angle_unwrap,1),1]);

diff_nadir_angle = angle_fit(:, nadir_idx) - zeros(size(angle_fit(:, 1)));
angle_fit = angle_fit - repmat(diff_nadir_angle, [1,size(angle_fit(1, :),2)]);

diff = mag_fit(ref_chan, :) - ones(size(mag_fit(ref_chan, :)));
mag_fit = mag_fit - repmat(diff, [size(mag_table,1),1]);

diff_nadir = mag_fit(:, nadir_idx) - ones(size(mag_fit(:, 1)));
mag_fit = mag_fit - repmat(diff_nadir, [1,size(mag_fit(1, :),2)]);

sv_table = mag_fit .* exp(1i * angle_fit);

save('H:/sv_lookuptable.mat','sv_table');


figure(1);clf;
for i = 1:size(mag_table,1)
    h = plot(roll_range,10*log10(abs(sv_table(i,:))));
    set(h,'Color',color_modes(mod(i-1,size(color_modes,1))+1,:));
    hold on;
    legend_txt{i} = sprintf('Antenna %d',i);
end
hold off;
xlim([-45 45])
ylim([-3/10000 +3/10000])
legend(legend_txt);
title('Magnitude distribution of measured steering vectors')
grid on;


figure(2);clf;
plot(roll_range, unwrap_midRef(angle(sv_table)).');
xlim([-45 45])
legend(legend_txt);
title('Angle distribution of measured steering vectors')
grid on
