% script read_ins_applanix.m
%
% Reads in Applanix POSPAC software display plot export file.
% Specifically designed for X,Y,Z,Total acceleration being selected in tabular form
% and then exported.
%
% Author: John Paden

fn = 'C:\tmp\Forward Processed Trajectory, Reference Frame.txt';
save_dir = 'C:\tmp';

fid = fopen(fn);

A = textscan(fid, '%f %f %f %f %f','HeaderLines',1);
feof(fid);

fclose(fid);

%% All events
time = A{1} - A{1}(1);
x = A{2};
y = A{3};
z = A{4};
total = A{5};

figure(1); clf;
plot(time,A{5},'b')
xlabel('Time (sec)');
ylabel('Acceleration (g)');
hold on;
ylims = ylim;
plot([1805 1831.5 1831.5 1805 1805], ylims([1 1 2 2 1]),'r');
plot([3280 3283 3283 3280 3280], ylims([1 1 2 2 1]),'g');
plot(time,A{5},'b')
legend('Total XYZ', 'Event 1', 'Event 2');
saveas(1,fullfile(save_dir,'Total.jpg'));

%% Event 1

good_mask = time>2800 & time<2980;

% figure(1); clf;
% plot(time(good_mask),x(good_mask))
% hold on;
% plot(time(good_mask),y(good_mask))
% plot(time(good_mask),z(good_mask))

dt = time(2)-time(1)
fs = 1/dt

figure(2); clf;
periodogram(x(good_mask),[],'onesided',512,fs)
title('X Vibration PSD');
figure(3); clf;
periodogram(y(good_mask),[],'onesided',512,fs)
title('Y Vibration PSD');
figure(4); clf;
periodogram(z(good_mask),[],'onesided',512,fs)
title('Z Vibration PSD');
saveas(2,fullfile(save_dir,'Background1_X.jpg'));
saveas(3,fullfile(save_dir,'Background1_Y.jpg'));
saveas(4,fullfile(save_dir,'Background1_Z.jpg'));

%% Event 1

good_mask = time>1805 & time<1831.5;

% figure(1); clf;
% plot(time(good_mask),x(good_mask))
% hold on;
% plot(time(good_mask),y(good_mask))
% plot(time(good_mask),z(good_mask))

dt = time(2)-time(1)
fs = 1/dt

figure(2); clf;
periodogram(x(good_mask),[],'onesided',512,fs)
title('X Vibration PSD');
figure(3); clf;
periodogram(y(good_mask),[],'onesided',512,fs)
title('Y Vibration PSD');
figure(4); clf;
periodogram(z(good_mask),[],'onesided',512,fs)
title('Z Vibration PSD');
saveas(2,fullfile(save_dir,'Event1_X.jpg'));
saveas(3,fullfile(save_dir,'Event1_Y.jpg'));
saveas(4,fullfile(save_dir,'Event1_Z.jpg'));

%% Event 2

good_mask = time>3280 & time<3283;

figure(1); clf;
plot(time(good_mask),x(good_mask))
hold on;
plot(time(good_mask),y(good_mask))
plot(time(good_mask),z(good_mask))

dt = time(2)-time(1)
fs = 1/dt

figure(2); clf;
periodogram(x(good_mask),[],'onesided',512,fs)
title('X Vibration PSD');
figure(3); clf;
periodogram(y(good_mask),[],'onesided',512,fs)
title('Y Vibration PSD');
figure(4); clf;
periodogram(z(good_mask),[],'onesided',512,fs)
title('Z Vibration PSD');
saveas(2,fullfile(save_dir,'Event2_X.jpg'));
saveas(3,fullfile(save_dir,'Event2_Y.jpg'));
saveas(4,fullfile(save_dir,'Event2_Z.jpg'));
