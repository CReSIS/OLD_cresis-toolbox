physical_constants;

clear fns;
fns{1} = 'D:\data\snow\S11_testflight\KUBandSnow\S11_FlightTest_KUBand_Tx.s1p';
fns{2} = 'D:\data\snow\S11_testflight\KUBandSnow\S11_FlightTest_KUBand_Rx.s1p';
VF = 0.8;

  
figure(1); clf;
[freq, data, freq_noise, data_noise, Zo] = SXPParse(fns{1});
S11 = squeeze(data(1,1,:));
df = freq(2)-freq(1);
Nt = length(S11);
BW = Nt*df;
dt = 1/BW;
time = dt*(0:(Nt-1)).';
range = time*(c/2*VF);

plot(freq/1e6, lp(S11));

figure(2); clf;
plot(range,lp(ifft(S11 .* hanning(length(S11)))))

figure(1);
hold on;
[freq, data, freq_noise, data_noise, Zo] = SXPParse(fns{2});
S11 = squeeze(data(1,1,:));
plot(freq/1e6, lp(S11),'r');
hold off;
title('Ku-band');
xlabel('frequency (MHz)')
ylabel('relative power (dB)');
grid on;
legend('Tx','Rx');

figure(2);
hold on;
plot(range,lp(ifft(S11 .* hanning(length(S11)))),'r')
hold off;
title(sprintf('Ku-band (velocity factor %.2f)',VF));
xlabel('range (m)')
ylabel('relative power (dB)');
grid on;
legend('Tx','Rx');




clear fns;
fns{1} = 'D:\data\snow\S11_testflight\KUBandSnow\S11_FlightTest_Snow_Tx.s1p';
fns{2} = 'D:\data\snow\S11_testflight\KUBandSnow\S11_FlightTest_Snow_Rx.s1p';

figure(3); clf;
[freq, data, freq_noise, data_noise, Zo] = SXPParse(fns{1});
df = freq(2)-freq(1);
Nt = length(S11);
BW = Nt*df;
dt = 1/BW;
time = dt*(0:(Nt-1)).';
range = time*(c/2*VF);

S11 = squeeze(data(1,1,:));
plot(freq/1e6, lp(S11));

figure(4); clf;
plot(range,lp(ifft(S11 .* hanning(length(S11)))))

figure(3);
hold on;
[freq, data, freq_noise, data_noise, Zo] = SXPParse(fns{2});
S11 = squeeze(data(1,1,:));
plot(freq/1e6, lp(S11),'r');
hold off;
title('Snow');
xlabel('frequency (MHz)')
ylabel('relative power (dB)');
grid on;
legend('Tx','Rx');

figure(4);
hold on;
plot(range,lp(ifft(S11 .* hanning(length(S11)))),'r')
hold off;
title(sprintf('Snow (velocity factor %.2f)',VF));
xlabel('range (m)')
ylabel('relative power (dB)');
grid on;
legend('Tx','Rx');

fn_dir = fileparts(fns{1});
out_fn = fullfile(fn_dir,'Kuband_freq_return_loss.fig');
saveas(1,out_fn);
out_fn = fullfile(fn_dir,'Kuband_freq_return_loss.jpg');
saveas(1,out_fn);
out_fn = fullfile(fn_dir,'Kuband_time_return_loss.fig');
saveas(2,out_fn);
out_fn = fullfile(fn_dir,'Kuband_time_return_loss.jpg');
saveas(2,out_fn);
out_fn = fullfile(fn_dir,'Snow_freq_return_loss.fig');
saveas(3,out_fn);
out_fn = fullfile(fn_dir,'Snow_freq_return_loss.jpg');
saveas(3,out_fn);
out_fn = fullfile(fn_dir,'Snow_time_return_loss.fig');
saveas(4,out_fn);
out_fn = fullfile(fn_dir,'Snow_time_return_loss.jpg');
saveas(4,out_fn);


return;
