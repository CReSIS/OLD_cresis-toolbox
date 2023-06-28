% Load Ku-band data
% freq = freq + 15e9;
% save('/cresis/snfs1/scratch1/paden/global_hawk/kuband_subband.mat','g_data', 'freq', 'time');
% dbup
% save('/cresis/snfs1/scratch1/paden/global_hawk/kuband_subband.mat','-append','out_records');
%
% Load Snow data
% freq = freq + 5e9;
% save('/cresis/snfs1/scratch1/paden/global_hawk/snow_subband.mat','g_data', 'freq', 'time');
% dbup
% save('/cresis/snfs1/scratch1/paden/global_hawk/snow_subband.mat','-append','out_records');
% load('/cresis/snfs1/scratch1/paden/global_hawk/snow_subband.mat');
% kuband = load('/cresis/snfs1/scratch1/paden/global_hawk/kuband_subband.mat');
load('S:/scratch1/paden/global_hawk/snow_subband.mat');
kuband = load('S:/scratch1/paden/global_hawk/kuband_subband.mat');

dd = fft(g_data);
%freq = freq + 5e9;

noise_bins = 10000:12000;
dd_nb = fft(g_data(noise_bins,:));
dt = time(2)-time(1);
Nt_nb = length(noise_bins);
T_nb = Nt_nb * dt;
df_nb = 1/T_nb;
freq_nb = 5e9 + df_nb * [-floor(Nt_nb/2) : floor((Nt_nb-1)/2)].';

% Polyfit to noise
x_axis = linspace(-1,1,size(dd_nb,1));
select_bins = find(freq_nb > 2.25e9 & freq_nb < 7.75e9);
x_axis_trunc = x_axis(select_bins);
p = polyfit(x_axis_trunc, (lp(mean(abs(dd_nb(select_bins,:)).^2,2))-16).', 4);
noise_fit_nb = polyval(p, x_axis).';

noise_fit = interp1(freq_nb,noise_fit_nb,freq);

figure(9); clf;
plot(freq_nb/1e9, lp(mean(abs(dd_nb).^2,2))-16 - noise_fit_nb);
xlabel('Frequency (GHz)');
ylim([16 40]-16);
hold on;
plot(freq/1e9, lp(mean(abs(dd).^2,2))-16 - noise_fit,'r');
xlabel('Frequency (GHz)');
ylim([16 40]-16);
xlim([2.25 7.75]);

figure(11); clf;
imagesc([], time*1e6, lp(fir_dec(abs(g_data).^2,5)));
ylabel('Time (us)');
colormap(1-gray(256))
ylim([3.08 3.16])
title(sprintf('%.1f-%.1f GHz\n', freq(1)/1e9, freq(end)/1e9));

figure(12); clf;
imagesc([], (kuband.time-2e-9)*1e6, lp(fir_dec(abs(kuband.g_data).^2,5)));
ylabel('Time (us)');
colormap(1-gray(256))
ylim([3.08 3.16])
title(sprintf('%.1f-%.1f GHz\n', kuband.freq(1)/1e9, kuband.freq(end)/1e9));
figure(11); aa = axis;
figure(12); axis(aa);
pos = get(11,'Position');
set(12,'Position',pos);

freq_bins_list = {[1001:6000], [9001:14000]};

for freq_bins_idx = 1:length(freq_bins_list)
  freq_bins = freq_bins_list{freq_bins_idx};
  
  Nt_sb = length(freq_bins);
  df_sb = freq(2)-freq(1);
  BW_sb = df_sb*Nt_sb;
  dt_sb = 1/BW_sb;
  freq_sb = freq(freq_bins);
  time_sb = time(1) + (0:Nt_sb-1)*dt_sb;
  
  dd_sb = ifft(dd(freq_bins,:));
  
  figure(freq_bins_idx); clf;
  imagesc([], time_sb*1e6, lp(fir_dec(abs(dd_sb).^2,5)));
  title(sprintf('%.1f-%.1f GHz\n', freq_sb(1)/1e9, freq_sb(end)/1e9));
  ylabel('Time (us)');
  colormap(1-gray(256))
  figure(11); aa = axis;
  figure(freq_bins_idx); axis(aa);
  pos = get(11,'Position');
  set(freq_bins_idx,'Position',pos);

end








