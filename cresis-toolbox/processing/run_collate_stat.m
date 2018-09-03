% figure(1); clf;
% title('Low Gain Not Synced 555 kHz PA off'); hold on;
% fn = '/scratch/accum/2018_Antarctica_TObas/CSARP_analysis_kx/stats_20180831_01_wf_1_adc_1.mat';

figure(17); hold on;
title('High Gain Alis unsynched vs. Haras Synched'); hold on;
fn = '/scratch/accum/2018_Antarctica_TObas/CSARP_analysis_kx/stats_20180831_15_wf_2_adc_1.mat';

% figure(5); clf;
% title('Low Gain ~Synced 555 kHz PA on'); hold on;
% fn = '/scratch/accum/2018_Antarctica_TObas/CSARP_analysis_kx/stats_20180831_02_wf_1_adc_1.mat';

% figure(6); hold on; clf;
% title('High Gain ~Synced 555 kHz PA on'); hold on;
% fn = '/scratch/accum/2018_Antarctica_TObas/CSARP_analysis_kx/stats_20180831_02_wf_2_adc_1.mat';

% figure(1); clf;
% title('Low Gain Synced 555 kHz PA off'); hold on;
% fn = '/scratch/accum/2018_Antarctica_TObas/CSARP_analysis_kx/stats_20180831_03_wf_1_adc_1.mat';

% figure(2); clf;
% title('High Gain Synced 555 kHz PA off'); hold on;
% fn = '/scratch/accum/2018_Antarctica_TObas/CSARP_analysis_kx/stats_20180831_03_wf_2_adc_1.mat';

% figure(3); clf;
% title('Low Gain Not Synced 555 kHz PA off'); hold on;
% fn = '/scratch/accum/2018_Antarctica_TObas/CSARP_analysis_kx/stats_20180831_04_wf_1_adc_1.mat';

% figure(4); hold on; clf;
% title('High Gain Not Synced 555 kHz PA off'); hold on;
% fn = '/scratch/accum/2018_Antarctica_TObas/CSARP_analysis_kx/stats_20180831_04_wf_2_adc_1.mat';

% figure(3); clf;
% title('Low Gain Synced 500 kHz PA On'); hold on;
% fn = '/scratch/accum/2018_Antarctica_TObas/CSARP_analysis_kx/stats_20180831_08_wf_1_adc_1.mat';

% figure(4); clf;
% title('High Gain Synced 500 kHz PA On'); hold on;
% fn = '/scratch/accum/2018_Antarctica_TObas/CSARP_analysis_kx/stats_20180831_08_wf_2_adc_1.mat';

% =========================================================================

stats = load(fn);

block = 1;

legend_str = {};
for stat_idx = 1:3

  fc = stats.freq{block}(1);
  dt = stats.time{block}(2)-stats.time{block}(1);
  Nt = size(stats.stats{block}{stat_idx},1);
  df = 1/(dt*Nt);
  freq = fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';

  legend_str{stat_idx} = sprintf('%d: %s\n',stat_idx,func2str(stats.param_analysis.analysis.cmd{3}.stats{stat_idx}));
  
  plot(fftshift(freq/1e6), fftshift(lp(stats.stats{block}{stat_idx}(:,1),1)));
  hold on;
  grid on
  xlabel('Frequency (MHz)');
  ylabel('Relative power (dB)');
  
end
legend(legend_str,'location','best','interpreter','none');
