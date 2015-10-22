
% base_path = 'D:\data\mcords\20120918\S11_antennas\Ground_S11_Measurements\PostInstall\';
base_path = 'D:\data\mcords\20120918\S11_antennas\Testflight_S11_Measurements';
colors = {'k','r','g','b','m'};

figure(1); clf;
fn = [];
for ant_chan = 1:5
%   fn{ant_chan} = fullfile(base_path,sprintf('S11_PostInstall_CH%i.s1p',ant_chan));
  fn{ant_chan} = fullfile(base_path,sprintf('S11_FlightTest1_CH%i.s1p',ant_chan));
  
  [freq, data, freq_noise, data_noise, Zo] = SXPParse(fn{ant_chan});

  S11 = squeeze(data(1,1,:));
  plot(freq/1e6, lp(S11), colors{ant_chan});
  hold on;
  xlabel('frequency (MHz)')
  ylabel('relative power (dB)');
  title_str = sprintf('channel %i S11', ant_chan);
  title(title_str);
  grid on;
  xlim([185 205])
  legends{ant_chan} = title_str;
end
legend(legends,'Location','SouthEast');
out_fn = fullfile(base_path,'S11_matlab.fig');
saveas(1,out_fn)
out_fn = fullfile(base_path,'S11_matlab.jpg');
saveas(1,out_fn)
