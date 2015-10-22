
figure(1); clf;

colors = {'k','r','g','c','b','k:','r:','g:'};
for ant = 1:8
  fn = sprintf('ANT%02d.s1p', ant);
  [freq, data, freq_noise, data_noise, Zo] = SXPParse(fn);
  S11 = squeeze(data(1,1,:));
  h(ant) = plot(freq/1e6, 20*log10(abs(S11)), colors{ant});
  h_label{ant} = sprintf('Ant %02d', ant);
  hold on;
end
hold off;
xlim([100 600]);
ylim([-25 0]);
xlabel('Frequency (MHz');
ylabel('S_1_1 (dB)');
grid on;
legend(h,h_label);
