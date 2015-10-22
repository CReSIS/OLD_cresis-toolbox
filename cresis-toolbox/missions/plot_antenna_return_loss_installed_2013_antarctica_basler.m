if 0
  %% MCORDS4
  figure(1); clf; f1_axis = axes;
  figure(2); clf; f2_axis = axes;
  
  colors = {'k','r','g','c','b','k:','r:','g:'};
  clear h h2 h_label;
  for ant = 1:8
    fn = sprintf('ANTENNA%d.s1p', ant);
    [freq, data, freq_noise, data_noise, Zo] = SXPParse(fn);
    S11 = squeeze(data(1,1,:));
    h(ant) = plot(f1_axis, freq/1e6, 20*log10(abs(S11)), colors{ant});
    hold(f1_axis,'on');
    h_label{ant} = sprintf('Ant %02d', ant);
    Nt = length(freq);
    df = freq(2)-freq(1);
    BW = Nt*df;
    dt = 1/BW;
    time = (0:Nt-1).' * dt;
    h2(ant) = plot(f2_axis, time*3e8/2, 20*log10(abs(ifft(S11))), colors{ant});
    hold(f2_axis,'on');
  end
  hold(f1_axis,'off');
  xlim(f1_axis, [150 525]);
  ylim(f1_axis, [-25 0]);
  xlabel(f1_axis, 'Frequency (MHz)');
  ylabel(f1_axis, 'S_1_1 (dB)');
  grid(f1_axis,'on');
  legend(f1_axis, h, h_label);
  
  hold(f2_axis,'off');
  xlim(f2_axis, [0 40]);
  ylim(f2_axis, [-60 0]);
  xlabel(f2_axis, 'Range, er=1 (m)');
  ylabel(f2_axis, 'S_1_1 (dB)');
  grid(f2_axis,'on');
  legend(f2_axis, h2, h_label);
  
  saveas(1,'mcords4_antenna_return_loss.fig');
  saveas(2,'mcords4_antenna_return_loss_time.fig');
  saveas(1,'mcords4_antenna_return_loss.jpg');
  saveas(2,'mcords4_antenna_return_loss_time.jpg');
  
end
%% Kuband
if 0
  
  figure(11); clf;
  figure(12); clf;
  figure(13); clf;
  figure(14); clf;
  
  fn = 'KU.s2p';
  [freq, data, freq_noise, data_noise, Zo] = SXPParse(fn);
  S11 = squeeze(data(1,1,:));
  S12 = squeeze(data(1,2,:));
  S21 = squeeze(data(2,1,:));
  S22 = squeeze(data(2,2,:));
  
  figure(11);
  h_S11 = plot(freq/1e9, 20*log10(abs(S11)), 'r');
  hold on;
  h_S22 = plot(freq/1e9, 20*log10(abs(S22)), 'b');
  hold off;
  xlabel('Frequency (GHz)');
  ylabel('Return Loss (dB)');
  grid on;
  legend([h_S11 h_S22],{'S11','S22'});
  xlim([11 19]);
  ylim([-55 0]);
  
  Nt = length(freq);
  df = freq(2)-freq(1);
  BW = Nt*df;
  dt = 1/BW;
  time = (0:Nt-1).' * dt;
  
  figure(12);
  h_S11 = plot(time*3e8/2, 20*log10(abs(ifft(S11))), 'r');
  hold on;
  h_S22 = plot(time*3e8/2, 20*log10(abs(ifft(S22))), 'b');
  hold off;
  xlabel('Range, er=1 (m)');
  ylabel('Return Loss (dB)');
  grid on;
  legend([h_S11 h_S22],{'S11','S22'});
  ylim([-90 0]);
  
  figure(13);
  h_S11 = plot(freq/1e9, 20*log10(abs(S12)), 'r');
  hold on;
  h_S22 = plot(freq/1e9, 20*log10(abs(S21)), 'b');
  hold off;
  xlabel('Frequency (GHz)');
  ylabel('Return Loss (dB)');
  grid on;
  legend([h_S11 h_S22],{'S12','S21'});
  
  Nt = length(freq);
  df = freq(2)-freq(1);
  BW = Nt*df;
  dt = 1/BW;
  time = (0:Nt-1).' * dt;
  
  figure(14);
  h_S11 = plot(time*3e8/2, 20*log10(abs(ifft(S12))), 'r');
  hold on;
  h_S22 = plot(time*3e8/2, 20*log10(abs(ifft(S21))), 'b');
  hold off;
  xlabel('Range, er=1 (m)');
  ylabel('Return Loss (dB)');
  grid on;
  legend([h_S11 h_S22],{'S12','S21'});
  
  saveas(11,'kuband_antenna_return_loss.fig');
  saveas(12,'kuband_antenna_return_loss_time.fig');
  saveas(13,'kuband_antenna_coupling.fig');
  saveas(14,'kuband_antenna_coupling_time.fig');
  
  saveas(11,'kuband_antenna_return_loss.jpg');
  saveas(12,'kuband_antenna_return_loss_time.jpg');
  saveas(13,'kuband_antenna_coupling.jpg');
  saveas(14,'kuband_antenna_coupling_time.jpg');
end

%% Snow

figure(21); clf;
figure(22); clf;
figure(23); clf;
figure(24); clf;

fn = 'SNOW.s2p';
[freq, data, freq_noise, data_noise, Zo] = SXPParse(fn);
S11 = squeeze(data(1,1,:));
S12 = squeeze(data(1,2,:));
S21 = squeeze(data(2,1,:));
S22 = squeeze(data(2,2,:));

figure(21);
h_S11 = plot(freq/1e9, 20*log10(abs(S11)), 'r');
hold on;
h_S22 = plot(freq/1e9, 20*log10(abs(S22)), 'b');
hold off;
xlabel('Frequency (GHz)');
ylabel('Return Loss (dB)');
grid on;
legend([h_S11 h_S22],{'S11','S22'});
xlim([1 9]);
ylim([-55 0]);

Nt = length(freq);
df = freq(2)-freq(1);
BW = Nt*df;
dt = 1/BW;
time = (0:Nt-1).' * dt;

figure(22);
h_S11 = plot(time*3e8/2, 20*log10(abs(ifft(S11))), 'r');
hold on;
h_S22 = plot(time*3e8/2, 20*log10(abs(ifft(S22))), 'b');
hold off;
xlabel('Range, er=1 (m)');
ylabel('Return Loss (dB)');
grid on;
legend([h_S11 h_S22],{'S11','S22'});
ylim([-90 0]);

% Time gate antennas
A = ifft(S21); B = zeros(size(A)); B(447:467) = A(447:467); S21_tg = fft(B);

figure(23);
h_S11 = plot(freq/1e9, 20*log10(abs(S12)), 'r');
hold on;
h_S22 = plot(freq/1e9, 20*log10(abs(S21)), 'b');
h_S22_tg = plot(freq/1e9, 20*log10(abs(S21_tg)), 'k');
hold off;
xlabel('Frequency (GHz)');
ylabel('Return Loss (dB)');
grid on;
legend([h_S11 h_S22 h_S22_tg],{'S12','S21','S21 Time Gate'});

Nt = length(freq);
df = freq(2)-freq(1);
BW = Nt*df;
dt = 1/BW;
time = (0:Nt-1).' * dt;

figure(24);
h_S11 = plot(time*3e8/2, 20*log10(abs(ifft(S12))), 'r');
hold on;
h_S22 = plot(time*3e8/2, 20*log10(abs(ifft(S21))), 'b');
hold off;
xlabel('Range, er=1 (m)');
ylabel('Return Loss (dB)');
grid on;
legend([h_S11 h_S22],{'S12','S21'});

saveas(21,'snow_antenna_return_loss.fig');
saveas(22,'snow_antenna_return_loss_time.fig');
saveas(23,'snow_antenna_coupling.fig');
saveas(24,'snow_antenna_coupling_time.fig');

saveas(21,'snow_antenna_return_loss.jpg');
saveas(22,'snow_antenna_return_loss_time.jpg');
saveas(23,'snow_antenna_coupling.jpg');
saveas(24,'snow_antenna_coupling_time.jpg');

return





