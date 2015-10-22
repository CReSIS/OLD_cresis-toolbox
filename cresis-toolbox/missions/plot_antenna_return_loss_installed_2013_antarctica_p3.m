if 1
  %% MCORDS4
  figure(1); clf; f1_axis = axes;
  figure(2); clf; f2_axis = axes;
  
  colors = {'k','r','g','c','b','k:','r:','g:'};
  clear h h2 h_label;
  antennas = 5:11;
  for ant_idx = 1:length(antennas)
    ant = antennas(ant_idx);
    fn = sprintf('p:\\metadata\\2013_Antarctica_P3\\install\\Antenna\\rds\\element%d.s1p', ant);
    [freq, data, freq_noise, data_noise, Zo] = SXPParse(fn);
    freq_idxs = find(freq >= 180e6 & freq <= 210e6);
    time_gate_window_func = inline('tukeywin_trim(N,0.2)');
    freq_window_func = inline('tukeywin_trim(N,0.2)');
    S11 = squeeze(data(1,1,:));
    h(ant_idx) = plot(f1_axis, freq/1e6, 20*log10(abs(S11)), colors{mod(ant_idx-1,length(colors))+1});
    hold(f1_axis,'on');
    h_label{ant_idx} = sprintf('Ant %02d', ant);
    
    Nt = length(freq_idxs);
    df = freq(2)-freq(1);
    BW = Nt*df;
    dt = 1/BW;
    time = (0:Nt-1).' * dt;
    h2(ant_idx) = plot(f2_axis, time*3e8/2, 20*log10(abs(ifft( S11(freq_idxs) .* freq_window_func(length(freq_idxs)) ))), ...
      colors{mod(ant_idx-1,length(colors))+1});
    hold(f2_axis,'on');
  end
  hold(f1_axis,'off');
  xlim(f1_axis, [100 300]);
  ylim(f1_axis, [-30 0]);
  xlabel(f1_axis, 'Frequency (MHz)');
  ylabel(f1_axis, 'S_1_1 (dB)');
  grid(f1_axis,'on');
  legend(f1_axis, h, h_label,'Location','southwest');
  
  hold(f2_axis,'off');
  xlim(f2_axis, [0 120]);
  ylim(f2_axis, [-60 0]);
  xlabel(f2_axis, 'Range, er=1 (m)');
  ylabel(f2_axis, 'S_1_1 (dB)');
  grid(f2_axis,'on');
  legend(f2_axis, h2, h_label);
  
  saveas(1,'rds_antenna_return_loss_center.fig');
  saveas(2,'rds_antenna_return_loss_time_center.fig');
  saveas(1,'rds_antenna_return_loss_center.jpg');
  saveas(2,'rds_antenna_return_loss_time_center.jpg');
  
  figure(3); clf; f1_axis = axes;
  figure(4); clf; f2_axis = axes;
  
  colors = {'k','r','g','c','b','k:','r:','g:'};
  clear h h2 h_label;
  antennas = [1:4 12:15];
  for ant_idx = 1:length(antennas)
    ant = antennas(ant_idx);
    fn = sprintf('p:\\metadata\\2013_Antarctica_P3\\install\\Antenna\\rds\\element%d.s1p', ant);
    [freq, data, freq_noise, data_noise, Zo] = SXPParse(fn);
    freq_idxs = find(freq >= 180e6 & freq <= 210e6);
    time_gate_window_func = inline('tukeywin_trim(N,0.2)');
    freq_window_func = inline('tukeywin_trim(N,0.2)');
    S11 = squeeze(data(1,1,:));
    h(ant_idx) = plot(f1_axis, freq/1e6, 20*log10(abs(S11)), colors{mod(ant_idx-1,length(colors))+1});
    hold(f1_axis,'on');
    h_label{ant_idx} = sprintf('Ant %02d', ant);
    
    Nt = length(freq_idxs);
    df = freq(2)-freq(1);
    BW = Nt*df;
    dt = 1/BW;
    time = (0:Nt-1).' * dt;
    h2(ant_idx) = plot(f2_axis, time*3e8/2, 20*log10(abs(ifft( S11(freq_idxs) .* freq_window_func(length(freq_idxs)) ))), ...
      colors{mod(ant_idx-1,length(colors))+1});
    hold(f2_axis,'on');
  end
  hold(f1_axis,'off');
  xlim(f1_axis, [100 300]);
  ylim(f1_axis, [-30 0]);
  xlabel(f1_axis, 'Frequency (MHz)');
  ylabel(f1_axis, 'S_1_1 (dB)');
  grid(f1_axis,'on');
  legend(f1_axis, h, h_label,'Location','southwest');
  
  hold(f2_axis,'off');
  xlim(f2_axis, [0 120]);
  ylim(f2_axis, [-60 0]);
  xlabel(f2_axis, 'Range, er=1 (m)');
  ylabel(f2_axis, 'S_1_1 (dB)');
  grid(f2_axis,'on');
  legend(f2_axis, h2, h_label);
  
  saveas(3,'rds_antenna_return_loss_center.fig');
  saveas(4,'rds_antenna_return_loss_time_center.fig');
  saveas(3,'rds_antenna_return_loss_center.jpg');
  saveas(4,'rds_antenna_return_loss_time_center.jpg');
  
end

%% Accumulation Radar
if 0
  
  S11 = [];
  for chan = 1:4
    % Only S11 is valid in these files, one file for each antenna channel
    fn = sprintf('p:\\metadata\\2013_Antarctica_P3\\install\\Antenna\\ACCUM-%d.s2p', chan);
    [freq, data, freq_noise, data_noise, Zo] = SXPParse(fn);
    S11{chan} = squeeze(data(1,1,:));
  end
  freq_idxs = find(freq >= 600e6 & freq <= 900e6);
  time_gate_window_func = inline('tukeywin_trim(N,0.2)');
  freq_window_func = inline('tukeywin_trim(N,0.2)');
  figure_offset = 30;
  title_string = 'Accumulation Radar Antennas';
    
  figure(figure_offset+1); clf;
  colors = {'k','r','g','b'};
  for chan = 1:4
    h_S11(chan) = plot(freq/1e6, 20*log10(abs(S11{chan})), colors{chan});
    legend_str{chan} = sprintf('Chan %d', chan);
    hold on;
  end
  hold off;
  xlabel('Frequency (MHz)');
  ylabel('S_1_1 Return Loss (dB)');
  grid on;
  title(title_string);
  legend(h_S11,legend_str);
  xlim([100 2000]);
  ylim([-30 0]);
  
  Nt = length(freq(freq_idxs));
  df = freq(2)-freq(1);
  BW = Nt*df;
  dt = 1/BW;
  time = (0:Nt-1).' * dt;
  
  figure(figure_offset+2); clf;
  for chan = 1:4
    h_S11(chan) = plot(time*3e8/2, 20*log10(abs(ifft(S11{chan}(freq_idxs)))), colors{chan});
    legend_str{chan} = sprintf('Chan %d', chan);
    hold on;
  end
  hold off;
  xlabel('Range, er=1 (m)');
  ylabel('Return Loss (dB)');
  grid on;
  title(title_string);
  legend(h_S11,legend_str);
  xlim([0 50]);
  ylim([-90 0]);
  
  saveas(figure_offset+1,'accum_antenna_return_loss.fig');
  saveas(figure_offset+2,'accum_antenna_return_loss_time.fig');
  
  saveas(figure_offset+1,'accum_antenna_return_loss.jpg');
  saveas(figure_offset+2,'accum_antenna_return_loss_time.jpg');
end


%% Kuband
if 0
  
  % Four measurement files are just repeats of themselves (to verify noise
  % statistics)
  fn = 'p:\metadata\2013_Antarctica_P3\install\Antenna\KU1.s2p';
  %   fn = 'p:\metadata\2013_Antarctica_P3\install\Antenna\KU2.s2p';
  %   fn = 'p:\metadata\2013_Antarctica_P3\install\Antenna\KU3.s2p';
  %   fn = 'p:\metadata\2013_Antarctica_P3\install\Antenna\KU4.s2p';
  [freq, data, freq_noise, data_noise, Zo] = SXPParse(fn);
  S11 = squeeze(data(1,1,:));
  S12 = squeeze(data(1,2,:));
  S21 = squeeze(data(2,1,:));
  S22 = squeeze(data(2,2,:));
  freq_idxs = find(freq >= 12e9 & freq <= 18e9);
  time_gate_window_func = @hanning;
  freq_window_func = @hanning;
  figure_offset = 10;
  title_string = 'Ku-band Altimeter Antennas';
  
  % Time gate antennas
  if 0
    % Run this to select the time gate bins
    figure(1+figure_offset); clf;
    A = ifft(S11(freq_idxs).* freq_window_func(length(freq_idxs))); plot(lp(A));
  end
  
  A = ifft(S11); B = zeros(size(A));
  time_gate_idxs = round(300*19/6):round(330*19/6);
  B(time_gate_idxs) = A(time_gate_idxs) .* time_gate_window_func(length(time_gate_idxs));
  %time_gate_idxs = 1:round(10*19/6);
  %B(time_gate_idxs) = A(time_gate_idxs) .* time_gate_window_func(length(time_gate_idxs));
  S11_tg = fft(B);
  
  A = ifft(S22); B = zeros(size(A));
  B(time_gate_idxs) = A(time_gate_idxs) .* time_gate_window_func(length(time_gate_idxs));
  S22_tg = fft(B);
  
  figure(figure_offset+1); clf;
  h_S11 = plot(freq/1e9, 20*log10(abs(S11)), 'r');
  hold on;
  h_S22 = plot(freq/1e9, 20*log10(abs(S22)), 'b');
  h_S11_tg = plot(freq/1e9, 20*log10(abs(S11_tg)), 'k','LineWidth',2);
  h_S22_tg = plot(freq/1e9, 20*log10(abs(S22_tg)), 'g','LineWidth',2);
  hold off;
  xlabel('Frequency (GHz)');
  ylabel('Return Loss (dB)');
  grid on;
  title(title_string);
  legend([h_S11 h_S22 h_S11_tg h_S22_tg],{'S11','S22','S11 Time Gate','S22 Time Gate'});
  xlim([11 19]);
  ylim([-30 0]);
  
  Nt = length(freq(freq_idxs));
  df = freq(2)-freq(1);
  BW = Nt*df;
  dt = 1/BW;
  time = (0:Nt-1).' * dt;
  
  figure(figure_offset+2); clf;
  h_S11 = plot(time*3e8/2, 20*log10(abs(ifft(S11(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'r');
  hold on;
  h_S22 = plot(time*3e8/2, 20*log10(abs(ifft(S22(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'b');
  h_S11_tg = plot(time*3e8/2, 20*log10(abs(ifft(S11_tg(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'k');
  h_S22_tg = plot(time*3e8/2, 20*log10(abs(ifft(S22_tg(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'g');
  hold off;
  xlabel('Range, er=1 (m)');
  ylabel('Return Loss (dB)');
  grid on;
  title(title_string);
  legend([h_S11 h_S22 h_S11_tg h_S22_tg],{'S11','S22','S11 Time Gate','S22 Time Gate'});
  xlim([0 18]);
  ylim([-90 0]);
  
  % Time gate antennas
  if 0
    % Run this to select the time gate bins
    figure(figure_offset+3); clf;
    A = ifft(S21(freq_idxs).* freq_window_func(length(freq_idxs))); plot(lp(A));
  end
  
  A = ifft(S21); B = zeros(size(A));
  time_gate_idxs = 1:length(B);
  B(time_gate_idxs) = A(time_gate_idxs) .* time_gate_window_func(length(time_gate_idxs));
  S21_tg = fft(B);
  
  A = ifft(S12); B = zeros(size(A));
  B(time_gate_idxs) = A(time_gate_idxs) .* time_gate_window_func(length(time_gate_idxs));
  S12_tg = fft(B);
  
  figure(figure_offset+3); clf;
  h_S12 = plot(freq/1e9, 20*log10(abs(S12)), 'r');
  hold on;
  h_S21 = plot(freq/1e9, 20*log10(abs(S21)), 'b');
  h_S12_tg = plot(freq/1e9, 20*log10(abs(S12_tg)), 'k','LineWidth',2);
  h_S21_tg = plot(freq/1e9, 20*log10(abs(S21_tg)), 'g','LineWidth',2);
  hold off;
  xlabel('Frequency (GHz)');
  ylabel('Isolation (dB)');
  grid on;
  title(title_string);
  legend([h_S12 h_S21 h_S12_tg h_S21_tg],{'S12','S21','S12 Time Gate','S21 Time Gate'});
  
  Nt = length(freq(freq_idxs));
  df = freq(2)-freq(1);
  BW = Nt*df;
  dt = 1/BW;
  time = (0:Nt-1).' * dt;
  
  figure(figure_offset+4); clf;
  h_S12 = plot(time*3e8/2, 20*log10(abs(ifft(S12(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'r');
  hold on;
  h_S21 = plot(time*3e8/2, 20*log10(abs(ifft(S21(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'b');
  h_S12_tg = plot(time*3e8/2, 20*log10(abs(ifft(S12_tg(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'k');
  h_S21_tg = plot(time*3e8/2, 20*log10(abs(ifft(S21_tg(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'g');
  hold off;
  xlabel('Range, er=1 (m)');
  ylabel('Isolation (dB)');
  grid on;
  title(title_string);
  legend([h_S12 h_S21 h_S12_tg h_S21_tg],{'S12','S21','S12 Time Gate','S21 Time Gate'});
  
  saveas(figure_offset+1,'kuband_antenna_return_loss.fig');
  saveas(figure_offset+2,'kuband_antenna_return_loss_time.fig');
  saveas(figure_offset+3,'kuband_antenna_coupling.fig');
  saveas(figure_offset+4,'kuband_antenna_coupling_time.fig');
  
  saveas(figure_offset+1,'kuband_antenna_return_loss.jpg');
  saveas(figure_offset+2,'kuband_antenna_return_loss_time.jpg');
  saveas(figure_offset+3,'kuband_antenna_coupling.jpg');
  saveas(figure_offset+4,'kuband_antenna_coupling_time.jpg');
end

%% Snow
if 0
  
  % Four measurement files are just repeats of themselves (to verify noise
  % statistics)
  fn = 'p:\metadata\2013_Antarctica_P3\install\Antenna\SNOW1.s2p';
  %   fn = 'p:\metadata\2013_Antarctica_P3\install\Antenna\SNOW2.s2p';
  %   fn = 'p:\metadata\2013_Antarctica_P3\install\Antenna\SNOW3.s2p';
  %   fn = 'p:\metadata\2013_Antarctica_P3\install\Antenna\SNOW4.s2p';
  [freq, data, freq_noise, data_noise, Zo] = SXPParse(fn);
  S11 = squeeze(data(1,1,:));
  S12 = squeeze(data(1,2,:));
  S21 = squeeze(data(2,1,:));
  S22 = squeeze(data(2,2,:));
  freq_idxs = find(freq >= 2e9 & freq <= 8e9);
  time_gate_window_func = inline('tukeywin_trim(N,0.2)');
  freq_window_func = @hanning;
  figure_offset = 20;
  title_string = 'Snow Radar Antennas';
  
  % Time gate antennas
  if 0
    % Run this to select the time gate bins
    figure(1+figure_offset); clf;
    A = ifft(S11(freq_idxs).* freq_window_func(length(freq_idxs))); plot(lp(A));
  end
  
  A = ifft(S11); B = zeros(size(A));
  time_gate_idxs = round(280*19/6):round(315*19/6);
  B(time_gate_idxs) = A(time_gate_idxs) .* time_gate_window_func(length(time_gate_idxs));
  %time_gate_idxs = 1:round(10*19/6);
  %B(time_gate_idxs) = A(time_gate_idxs) .* time_gate_window_func(length(time_gate_idxs));
  S11_tg = fft(B);
  
  A = ifft(S22); B = zeros(size(A));
  B(time_gate_idxs) = A(time_gate_idxs) .* time_gate_window_func(length(time_gate_idxs));
  S22_tg = fft(B);
  
  figure(figure_offset+1); clf;
  h_S11 = plot(freq/1e9, 20*log10(abs(S11)), 'r');
  hold on;
  h_S22 = plot(freq/1e9, 20*log10(abs(S22)), 'b');
  h_S11_tg = plot(freq/1e9, 20*log10(abs(S11_tg)), 'k','LineWidth',2);
  h_S22_tg = plot(freq/1e9, 20*log10(abs(S22_tg)), 'g','LineWidth',2);
  hold off;
  xlabel('Frequency (GHz)');
  ylabel('Return Loss (dB)');
  grid on;
  title(title_string);
  legend([h_S11 h_S22 h_S11_tg h_S22_tg],{'S11','S22','S11 Time Gate','S22 Time Gate'});
  xlim([1 19]);
  ylim([-30 0]);
  
  Nt = length(freq(freq_idxs));
  df = freq(2)-freq(1);
  BW = Nt*df;
  dt = 1/BW;
  time = (0:Nt-1).' * dt;
  
  figure(figure_offset+2); clf;
  h_S11 = plot(time*3e8/2, 20*log10(abs(ifft(S11(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'r');
  hold on;
  h_S22 = plot(time*3e8/2, 20*log10(abs(ifft(S22(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'b');
  h_S11_tg = plot(time*3e8/2, 20*log10(abs(ifft(S11_tg(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'k');
  h_S22_tg = plot(time*3e8/2, 20*log10(abs(ifft(S22_tg(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'g');
  hold off;
  xlabel('Range, er=1 (m)');
  ylabel('Return Loss (dB)');
  grid on;
  title(title_string);
  legend([h_S11 h_S22 h_S11_tg h_S22_tg],{'S11','S22','S11 Time Gate','S22 Time Gate'});
  xlim([0 18]);
  ylim([-90 0]);
  
  % Time gate antennas
  if 0
    % Run this to select the time gate bins
    figure(figure_offset+3); clf;
    A = ifft(S21(freq_idxs).* freq_window_func(length(freq_idxs))); plot(lp(A));
  end
  
  A = ifft(S21); B = zeros(size(A));
  time_gate_idxs = round(280*19/6):round(400*19/6);
  B(time_gate_idxs) = A(time_gate_idxs) .* time_gate_window_func(length(time_gate_idxs));
  S21_tg = fft(B);
  
  A = ifft(S12); B = zeros(size(A));
  B(time_gate_idxs) = A(time_gate_idxs) .* time_gate_window_func(length(time_gate_idxs));
  S12_tg = fft(B);
  
  figure(figure_offset+3); clf;
  h_S12 = plot(freq/1e9, 20*log10(abs(S12)), 'r');
  hold on;
  h_S21 = plot(freq/1e9, 20*log10(abs(S21)), 'b');
  h_S12_tg = plot(freq/1e9, 20*log10(abs(S12_tg)), 'k','LineWidth',2);
  h_S21_tg = plot(freq/1e9, 20*log10(abs(S21_tg)), 'g','LineWidth',2);
  hold off;
  xlabel('Frequency (GHz)');
  ylabel('Isolation (dB)');
  grid on;
  title(title_string);
  legend([h_S12 h_S21 h_S12_tg h_S21_tg],{'S12','S21','S12 Time Gate','S21 Time Gate'});
  
  Nt = length(freq(freq_idxs));
  df = freq(2)-freq(1);
  BW = Nt*df;
  dt = 1/BW;
  time = (0:Nt-1).' * dt;
  
  figure(figure_offset+4); clf;
  h_S12 = plot(time*3e8/2, 20*log10(abs(ifft(S12(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'r');
  hold on;
  h_S21 = plot(time*3e8/2, 20*log10(abs(ifft(S21(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'b');
  h_S12_tg = plot(time*3e8/2, 20*log10(abs(ifft(S12_tg(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'k');
  h_S21_tg = plot(time*3e8/2, 20*log10(abs(ifft(S21_tg(freq_idxs) .* freq_window_func(length(freq_idxs))))), 'g');
  hold off;
  xlabel('Range, er=1 (m)');
  ylabel('Isolation (dB)');
  grid on;
  title(title_string);
  legend([h_S12 h_S21 h_S12_tg h_S21_tg],{'S12','S21','S12 Time Gate','S21 Time Gate'});
  
  saveas(figure_offset+1,'snow_antenna_return_loss.fig');
  saveas(figure_offset+2,'snow_antenna_return_loss_time.fig');
  saveas(figure_offset+3,'snow_antenna_coupling.fig');
  saveas(figure_offset+4,'snow_antenna_coupling_time.fig');
  
  saveas(figure_offset+1,'snow_antenna_return_loss.jpg');
  saveas(figure_offset+2,'snow_antenna_return_loss_time.jpg');
  saveas(figure_offset+3,'snow_antenna_coupling.jpg');
  saveas(figure_offset+4,'snow_antenna_coupling_time.jpg');
end
