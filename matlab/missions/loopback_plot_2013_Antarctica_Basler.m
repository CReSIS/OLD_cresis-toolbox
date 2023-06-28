base_path = 'C:\tmp\mcords4\ground_tests';

chan = 1;
fns = get_filenames(fullfile(base_path,sprintf('chan%d',chan)),'mcords4','','0001.bin');

% [8 8 4 7 3 6 5 1]
% The fns index to use for each channel is:
chan_to_fns = [8 NaN 5 3 7 6 4 1];

colors = {'k' 'r' 'g' 'b' 'k:' 'r:' 'g:' 'b:'};

figure(100); clf;
h_plot = [];
h_label = {};
for chan = 1:8
  fns = get_filenames(fullfile(base_path,sprintf('chan%d',chan)),'mcords4','','0001.bin');
  
  if isnan(chan_to_fns(chan))
    fprintf('No data for chan %d\n', chan);
    continue
  end
  
  fn = fns{chan_to_fns(chan)};
  
  [hdr,data] = basic_load_mcords4(fn,struct('clk',1e9/8));
  
  clear pc_param;
  signal = mean(data{5}(:,1:end) - j*data{6}(:,1:end),2);
  pc_param.Tpd = 10e-6;
%   signal = mean(data{3}(:,1:end) - j*data{4}(:,1:end),2);
%   pc_param.Tpd = 3e-6;
  
  pc_param.f0 = 200e6;
  pc_param.f1 = 450e6;
  fs = 1e9/2;
  Nt = length(signal);
  t0 = -10e-6;
  dt = 1/fs;
  BW = pc_param.f1-pc_param.f0;
  alpha = BW/pc_param.Tpd;
  pc_param.time = t0 + (0:dt:(Nt-1)*dt).';
  
  pc_param.tukey = 0.2;
  pc_param.rx_gain = 1;
  pc_param.window_func = @hanning;
  [pc_signal,pc_time] = pulse_compress(signal,pc_param);
  
  figure(chan); clf;
  plot(pc_time*1e6, lp(pc_signal,2));
  title(sprintf('Chan %d', chan));
  %   plot(pc_param.time*1e6, lp(signal,2),'r');
  xlabel('Time (us)');
  ylabel('Relative power (dB)');
  xlim([-7.7 -7.5]);
  grid on;
  
  figure(100);
  h_plot(chan) = plot(pc_time*1e6, lp(pc_signal,2), colors{chan});
  hold on;
  h_label{chan} = sprintf('Chan %d', chan);
  title(h_label{chan});
  %   plot(pc_param.time*1e6, lp(signal,2),'r');
  xlabel('Time (us)');
  ylabel('Relative power (dB)');
  xlim([-7.7 -7.5]);
  grid on;
end
hold off;

h_plot(isnan(chan_to_fns)) = [];
h_label(isnan(chan_to_fns)) = [];
legend(h_plot,h_label);



return



% OLD CODE TO FIGURE OUT DIRECTORY MAPPING
base_path = 'C:\tmp\mcords4';

chan = 1;
fns = get_filenames(fullfile(base_path,sprintf('chan%d',chan)),'mcords4','','0001.bin');

% [8 8 4 7 3 6 5 1]
% The fns index to use for each channel is:
[8 NaN 5 3 7 6 4 1]

for fn_idx = 8:length(fns)
  
  for chan = 1:8
    fns = get_filenames(fullfile(base_path,sprintf('chan%d',chan)),'mcords4','','0001.bin');
    
    fn = fns{fn_idx};
    
    [hdr,data] = basic_load_mcords4(fn,struct('clk',1e9/8));
    
    signal = mean(data{5}(:,1:1000) - j*data{6}(:,1:1000),2);
    
    clear pc_param;
    pc_param.f0 = 200e6;
    pc_param.f1 = 450e6;
    pc_param.Tpd = 10e-6;
    fs = 1e9/2;
    Nt = length(signal);
    t0 = -10e-6;
    dt = 1/fs;
    BW = pc_param.f1-pc_param.f0;
    alpha = BW/pc_param.Tpd;
    pc_param.time = t0 + (0:dt:(Nt-1)*dt).';
    
    pc_param.tukey = 0.2;
    pc_param.rx_gain = 1;
    [pc_signal,pc_time] = pulse_compress(signal,pc_param);
    
    plot(pc_param.time*1e6, lp(signal,2),'r');
    title(sprintf('Chan %d', chan));
    hold on;
    plot(pc_time*1e6, lp(pc_signal,2));
    hold off;
    xlabel('Time (us)');
    ylabel('Relative power (dB)');
    pause;
  end
  
end

