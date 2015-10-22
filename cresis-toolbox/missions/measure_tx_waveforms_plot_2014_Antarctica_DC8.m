base_dir = 'D:\tx_Waveforms\';

plot_TTL_control_signals_en = false;
plot_WG_en = true;
plot_TxAnt_en = false;

TTL_control_signal_channels = [1 3];
signal_channel = 2;

waveforms = {'10000ns','3000ns','1000ns','30ns'};
xlims = {[8.8 19], [8.8 12], [8.8 10], [8.82 8.87]}; % WG
% xlims = {[8.8 19], [8.8 12], [8.8 10], [8.86 8.92]}; % TxAnt

if 0
  % Plot WG output and Tx Antenna output
  for chan = 1:6
    
    for wf_idx = 1:length(waveforms)
      figure(chan*10 + wf_idx); clf;
      if plot_TTL_control_signals_en
        hold on;
        fn = fullfile(base_dir,sprintf('chan%d_%s.mat',chan,waveforms{wf_idx}));
        if exist(fn,'file')
          for idx = TTL_control_signal_channels
            plot(time*1e6,data(:,idx),'m');
          end
        end
      end
      if plot_TxAnt_en
        hold on;
        fn = fullfile(base_dir,sprintf('txamp_chan%d_%s.mat',chan,waveforms{wf_idx}));
        if exist(fn,'file')
          load(fn,'time','data','notes');
          plot(time*1e6,data(:,signal_channel),'r');
        end
      end
      if plot_WG_en
        hold on;
        fn = fullfile(base_dir,sprintf('chan%d_%s.mat',chan,waveforms{wf_idx}));
        if exist(fn,'file')
          load(fn,'time','data','notes');
          plot(time*1e6,data(:,signal_channel),'b');
        end
      end
      hold off;
      %title(notes);
      title(waveforms{wf_idx});
      grid on;
      xlabel('Time (us)')'
      ylabel('Voltage (V)')'
      set(chan*10 + wf_idx, 'WindowStyle', 'docked');
      xlim(xlims{wf_idx});
    end
    
  end
end

if 1
  % Plot all WG outputs on top of each other
  chan_to_plot = [1:6];
  for wf_idx = 1:length(waveforms)
    figure(1 + wf_idx); clf;
    
    h_plot = [];
    legend_str = {};
    for chan = chan_to_plot
      
      if plot_TTL_control_signals_en
        hold on;
        fn = fullfile(base_dir,sprintf('chan%d_%s.mat',chan,waveforms{wf_idx}));
        if exist(fn,'file')
          for idx = TTL_control_signal_channels
            plot(time*1e6,data(:,idx),'m');
          end
        end
      end
      if plot_TxAnt_en
        hold on;
        fn = fullfile(base_dir,sprintf('txamp_chan%d_%s.mat',chan,waveforms{wf_idx}));
        if exist(fn,'file')
          load(fn,'time','data','notes');
          h_plot(end+1) = plot(time*1e6,data(:,signal_channel),'r');
          plot_color(chan,h_plot(end));
          legend_str{end+1} = sprintf('Chan %d', chan);
        end
      end
      if plot_WG_en
        hold on;
        fn = fullfile(base_dir,sprintf('chan%d_%s.mat',chan,waveforms{wf_idx}));
        if exist(fn,'file')
          load(fn,'time','data','notes');
          h_plot(end+1) = plot(time*1e6,data(:,signal_channel),'b');
          plot_color(chan,h_plot(end));
          legend_str{end+1} = sprintf('Chan %d', chan);
        end
      end
      hold off;
      %title(notes);
      title(waveforms{wf_idx});
      grid on;
      xlabel('Time (us)');
      ylabel('Voltage (V)');
      set(1 + wf_idx, 'WindowStyle', 'docked');
      xlim(xlims{wf_idx});
    end
    legend(h_plot, legend_str);
    
  end
end




