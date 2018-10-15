if 1
  % Plot analysis_mean
  param = read_param_xls(ct_filename_param('rds_param_2018_Antarctica_Ground.xls'),'20181014_02',{'analysis','analysis'});
  img = 2;
  wf_adcs = param.analysis.imgs{img};
  h_fig(1) = figure(1); clf(h_fig(1));
  h_axes(1) = axes('parent',h_fig(1));
  h_fig(2) = figure(2); clf(h_fig(2));
  h_axes(2) = axes('parent',h_fig(2));
  legend_str = {};
  for wf_adc = 1:length(wf_adcs)
    wf = wf_adcs(wf_adc,1);
    adc = wf_adcs(wf_adc,2);
    stat = load(fullfile(ct_filename_out(param,'analysis_mean','',1),sprintf('stats_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc)));
    plot(h_axes(1),stat.time{1}, 10*log10(stat.stats{1}{1}(:,1)/2 / 50) + 30 + stat.param_analysis.radar.wfs(wf).adc_gains_dB(adc));
    hold(h_axes(1),'on');
    plot(h_axes(2),stat.time{1}, 10*log10(stat.stats{1}{1}(:,1) / 50 * param.radar.wfs(wf).presums) + 30 + stat.param_analysis.radar.wfs(wf).adc_gains_dB(adc));
    hold(h_axes(2),'on');
    legend_str{end+1} = sprintf('%d-%d',wf,adc);
  end
  legend(h_axes(1),legend_str,'location','best');
  legend(h_axes(2),legend_str,'location','best');
  grid(h_axes(1),'on');
  grid(h_axes(2),'on');
  title(h_axes(1),'Coherent signals correct power');
  ylabel(h_axes(1),'Power (dBm)');
  xlabel(h_axes(1),'Time (sec)');
  title(h_axes(2),'Incoherent signals correct power');
  ylabel(h_axes(2),'Power (dBm)');
  xlabel(h_axes(2),'Time (sec)');
  
  return
end

if 0
  % Plot analysis_freq
  param = read_param_xls(ct_filename_param('rds_param_2018_Antarctica_Ground.xls'),'20181014_02',{'analysis','analysis'});
  img = 2;
  wf_adcs = param.analysis.imgs{img};
  num_ave = [1 10 100];
  h_fig(1) = figure(1); clf(h_fig(1));
  h_axes(1) = axes('parent',h_fig(1));
  h_fig(2) = figure(2); clf(h_fig(2));
  h_axes(2) = axes('parent',h_fig(2));
  legend_str = {};
  for idx = 1:length(num_ave)
    for wf_adc = 1:length(wf_adcs)
      wf = wf_adcs(wf_adc,1);
      adc = wf_adcs(wf_adc,2);
      stat = load(fullfile(ct_filename_out(param,'analysis_freq','',1),sprintf('stats_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc)));
      Nt = size(stat.stats{1}{1},1);
      dt = stat.time{1}(2) - stat.time{1}(1);
      fc = stat.freq{1}(1);
      T = Nt*dt;
      df = 1/T;
      freq = fftshift(fc + df * ifftshift(-floor(Nt/2) : floor((Nt-1)/2)).');
      if idx == 1
        h_plot(wf_adc) = plot(h_axes(1),freq, 10*log10(fftshift(stat.stats{1}{idx}(:,1))/2 / 50 / Nt) + 30 + stat.param_analysis.radar.wfs(wf).adc_gains_dB(adc));
        legend_str{end+1} = sprintf('%d-%d',wf,adc);
      else
        plot(h_axes(1),freq, 10*log10(fftshift(stat.stats{1}{idx}(:,1))/2 / 50 / Nt) + 30 + stat.param_analysis.radar.wfs(wf).adc_gains_dB(adc),'Color',get(h_plot(wf_adc),'Color'));
      end
      hold(h_axes(1),'on');
      plot(h_axes(2),freq, 10*log10(fftshift(stat.stats{1}{idx}(:,1)) / 50 * param.radar.wfs(wf).presums / Nt) + 30 + stat.param_analysis.radar.wfs(wf).adc_gains_dB(adc),'Color',get(h_plot(wf_adc),'Color'));
      hold(h_axes(2),'on');
    end
  end
  legend(h_axes(1),legend_str,'location','best');
  legend(h_axes(2),legend_str,'location','best');
  grid(h_axes(1),'on');
  grid(h_axes(2),'on');
  title(h_axes(1),'Coherent signals correct power');
  ylabel(h_axes(1),'Power (dBm)');
  xlabel(h_axes(1),'Frequency (Hz)');
  title(h_axes(2),'Incoherent signals correct power');
  ylabel(h_axes(2),'Power (dBm)');
  xlabel(h_axes(2),'Frequency (Hz)');
  
  return
end

if 0
  % Plot analysis_max
  param = read_param_xls(ct_filename_param('rds_param_2018_Antarctica_Ground.xls'),'20181014_02',{'analysis','analysis'});
  img = 2;
  wf_adcs = param.analysis.imgs{img};
  h_fig(1) = figure(1); clf(h_fig(1));
  h_axes(1) = axes('parent',h_fig(1));
  h_fig(2) = figure(2); clf(h_fig(2));
  h_axes(2) = axes('parent',h_fig(2));
  legend_str = {};
  for wf_adc = 1:length(wf_adcs)
    wf = wf_adcs(wf_adc,1);
    adc = wf_adcs(wf_adc,2);
    stat = load(fullfile(ct_filename_out(param,'analysis_max','',1),sprintf('stats_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc)));
    plot(h_axes(1),10*log10(abs(stat.stats{1}{1}(1,:)).^2/2 / 50) + 30 + stat.param_analysis.radar.wfs(wf).adc_gains_dB(adc));
    hold(h_axes(1),'on');
    plot(h_axes(2),10*log10(abs(stat.stats{1}{1}(1,:)).^2 / 50 * param.radar.wfs(wf).presums) + 30 + stat.param_analysis.radar.wfs(wf).adc_gains_dB(adc));
    hold(h_axes(2),'on');
    legend_str{end+1} = sprintf('%d-%d',wf,adc);
  end
  legend(h_axes(1),legend_str,'location','best');
  legend(h_axes(2),legend_str,'location','best');
  grid(h_axes(1),'on');
  grid(h_axes(2),'on');
  title(h_axes(1),'Coherent signals correct power');
  ylabel(h_axes(1),'Power (dBm)');
  xlabel(h_axes(1),'Time (sec)');
  title(h_axes(2),'Incoherent signals correct power');
  ylabel(h_axes(2),'Power (dBm)');
  xlabel(h_axes(2),'Time (sec)');
  
  return
end

if 0
  % Plot analysis_kx
  param = read_param_xls(ct_filename_param('rds_param_2018_Antarctica_Ground.xls'),'20181014_02',{'analysis','analysis'});
  img = 2;
  wf_adcs = param.analysis.imgs{img};
  h_fig = []; h_axes = [];
  for wf_adc = 1:size(wf_adcs,1)
    wf = wf_adcs(wf_adc,1);
    adc = wf_adcs(wf_adc,2);
    stat = load(fullfile(ct_filename_out(param,'analysis_kx','',1),sprintf('stats_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc)));
    h_fig(wf_adc) = figure(wf_adc); clf(h_fig(wf_adc));
    set(h_fig(wf_adc),'WindowStyle','docked','NumberTitle','off','Name',sprintf('%d-%d',wf,adc));
    if wf_adc == 1
      h_axes = axes('parent',h_fig(wf_adc));
    else
      h_axes(wf_adc) = axes('parent',h_fig(wf_adc));
    end
    imagesc(lp(stat.stats{1}{1}),'parent',h_axes(wf_adc));
    grid(h_axes(wf_adc),'on');
    title(h_axes(wf_adc),'Alongtrack wavenumber');
    ylabel(h_axes(wf_adc),'Range bin (bin)');
    xlabel(h_axes(wf_adc),'Spatial frequency bin (bin)');
    h_color = colorbar(h_axes(wf_adc));
    set(get(h_color,'YLabel'),'String','Relative power (dB)');
  end
  linkaxes(h_axes);
  
  return
end
