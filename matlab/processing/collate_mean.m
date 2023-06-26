% function collate_mean(param,param_override)
%
% Collects statistics results from run_analysis.m specified in analysis sheet
% Requires param.analysis.cmd{1} to be analysis_mean
% Generates plots to verify mean power levels for specified/all wf-adc pairs
% Magnitude (dBm/Hz) vs Range_line for individual plots for wf-adc pair
% Magnitude (dBm/Hz) vs adc for combined plot for all waveforms
%
% Example:
%   See run_collate_mean for how to run.
%
% Authors: John Paden, Hara Madhav Talasila

%% General Setup
% =========================================================================

param = merge_structs(param, param_override);
physical_constants;

fprintf('=============================================================\n');
fprintf('%s: %s (%s)\n', param.season_name, param.day_seg, datestr(now));
fprintf('=============================================================\n');

%% Input Checks
% =========================================================================
if ~isfield(param.analysis, 'plot_individual')
  param.analysis.plot_individual = 1;
end

if ~isfield(param.analysis, 'enable_visible_plot_individual')
  param.analysis.enable_visible_plot_individual = 0;
end

if ~isfield(param.analysis, 'plot_combined')
  param.analysis.plot_combined = 0;
end

if ~isfield(param.analysis, 'enable_visible_plot_combined')
  param.analysis.enable_visible_plot_combined = 0;
end

stat_dir = ct_filename_out(param,param.analysis.cmd{1}.out_path,'',1);
if ~exist(stat_dir, 'dir')
  fprintf('Empty directory --> run_analysis results (%s)\n',stat_dir);
  return;
end

out_dir = ct_filename_ct_tmp(param,'','stat','mean');
if ~exist(out_dir, 'dir')
  mkdir(out_dir);
end

%% Collate the mean reults
% =========================================================================

% For figure handles
if param.analysis.plot_individual
  fig_count = 0;
  for idx =1:length(param.analysis.imgs)
  fig_count = fig_count + size(param.analysis.imgs{idx},1);
  end
  h_fig = get_figures(fig_count,param.analysis.enable_visible_plot_individual);
end
fig_idx=0;

% The big loop for each img-wf_adc pair
for img = 1:length(param.analysis.imgs)
  pow_spec_cmd = [];
  for wf_adc = 1:size(param.analysis.imgs{img},1)
    wf = param.analysis.imgs{img}(wf_adc,1);
    adc = param.analysis.imgs{img}(wf_adc,2);
    fig_idx = fig_idx +1;
    %% Load the statistics file
    % =====================================================================
    dd = load(fullfile(stat_dir, sprintf('stats_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc)));
    % Initialization
    presums = dd.param_analysis.radar.wfs(wf).presums;
    rline_sum=0;
    rline_pow = [];
    num_blocks=0;
    
    % Collate
    for b_size = 1:length(dd.stats)
      rline_sum = rline_sum + sum(dd.stats{b_size}{1});
      rline_pow =[rline_pow dd.stats{b_size}{1}'];
      num_blocks = num_blocks + size(dd.stats{b_size}{1},1);
    end
    
    % Average the power spectrum over all the blocks
    rline_mean = rline_sum/num_blocks/50; %Check Z0=50 Ohm in param sheet 
    rline_pow = rline_pow/50;
    rline_pow_mean = sum(rline_pow)/num_blocks;

    % Calculate the expected Noise levels
    Noise_plot1 = 10*log10(BoltzmannConst*290*param.radar.noise_fig(adc)) + 30 - 10*log10(presums);
    Noise_plot2 = 10*log10(10.^(-9.0)/30e6) + 30 - 10*log10(presums) - dd.param_analysis.radar.wfs(wf).adc_gains_dB(adc);
    Noise_plot = 10*log10(1*10.^(Noise_plot1/10) + 1*10.^(Noise_plot2/10));
    Noise_plots = Noise_plot *ones(1,length(rline_pow));
    
    % Calculate mean values for Measured power vs Expected noise power
    mean_power(wf,adc) = 10*log10(rline_mean/30e6) + 30;
    mean_noise_plot(wf,adc) = mean(Noise_plot);
    
    %% Generate the plots
    % =====================================================================
    if param.analysis.plot_individual
      cur_fig = h_fig(fig_idx);
      clf(cur_fig);
      if param.analysis.enable_visible_plot_individual
        figure(cur_fig); % Brings the current figure to the top
      end
      h_axes = axes('parent',cur_fig);
      plot(h_axes, 10*log10(rline_pow/30e6)+30,'LineWidth',2);
      hold(h_axes, 'on');
      xlabel(h_axes,'Range lines');
      ylabel(h_axes,'Magnitude (dBm/Hz)');
      colors = get(h_axes,'ColorOrder');
      grid(h_axes, 'on');
      leg=[];
      for leg_idx = 1
        leg{leg_idx} = sprintf('( %.2f vs %.2f )',mean_power(wf,adc), mean_noise_plot(wf,adc)) ;
        plot(h_axes, Noise_plots,'Color',colors(leg_idx,:));
      end
      axis(h_axes, 'tight');
      hold(h_axes, 'off');
      title(h_axes,sprintf('For [ %d-%d ] Mean Power  (Measured[Thick] vs Expected[Line])',wf,adc));
      legend(h_axes,leg, 'Location', 'South');

      %% Save the individual plots
      saveas(cur_fig, fullfile(out_dir,sprintf('mean_wf_%d_adc_%d.fig',wf,adc)));
      saveas(cur_fig, fullfile(out_dir,sprintf('mean_wf_%d_adc_%d.jpg',wf,adc)));
    end
    % To see the progress of current segment
    if wf_adc == 1
      fprintf('wf-adc %d-%d',wf,adc);
    else
      fprintf(' %d-%d',wf,adc);
    end

  end
  fprintf('\n');
end

if param.analysis.plot_combined
  h_fig = get_figures(1,param.analysis.enable_visible_plot_combined);
  cur_fig = h_fig;
  clf(cur_fig);
  if param.analysis.enable_visible_plot_combined
    figure(cur_fig); % Brings the current figure to the top
  end
  h_axes = axes('parent',cur_fig);
  hold(h_axes, 'on');
  for img = 1:length(param.analysis.imgs)
    plot(h_axes, mean_power(img,:),'LineWidth',3);
  end
  xlabel(h_axes,sprintf('adc [1:%d]',adc));
  ylabel(h_axes,'Magnitude (dBm/Hz)');
  colors = get(h_axes,'ColorOrder');
  grid(h_axes, 'on');
  leg=[];
  for leg_idx = 1:length(param.analysis.imgs)
    leg{leg_idx} = sprintf('( %.2f vs %.2f ) wf-%d',mean(mean_power(leg_idx,:)), mean(mean_noise_plot(leg_idx,:)), leg_idx ) ;
    plot(h_axes, mean_noise_plot(leg_idx,:),'--','Color',colors(leg_idx,:));
  end
  axis(h_axes, 'tight');
  hold(h_axes, 'off');
  title(h_axes,sprintf('For %d imgs Mean Power (Measured[Thick] vs Expected[Thin])',img));
  legend(h_axes,leg, 'Location', 'Best');

  saveas(cur_fig, fullfile(out_dir,sprintf('mean_all.fig')));
  saveas(cur_fig, fullfile(out_dir,sprintf('mean_all.jpg')));
end

  
fprintf('=============================================================\n');
fprintf('%s: %s (%s)\n', param.season_name, param.day_seg, datestr(now));
fprintf('Output directory -->\n %s\n',out_dir);
fprintf('=============================================================\n');

