% function collate_fast_time(param,param_override)
%
% Collects statistics results from run_analysis.m specified in analysis sheet
% Requires the param.analysis.cmd{4} is analysis_fft
% Generates plots to verify noise levels for specified/all wf-adc pairs
% Magnitude (dBm/Hz) vs Frequency for each wf-adc pair
%
% Example:
%   See run_collate_fast_time for how to run.
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
if ~isfield(param.analysis, 'enable_visible_plot')
  param.analysis.enable_visible_plot = 0;
end

stat_dir = ct_filename_out(param,param.analysis.cmd{4}.out_path,'',1);
if ~exist(stat_dir, 'dir')
  fprintf('Empty directory --> run_analysis results (%s)\n',stat_dir);
  return;
end

out_dir = ct_filename_ct_tmp(param,'','stat','ft');
if ~exist(out_dir, 'dir')
  mkdir(out_dir);
end

%% Collate the noise reults
% =========================================================================

% For figure handles
fig_count = 0;
for idx =1:length(param.analysis.imgs)
  fig_count = fig_count + size(param.analysis.imgs{idx},1);
end
h_fig = get_figures(fig_count,param.analysis.enable_visible_plot);
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
    averages = dd.param_analysis.analysis.cmd{3}.ave;
    presums = dd.param_analysis.radar.wfs(wf).presums;
    pow_spec=[];
    fft_len = size(dd.stats{1}{1},1);
    fs = dd.param_analysis.radar.fs;
    num_blocks=[];
    
    % Collate
    for b_size = 1:length(dd.stats)
      num_blocks(b_size,:) = zeros(1,length(dd.stats{b_size}));
      for b_cmd = 1:length(dd.stats{b_size})
        pow_spec(:,b_cmd,b_size) = sum(dd.stats{b_size}{b_cmd}, 2);
        num_blocks(b_size,b_cmd) = num_blocks(b_size,b_cmd) + size(dd.stats{b_size}{b_cmd},2);
      end
    end
    
    % Average the power spectrum over all the blocks
    num_blocks = sum(num_blocks,1);
    pow_spec_cmd(:,:,adc) = bsxfun(@rdivide,sum(pow_spec,3),num_blocks);

%     % Calculate the expected Noise levels
    Noise_plot1 = bsxfun(@times, 10*log10(BoltzmannConst*290*param.radar.noise_fig(adc)./averages), ones(fft_len,1)) + 30 - 10*log10(presums);
    Noise_plot2 = bsxfun(@times, 10*log10(10.^(-9.0)/30e6./averages), ones(fft_len,1)) + 30 - 10*log10(presums) - dd.param_analysis.radar.wfs(wf).adc_gains_dB(adc);
    Noise_plot = 10*log10(1*10.^(Noise_plot1/10) + 1*10.^(Noise_plot2/10));
%     Noise_plot1 =  10*log10(BoltzmannConst*290*param.radar.noise_fig(adc)./averages) + 30 - 10*log10(presums);
%     Noise_plot2 =  10*log10(10.^(-9.0)/30e6./averages) + 30 - 10*log10(presums) - dd.param_analysis.radar.wfs(wf).adc_gains_dB(adc);
%     Noise_plot = 10*log10(1*10.^(Noise_plot1/10) + 1*10.^(Noise_plot2/10));
    
    df = dd.freq{1}(2)-dd.freq{1}(1);
    dt = dd.time{1}(2)-dd.time{1}(1);
    Nt = fft_len;
    df = 1/(Nt*dt);
    
    % Calculate mean values for Measured power vs Expected noise power
    mean_power = 10*log10(mean(pow_spec_cmd(2:end,:,adc),1)/30e6) + 30;
    mean_noise_plot = mean(Noise_plot);
    
    %% Generate the plots
    % =====================================================================
    cur_fig = h_fig(fig_idx);
    clf(cur_fig);
    if param.analysis.enable_visible_plot
      figure(cur_fig); % Brings the current figure to the top
    end
    h_axes = axes('parent',cur_fig);
    plot(h_axes, 10*log10(pow_spec_cmd(2:end,:,adc)/30e6)+30,'LineWidth',3);
    hold(h_axes, 'on');
    xlabel(h_axes,'Frequency');
    ylabel(h_axes,'Magnitude (dBm/Hz)');
    colors = get(h_axes,'ColorOrder');
    grid(h_axes, 'on');
    leg=[];
    for leg_idx = 1:length(averages)
      leg{leg_idx} = sprintf('( %.2f vs %.2f ) %d',mean_power(leg_idx), mean_noise_plot(leg_idx),averages(leg_idx)) ;
      plot(h_axes, Noise_plot(2:end,leg_idx),'Color',colors(leg_idx,:));
    end
    axis(h_axes, 'tight');
    hold(h_axes, 'off');
    title(h_axes,sprintf('For [ %d-%d ] Noise Power  (Measured[Thick] vs Expected[Line]) for ''n'' Averages',wf,adc));
    legend(h_axes,leg, 'Location', 'South');

    %% Save the plots
    saveas(cur_fig, fullfile(out_dir,sprintf('ft_wf_%d_adc_%d.fig',wf,adc)));
    saveas(cur_fig, fullfile(out_dir,sprintf('ft_wf_%d_adc_%d.jpg',wf,adc)));

    % To see the progress of current segment
    if wf_adc == 1
      fprintf('wf-adc %d-%d',wf,adc);
    else
      fprintf(' %d-%d',wf,adc);
    end

  end
  fprintf('\n');
end

fprintf('=============================================================\n');
fprintf('%s: %s (%s)\n', param.season_name, param.day_seg, datestr(now));
fprintf('Output directory -->\n %s\n',out_dir);
fprintf('=============================================================\n');

