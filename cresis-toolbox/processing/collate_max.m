% function collate_max(param,param_override)
%
% Collects statistics results from run_analysis.m specified in analysis sheet
% Requires the param.analysis.cmd{2} is analysis_max
% Generates plots to verify max power levels for specified/all wf-adc pairs
% Normalized distribution vs Magnitude dBm for individual plots for each wf-adc pair
% Normalized distribution vs Magnitude dBm for combined plot for each waveform
%
% Example:
%   See run_collate_max for how to run.
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

if ~isfield(param.analysis, 'BinWidth')
  param.analysis.BinWidth = 1; %dB
end

stat_dir = ct_filename_out(param,param.analysis.cmd{2}.out_path,'',1);
if ~exist(stat_dir, 'dir')
  fprintf('Empty directory --> run_analysis results (%s)\n',stat_dir);
  return;
end

out_dir = ct_filename_ct_tmp(param,'','stat','max');
if ~exist(out_dir, 'dir')
  mkdir(out_dir);
end

%% Collate the max reults
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
max_wf = [];
max_adc = [];

% The big loop for each img-wf_adc pair
for img = 1:length(param.analysis.imgs)
  for wf_adc = 1:size(param.analysis.imgs{img},1)
    wf = param.analysis.imgs{img}(wf_adc,1);
    adc = param.analysis.imgs{img}(wf_adc,2);
    max_wf = max([wf,max_wf]);
    max_adc = max([adc,max_adc]);
    fig_idx = fig_idx +1;
    %% Load the statistics file
    % =====================================================================
    dd = load(fullfile(stat_dir, sprintf('stats_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc)));
    % Initialization
    rline_pow = [];
    num_blocks=0;
    
    % Collate
    for b_size = 1:length(dd.stats)
      rline_pow = [rline_pow dd.stats{b_size}{1}(1,:)];
      num_blocks = num_blocks + size(dd.stats{b_size}{1}(1,:),2);
    end
    
    % Calculate the Max power for all the blocks
    rline_max = rline_pow.*conj(rline_pow) / 50; %Check Z0=50 Ohm in param sheet 
    rline_max_dB(wf,adc,:) = 10*log10(abs(rline_max)) + 30; % in dBm  
    temp = squeeze(rline_max_dB(wf,adc,:));
    temp2 = temp(~isinf(temp));
    rline_max_dB_mean(wf,adc) = mean(temp2,1,'omitnan');
    rline_max_dB_std(wf,adc) = std(temp2,1,'omitnan');
    %% Generate the plots
    % =====================================================================
    if param.analysis.plot_individual
      cur_fig = h_fig(fig_idx);
      clf(cur_fig);
      if param.analysis.enable_visible_plot_individual
        figure(cur_fig); % Brings the current figure to the top
      end
      h_axes = axes('parent',cur_fig);
      histogram(h_axes, rline_max_dB(wf,adc,:),...
        'BinWidth', param.analysis.BinWidth,...
        'Normalization', 'probability','EdgeColor', 'none');
      hold(h_axes, 'on');
      xlabel(h_axes,'Magnitude, dBm');
      ylabel(h_axes,'Normalized distribution');
      colors = get(h_axes,'ColorOrder');
      grid(h_axes, 'on');
      leg=[];
      leg = sprintf('rlines = %d\n \\mu = %.2f dBm\n \\sigma = %.2f dBm',...
        num_blocks,rline_max_dB_mean(wf,adc),rline_max_dB_std(wf,adc)) ;
%       axis(h_axes, 'tight');
      hold(h_axes, 'off');
      title(h_axes,sprintf('%s : For [ %d-%d ] Max Power',...
        param.day_seg,wf,adc), 'Interpreter', 'none');
      legend(h_axes,leg, 'Location', 'NorthEast','FontSize',12,'FontWeight', 'bold');

      %% Save the individual plots
      saveas(cur_fig, fullfile(out_dir,sprintf('max_%s_wf_%d_adc_%d.fig',param.day_seg,wf,adc)));
      saveas(cur_fig, fullfile(out_dir,sprintf('max_%s_wf_%d_adc_%d.jpg',param.day_seg,wf,adc)));
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

%Table Output
fprintf('=============================================================\n');
fprintf('Analysis_MAX for wf-adc pairs: Mean (Standard_deviation) dBm\n');
table_rows = max_wf;
table_cols = max_adc;
table = [];
fprintf('| mu(std)\t');
for col_idx = 1:table_cols
    fprintf('| \tadc-%d\t\t',col_idx);
end
fprintf('|\n');
for row_idx = 1:table_rows
  for col_idx = 1:table_cols
    if col_idx == 1
      fprintf('| wf-%d\t\t',row_idx);
    end
    if ~rline_max_dB_mean(row_idx,col_idx)==0
      fprintf('| %.2f(%.2f)\t',rline_max_dB_mean(row_idx,col_idx),rline_max_dB_std(row_idx,col_idx));
    else
      fprintf('| \t\t\t\t');
    end
  end
  fprintf('|\n');
end


% For combined plot (Not used generally. Difficult to distinguish)
if param.analysis.plot_combined
  h_fig = get_figures(1,param.analysis.enable_visible_plot_combined);
  cur_fig = h_fig;
  
  for img = 1:length(param.analysis.imgs)
    
    clf(cur_fig);
    if param.analysis.enable_visible_plot_combined
      figure(cur_fig); % Brings the current figure to the top
    end
    h_axes = axes('parent',cur_fig);
    hold(h_axes, 'on');
    leg=[];
    leg_idx=0;
    
    for wf_adc = 1:size(param.analysis.imgs{img},1)
      wf = param.analysis.imgs{img}(wf_adc,1);
      adc = param.analysis.imgs{img}(wf_adc,2);
      histogram(h_axes, rline_max_dB(wf,adc,:),...
        'BinWidth', param.analysis.BinWidth,...
        'Normalization', 'probability','EdgeColor', 'none');
      leg_idx=leg_idx+1;
      leg{leg_idx} = sprintf('%d-%d',wf,adc ) ;
    end
    
    xlabel(h_axes,'Magnitude, dBm');
    ylabel(h_axes,sprintf('Normalized distribution'));
    colors = get(h_axes,'ColorOrder');
    grid(h_axes, 'on');
%     axis(h_axes, 'tight');
    legend(h_axes,leg, 'Location', 'NorthEast');
    title(h_axes,sprintf('For img - %d Max Power',wf));
    saveas(cur_fig, fullfile(out_dir,sprintf('max_%s_all_wf_%d.fig',param.day_seg,wf)));
    saveas(cur_fig, fullfile(out_dir,sprintf('max_%s_all_wf_%d.jpg',param.day_seg,wf)));
    
  end
end

  
fprintf('=============================================================\n');
fprintf('%s: %s (%s)\n', param.season_name, param.day_seg, datestr(now));
fprintf('Output directory -->\n %s\n',out_dir);
fprintf('=============================================================\n');

