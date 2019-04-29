% function collate_gain(param,param_override)
%
% Collects waveform results from run_analysis.m specified in analysis sheet
% Requires the param.analysis.cmd{1} is waveform
% Generates plots to verify fast-time gain for specified/all wf-adc pairs
% Normalized gain (dB) vs relative time for wf-adc pair
%
% Example:
%   See run_collate_gain for how to run.
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

if ~isfield(param.analysis, 'raw_plot_en')
  param.analysis.raw_plot_en = 0;
end

if ~isfield(param.analysis, 'ftg_plot_en')
  param.analysis.ftg_plot_en = 1;
end

wf_dir = ct_filename_out(param,param.analysis.cmd{1}.out_path,'',1);
if ~exist(wf_dir, 'dir')
  fprintf('Empty directory --> run_analysis results (%s)\n',wf_dir);
  return;
end

out_dir = ct_filename_ct_tmp(param,'','waveform','gain');
if ~exist(out_dir, 'dir')
  mkdir(out_dir);
end

%% Collate the ftg reults
% =========================================================================

% For figure handles
fig_count = 0;
leg=[];
for idx =1:length(param.analysis.imgs)
  fig_count = fig_count + size(param.analysis.imgs{idx},1);
end
h_fig = get_figures(fig_count+1,param.analysis.enable_visible_plot);
fig_idx=0;
comb_plot = fig_count+1;
clf(h_fig(comb_plot));

if param.analysis.raw_plot_en
  fig2_count = 0;
  leg2=[];
  for idx =1:length(param.analysis.imgs)
    fig2_count = fig2_count + size(param.analysis.imgs{idx},1);
  end
  h_fig2 = get_figures(fig2_count+1,param.analysis.enable_visible_plot);
  fig2_idx=0;
  comb_plot2 = fig2_count+1;
  clf(h_fig2(comb_plot2));
end

% The big loop for each img-wf_adc pair
for img = 1:length(param.analysis.imgs)
  for wf_adc = 1:size(param.analysis.imgs{img},1)
    wf = param.analysis.imgs{img}(wf_adc,1);
    adc = param.analysis.imgs{img}(wf_adc,2);
    fig_idx = fig_idx +1;
    %% Load the waveform file
    % =====================================================================
    dd = load(fullfile(wf_dir, sprintf('waveform_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc)));
    % Initialization
    Nt = size(dd.wf_data,1);
    r_lines = size(dd.wf_data,2);
%     if size(dd.time_rng,1) == 2 && size(dd.time_rng,2) == r_lines && all(all(diff(dd.time_rng,1,2)))==0
%       time_axis = dd.time_rng(1,1):dd.dt:dd.time_rng(2,1);
%     else
%       fprintf('Time rannge not equal in all range lines');
%       continue;
%     end

%     [time_zero, time_zero_idx] = min(abs(time_axis));

    time_axis = dd.param_analysis.radar.wfs(wf).time;
    time_zero_idx = find(time_axis==0);
    time_axis_length = time_axis(end)-time_axis(1);
    
    
    TTL_start =   param.radar.wfs(wf).TTL_start * 1/param.radar.TTL_clock ;
    TTL_length =   param.radar.wfs(wf).TTL_length * 1/param.radar.TTL_clock ;
    TTL_end =  (param.radar.wfs(wf).TTL_start + param.radar.wfs(wf).TTL_length) * 1/param.radar.TTL_clock ;

    rec_start = param.radar.wfs(1).record_start*1/param.radar.fs ;
    rec_stop = param.radar.wfs(1).record_stop*1/param.radar.fs ;
    rec_length = rec_stop-rec_start;
    
    % Collate
    mean_wf_data = mean(abs(dd.wf_data),2);
    gain_dB = lp(mean_wf_data,2);
    
    low_Tlim = TTL_end - rec_start - dd.param_analysis.radar.wfs(wf).Tadc_adjust...
      -dd.param_analysis.radar.wfs(wf).time_correction;
    upper_Tlim = rec_length*0.8;
    low_idxs = find(dd.param_analysis.radar.wfs(wf).time>low_Tlim);
    upper_idxs = find(dd.param_analysis.radar.wfs(wf).time>upper_Tlim);
    
    Gain_trunc = mean_wf_data(low_idxs(1):upper_idxs(1));
    Gain = ( Gain_trunc ./ max(Gain_trunc) );
    Time = dd.param_analysis.radar.wfs(wf).time(low_idxs(1):upper_idxs(1));
    % Gain, Gain_trunc and Time have same size
    % Gain_raw and time_axis ( or dd.param_analysis.radar.wfs(wf).time) have same size
    
    % Below equation works. 
    % In data_load, 1./Inf = zero before low_Tlim
%     Gain_raw = [Gain(1)*Inf*ones(low_idxs(1)-1,1); Gain; Gain(end)*ones(upper_idxs(end)-upper_idxs(1),1) ]; 

    % Below equation doesnot work. 
    % In data_load, 1./Gain = amplified noise before low_Tlim
%     Gain_raw = [mean_wf_data(1:low_idxs(1)-1)./max(Gain_trunc); Gain; Gain(end)*ones(upper_idxs(end)-upper_idxs(1),1) ];

    Gain_raw = [ones(low_idxs(1)-1,1); Gain; Gain(end)*ones(upper_idxs(end)-upper_idxs(1),1) ];
    
    
    % dB
    gain_max_dB = lp(max(Gain_trunc),2);
    norm_gain_dB = gain_dB-gain_max_dB;
    
    % To save and use
    param_analysis = dd.param_analysis;
    param_records = dd.param_records;
    param_ftg = param;
    %% Generate the plots
    % =====================================================================
    cur_fig = h_fig(fig_idx);
    clf(cur_fig);
    if param.analysis.enable_visible_plot
      figure(cur_fig); % Brings the current figure to the top
    end
    h_axes = axes('parent',cur_fig);
    plot(h_axes, time_axis/1e-6,mean_wf_data, 'LineWidth', 1)
    hold(h_axes, 'on');
    plot(h_axes, Time/1e-6, Gain_trunc,'--', 'LineWidth', 3);
    xlabel(h_axes,'Fast Time (relative), us');
    ylabel(h_axes,'Gain');
    colors = get(h_axes,'ColorOrder');
    grid(h_axes, 'on');
    axis(h_axes, 'tight');
    line(h_axes, [low_Tlim low_Tlim]/1e-6, h_axes.YLim, 'Color','g', 'LineWidth',1);
    line(h_axes, [upper_Tlim upper_Tlim]/1e-6, h_axes.YLim, 'Color','g', 'LineWidth',1);
    legend(h_axes,sprintf('[%d-%d] G_m_a_x=%.2f',wf,adc,gain_max_dB), 'Location', 'South');
    hold(h_axes, 'off');
    if ~param_analysis.radar.wfs(wf).gain_en
      title(h_axes,sprintf('For [ %d-%d ] Fast Time Gain Curve ORIGINAL',wf,adc));
    else
      title(h_axes,sprintf('For [ %d-%d ] Fast Time Gain Curve CORRECTED',wf,adc));
    end

    %% Save the plots
    if ~dd.param_analysis.radar.wfs(wf).gain_en
      saveas(cur_fig, fullfile(out_dir,sprintf('gain_wf_%d_adc_%d.fig',wf,adc)));
      saveas(cur_fig, fullfile(out_dir,sprintf('gain_wf_%d_adc_%d.jpg',wf,adc)));
    else
      saveas(cur_fig, fullfile(out_dir,sprintf('gain_wf_%d_adc_%d_2.fig',wf,adc)));
      saveas(cur_fig, fullfile(out_dir,sprintf('gain_wf_%d_adc_%d_2.jpg',wf,adc)));
    end

    if ~param_analysis.radar.wfs(wf).gain_en
      save( fullfile(out_dir,sprintf('gain_wf_%d_adc_%d',wf,adc)), '-v7.3', ...
        'Gain', 'Time', 'low_Tlim', 'upper_Tlim', 'low_idxs', 'upper_idxs',...
        'param_ftg', 'param_analysis', 'param_records', 'Gain_raw', 'time_axis');
    else
      fprintf('Data was compensated for ftg. Therefore, new ftg datafile is not saved for wf-adc %d-%d\n',wf,adc);
    end
    %% Generate the RAW plots
    % =====================================================================
    if param.analysis.raw_plot_en
      fig2_idx = fig2_idx +1;
      cur_fig = h_fig2(fig2_idx);
      clf(cur_fig);
      if param.analysis.enable_visible_plot
        figure(cur_fig); % Brings the current figure to the top
      end
      h_axes = axes('parent',cur_fig);
      plot(h_axes, time_axis/1e-6,gain_dB,'LineWidth',2);
      hold(h_axes, 'on');
      xlabel(h_axes,'Fast Time (relative), us');
      ylabel(h_axes,'Gain, dB');
      colors = get(h_axes,'ColorOrder');
      grid(h_axes, 'on');
      axis(h_axes, 'tight');
      hold(h_axes, 'off');
      if ~param_analysis.radar.wfs(wf).gain_en
        title(h_axes,sprintf('For [ %d-%d ] Fast Time Gain Curve ORIGINAL',wf,adc));
      else
        title(h_axes,sprintf('For [ %d-%d ] Fast Time Gain Curve CORRECTED',wf,adc));
      end
    %% Save the RAW plots
      saveas(cur_fig, fullfile(out_dir,sprintf('ftg2_wf_%d_adc_%d.fig',wf,adc)));
      saveas(cur_fig, fullfile(out_dir,sprintf('ftg2_wf_%d_adc_%d.jpg',wf,adc)));
    end
    
    %% Generate the COMBINED plot
    % =====================================================================
    cur_fig = h_fig(comb_plot);
    if param.analysis.enable_visible_plot
      figure(cur_fig); % Brings the current figure to the top
    end
    h_axes = gca(cur_fig);
    plot(h_axes, time_axis/1e-6,norm_gain_dB,'LineWidth',2);
    hold(h_axes, 'on');
    leg{fig_idx} = sprintf('[%d-%d] G_m_a_x=%.2f',wf,adc,gain_max_dB) ; 

    %% Generate the COMBINED RAW plot
    % =====================================================================
    if param.analysis.raw_plot_en
      cur_fig = h_fig2(comb_plot2);
      if param.analysis.enable_visible_plot
        figure(cur_fig); % Brings the current figure to the top
      end
      h_axes = gca(cur_fig);
      plot(h_axes, time_axis/1e-6,gain_dB,'LineWidth',2);
      hold(h_axes, 'on');
      leg2{fig2_idx} = sprintf('[%d-%d] G_m_a_x=%.2f',wf,adc,gain_max) ; 
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

%% Wrap up the COMBINED plot
% =========================================================================
cur_fig = h_fig(comb_plot);
if param.analysis.enable_visible_plot
  figure(cur_fig); % Brings the current figure to the top
end
h_axes = gca(cur_fig);
colors = get(h_axes,'ColorOrder');
xlabel(h_axes,'Fast Time (relative), us');
ylabel(h_axes,'Normalized Gain');
grid(h_axes, 'on');
axis(h_axes, 'tight');
hold(h_axes, 'off');
if ~param_analysis.radar.wfs(wf).gain_en
  title(h_axes,sprintf('Combined Fast Time Gain Curve ORIGINAL'));
else
  title(h_axes,sprintf('Combined Fast Time Gain Curve CORRECTED'));
end
legend(h_axes,leg, 'Location', 'South');

%% Save the combined plot plot
if ~dd.param_analysis.radar.wfs(wf).gain_en
  saveas(cur_fig, fullfile(out_dir,sprintf('gain_combined_ch_%d.fig',adc)));
  saveas(cur_fig, fullfile(out_dir,sprintf('gain_combined_ch_%d.jpg',adc)));
else
  saveas(cur_fig, fullfile(out_dir,sprintf('gain_combined_ch_%d_2.fig',adc)));
  saveas(cur_fig, fullfile(out_dir,sprintf('gain_combined_ch_%d_2.jpg',adc)));
end

%% Wrap up the COMBINED RAW plot
% =========================================================================
if param.analysis.raw_plot_en
  cur_fig = h_fig2(comb_plot2);
  if param.analysis.enable_visible_plot
    figure(cur_fig); % Brings the current figure to the top
  end
  h_axes = gca(cur_fig);
  colors = get(h_axes,'ColorOrder');
  xlabel(h_axes,'Fast Time (relative), us');
  ylabel(h_axes,'Gain, dB');
  grid(h_axes, 'on');
  axis(h_axes, 'tight');
  hold(h_axes, 'off');
  if ~param_analysis.radar.wfs(wf).gain_en
    title(h_axes,sprintf('Combined Fast Time Gain Curve ORIGINAL'));
  else
    title(h_axes,sprintf('Combined Fast Time Gain Curve CORRECTED'));
  end
  legend(h_axes,leg2, 'Location', 'South');

  %% Save the combined RAW plot plot
  saveas(cur_fig, fullfile(out_dir,sprintf('ftg2_combined_ch_%d.fig',adc)));
  saveas(cur_fig, fullfile(out_dir,sprintf('ftg2_combined_ch_%d.jpg',adc)));
end

fprintf('=============================================================\n');
fprintf('%s: %s (%s)\n', param.season_name, param.day_seg, datestr(now));
fprintf('Output directory -->\n %s\n',out_dir);
fprintf('=============================================================\n');

% 


%%

%Surprise!
% plot_ftg;
if param.analysis.ftg_plot_en
% script plot_ftg
%
% Plots results from run_collate_ftg
%
% Authors: John Paden, Hara Madhav Talasila

%% USER SETTINGS
% =========================================================================

% param_override = [];
% param_sheet_name = 'rds_param_2019_Greenland_P3.xls';
% param_fn = ct_filename_param(param_sheet_name);
% params = read_param_xls(param_fn,'',{'analysis_gain' 'analysis'});

% param_override.analysis.enable_visible_plot = 0; % Set 1/0; Default: 1

%% Automated Section
% =========================================================================
% 
% % Input checking
% global gRadar;
% if exist('param_override','var')
%   param_override = merge_structs(gRadar,param_override);
% else
%   param_override = gRadar;
% end

% Process each of the segments
for param_idx =param_idx
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  %% General Setup
  % =======================================================================

%   param = merge_structs(param, param_override);
  physical_constants;

  fprintf('=============================================================\n');
  fprintf('%s: %s (%s)\n', param.season_name, param.day_seg, datestr(now));
  fprintf('=============================================================\n');
  %% Input Checks
  % =======================================================================

  if ~isfield(param.analysis, 'enable_visible_plot')
    param.analysis.enable_visible_plot = 1;
  end

  ftg_dir = ct_filename_ct_tmp(param,'','waveform','gain');
  if ~exist(ftg_dir, 'dir')
    fprintf('Empty directory --> run_collate_ftg (%s)\n',ftg_dir);
    return;
  end
  %% PLOT
  % For figure handles
  leg=[];
  leg_idx=0;
  color_idx=0;
  h_fig = get_figures(1,param.analysis.enable_visible_plot);
  cur_fig = h_fig(1);
  clf(cur_fig);
  if param.analysis.enable_visible_plot
    figure(cur_fig); % Brings the current figure to the top
  end
  h_axes = axes('parent',cur_fig);
  hold(h_axes, 'on');
  colors = get(h_axes,'ColorOrder');
  dd={};
  % The loop for each img-wf_adc pair
  for img = 1:length(param.analysis.imgs)
    for wf_adc = 1:size(param.analysis.imgs{img},1)
      wf = param.analysis.imgs{img}(wf_adc,1);
      adc = param.analysis.imgs{img}(wf_adc,2);
      %% Load the ftg file
      % =====================================================================
      dd{wf,adc} = load(fullfile(ftg_dir, sprintf('gain_wf_%d_adc_%d.mat',wf,adc)));
      color_idx = color_idx+1;
      plot(h_axes, dd{wf,adc}.Time/1e-6,dd{wf,adc}.Gain,'LineWidth', 2, 'Color', colors(color_idx,:));
      leg_idx = leg_idx+1;
      leg{leg_idx} = sprintf('[%d-%d] Truncated',wf,adc) ;
      plot(h_axes, dd{wf,adc}.param_analysis.radar.wfs(wf).time/1e-6, dd{wf,adc}.Gain_raw,':', 'LineWidth', 2, 'Color', colors(color_idx,:));
      leg_idx = leg_idx+1;
      leg{leg_idx} = sprintf('[%d-%d] Extrapolated',wf,adc) ;
      line(h_axes, [dd{wf,adc}.low_Tlim dd{wf,adc}.low_Tlim]/1e-6, [dd{wf,adc}.Gain(1), 0.00],'Color', 'k', 'LineWidth',1);
      leg_idx = leg_idx+1;
      leg{leg_idx} = sprintf('[%d-%d] Limits',wf,adc) ;
      line(h_axes, [dd{wf,adc}.upper_Tlim dd{wf,adc}.upper_Tlim]/1e-6, [dd{wf,adc}.Gain(end), 0.94],'Color', 'g', 'LineWidth',2);
      leg_idx = leg_idx+1;
      leg{leg_idx} = sprintf('[%d-%d] Limits',wf,adc) ;
      
    end
  end
  xlabel(h_axes,'Fast Time (relative), us');
  ylabel(h_axes,'Normalized Gain');
  colors = get(h_axes,'ColorOrder');
  grid(h_axes, 'on');
  axis(h_axes, 'tight');


  legend(h_axes,leg, 'Location', 'South');
  hold(h_axes, 'off');
  title(h_axes,sprintf('For chan [ %d ] Fast Time Gain Curve',adc));

  %% Save the plots
  saveas(cur_fig, fullfile(ftg_dir,sprintf('gain_norm_chan_%d.fig',adc)));
  saveas(cur_fig, fullfile(ftg_dir,sprintf('gain_norm_chan_%d.jpg',adc)));
end

end % if param.analysis.ftg_plot_en