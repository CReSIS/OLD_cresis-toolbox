function pf_plot(param)
% This function is called from the main particle filtr scripts.
% It generates few plots. You can add more if you need to.
%
% Author: Mohanad Al-Ibadi
%
%%                               Define parameters
% =========================================================================
slice          = param.slice;          % Required slices to plot
est_doa        = param.est_doa;        % Estimated DoA
bin_rng        = param.bin_rng;        % Range-bins
Nsig           = param.Nsig;           % Maximum number of allowed targets
doa_lim_st     = param.plot_lim_st;     % Start of the doa (usually this is the start of the 3dB beamwidth)
doa_lim_end    = param.plot_lim_end;    % End of the doa (usually this is the end of the 3dB beamwidth)
est_method     = param.est_method;     % Type of DoA estimation: MAP, MMSE, or averaged MAP
pf_method      = param.pf_method;      % Type of particle filter: 'standard', RPF, or MCMC
rbins          = param.rbins;          % Actual range-bins
rlines         = param.rlines;         % Actual range-lines
Ntrials        = param.Ntrials;        % Number of Monte Carlo runs used to generate the RMSE;

%%         Determine the first and last active range-bins of the surface
% =========================================================================
if isfield(param,'actual_doa') && ~isempty(param.actual_doa)
  % Actual DoA, if exists
  actual_doa = param.actual_doa;
  idx_1_max = find(~isnan(squeeze(actual_doa(:,slice,1))),1,'last');
  idx_2_max = find(~isnan(squeeze(actual_doa(:,slice,2))),1,'last');
  
  idx_1_min = find(~isnan(squeeze(actual_doa(:,slice,1))),1,'first');
  idx_2_min = find(~isnan(squeeze(actual_doa(:,slice,2))),1,'first');
else
  % Estimated DoA
  idx_1_max = find(~isnan(squeeze(est_doa(:,slice,1))),1,'last');
  idx_2_max = find(~isnan(squeeze(est_doa(:,slice,2))),1,'last');
  
   idx_1_min = find(~isnan(squeeze(est_doa(:,slice,1))),1,'first');
   idx_2_min = find(~isnan(squeeze(est_doa(:,slice,2))),1,'first');
end

if ~isempty(idx_1_max) && ~isempty(idx_2_max)
  max_y_lim = max(idx_1_max,idx_2_max);
elseif isempty(idx_1_max) || isempty(idx_2_max)
  if isempty(idx_1_max)
    max_y_lim = idx_2_max;
  else
    max_y_lim = idx_1_max;
  end
else
  max_y_lim = param.Nt_end-max(bin_rng);
end
% max_y_lim = max_y_lim + param.Nt_st;

if ~isempty(idx_1_min) && ~isempty(idx_2_min)
  min_y_lim = min(idx_1_min,idx_2_min);
elseif isempty(idx_1_min) || isempty(idx_2_min)
  if isempty(idx_1_min)
    min_y_lim = idx_2_min;
  else
    min_y_lim = idx_1_min;
  end
else
  min_y_lim = param.Nt_st+max(bin_rng);
end
% min_y_lim = min_y_lim + param.Nt_st;

% if ~isempty(idx_1_max) || ~isempty(idx_2_max) 
%   if isempty(idx_1_max)
%     max_y_lim = idx_2_max+max(bin_rng);
%   else
%     max_y_lim = idx_1_max+max(bin_rng);
%   end
% else
%   max_y_lim = param.Nt_end-max(bin_rng);
% end
% 
% if ~isempty(idx_1_min) || ~isempty(idx_2_min)
%   if isempty(idx_1_min)
%     min_y_lim = idx_2_min+max(bin_rng);
%   else
%     min_y_lim = idx_1_min+max(bin_rng);
%   end
% else
%   min_y_lim = param.Nt_st+max(bin_rng);
% end

%%                               Plot
% =========================================================================
%% DoA vs Range
% -------------
if 0
figure(999);clf;
hold on
for signal_idx=1:Nsig
  if isfield(param,'actual_doa') && ~isempty(param.actual_doa)
    % Actual DoA, if exists
    actual_doa = param.actual_doa;
%     actual_doa_plot = [NaN(param.Nt_st,1);squeeze(actual_doa(:,slice,signal_idx));NaN(max(bin_rng),1)];
    actual_doa_plot = squeeze(actual_doa(:,slice,signal_idx));
    scatter(actual_doa_plot,1:length(actual_doa_plot), 80,'+','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
  end
  
  % Estimated DoA, always exists
%   doa_tmp_plot = [NaN(param.Nt_st,1);squeeze(est_doa(:,slice,signal_idx));NaN(max(bin_rng),1)];
  doa_tmp_plot = squeeze(est_doa(:,slice,signal_idx));
  scatter(doa_tmp_plot,1:length(doa_tmp_plot), 20,'s','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',2);
  hold on
  
  if isfield(param,'smoothing_method') && ~isempty(param.smoothing_method)
    % Smoothed DoA, if exists
    smoothed_est_doa = param.smoothed_est_doa;
%     smoothed_tmp_plot = [NaN(param.Nt_st,1);squeeze(smoothed_est_doa(:,slice,signal_idx));NaN(max(bin_rng),1)];
    smoothed_tmp_plot = squeeze(smoothed_est_doa(:,slice,signal_idx));
    scatter(smoothed_tmp_plot,1:length(smoothed_tmp_plot), 20,'^','MarkerFaceColor','m','MarkerEdgeColor','m','LineWidth',2);
  end
end

lim_guard = round((max_y_lim-min_y_lim)/10);
lim_guard = max(lim_guard,5);
  
xlim([doa_lim_st doa_lim_end]*180/pi + [-5 5])
ylim([min_y_lim max_y_lim] + [-lim_guard lim_guard])

set(gca,'YDir','reverse');
xlabel('DoA (deg.)')
ylabel('Range bin')
title(sprintf('%s (%s particle filter): slice # %u',upper(est_method),pf_method,slice))
grid on

if  isfield(param,'actual_doa') && ~isempty(param.actual_doa) && ...
    isfield(param,'smoothing_method') && ~isempty(param.smoothing_method)
   legend('Actual DOA',sprintf('PF:%s-before %s smoothing',param.pf_method,param.smoothing_method),sprintf('PF:%s-after %s smoothing',param.pf_method,param.smoothing_method),'Location','best');    
%   legend('Est. DoA-before smoothing',sprintf('Est. DoA-after %s smoothing',param.smoothing_method), 'Actual DoA','Location','best');
elseif isfield(param,'actual_doa') && ~isempty(param.actual_doa) && ...
    (~isfield(param,'smoothing_method') || isempty(param.smoothing_method))
%   legend('Est. DoA', 'Actual DoA','Location','best')
    legend('Actual DOA',sprintf('PF: %s',param.pf_method) ,'Location','best')
elseif (~isfield(param,'actual_doa') || isempty(param.actual_doa)) && ...
    isfield(param,'smoothing_method') && ~isempty(param.smoothing_method)
   legend(sprintf('PF:%s-before %s smoothing',param.pf_method,param.smoothing_method),sprintf('PF:%s-after %s smoothing',param.pf_method,param.smoothing_method),'Location','best');    
%   legend(sprintf('Est. DoA-before %s smoothing',param.smoothing_method),sprintf('Est. DoA-after %s smoothing',param.smoothing_method),'Location','best');
else
  legend('Est. DoA','Location','best');
end
end

%% True vs estimated DoA
% ----------------------
if 0
if isfield(param,'actual_doa') && ~isempty(param.actual_doa)
  actual_doa = param.actual_doa;
  figure(9990);clf
  for signal_idx=1:Nsig
    scatter(actual_doa(:,slice,signal_idx),est_doa(:,slice,signal_idx),80,'x','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',2);
    hold on
    
    if isfield(param,'smoothing_method') && ~isempty(param.smoothing_method)
      smoothed_est_doa = param.smoothed_est_doa;
      scatter(actual_doa(:,slice,signal_idx),smoothed_est_doa(:,slice,signal_idx),30,'s','fill','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',2);
    end
  end
  
  xlim([doa_lim_st doa_lim_end]*180/pi + [-5 5])
  ylim([doa_lim_st doa_lim_end]*180/pi + [-5 5])
  
  xlabel('True DoA (deg.)')
  ylabel('Estimated DoA (deg.)')
  title(sprintf('%s (%s particle filter): slice # %u',upper(est_method),pf_method,slice))
  grid on
  
  if isfield(param,'smoothing_method') && ~isempty(param.smoothing_method)
    legend(sprintf('Before %s smoothing',param.smoothing_method),sprintf('After %s smoothing',param.smoothing_method),'Location','best')
  end
end
end
%% RMSE before smoothing is applied
% ---------------------------------
if 0
if isfield(param,'actual_doa') && ~isempty(param.actual_doa)
  rmse           = param.rmse;           % RMSE before smoothing
  avg_rbin_rmse  = param.avg_rbin_rmse;  % RMSE per range-bin before smoothing
  avg_rline_rmse = param.avg_rline_rmse; % RMSE per range-line before smoothing
  
  figure(9991);clf;
  subplot(311)
  imagesc(rmse)
  h_colorbar = colorbar;
  set(get(h_colorbar,'YLabel'),'String','RMSE (deg.)');
  xlabel('Range line');
  ylabel('Range bin');
  title(sprintf('%s: Average RMSE over %u runs',upper(est_method),Ntrials))
  ylim([min_y_lim max_y_lim] + [-5 5])
  caxis([min(rmse(:)) max(rmse(:))])
  
  subplot(312)
  scatter(1:length(avg_rbin_rmse),avg_rbin_rmse,20,'fill', 'MarkerFaceColor','b');
  xlabel('Range-bin')
  ylabel('RMSE (deg.)')
  title(sprintf('%s: Average RMSE (over range-lines)',upper(est_method)))
  grid on
  
  subplot(313)
  scatter(1:length(avg_rline_rmse),avg_rline_rmse,20,'fill', 'MarkerFaceColor','b');
  xlabel('Range-line')
  ylabel('RMSE (deg.)')
  title(sprintf('%s: Average RMSE (over range-bins)',upper(est_method)))
  grid on
end
end
%% RMSE after smoothing is applied
% ---------------------------------
if 0
if isfield(param,'smoothing_method') && ~isempty(param.smoothing_method)
  rmse_smoothed_doa = param.rmse_smoothed_doa;
  avg_rline_rmse_smoothed_doa = param.avg_rline_rmse_smoothed_doa;
  avg_rbin_rmse_smoothed_doa = param.avg_rbin_rmse_smoothed_doa;
  
  figure(9992);clf;
  subplot(311)
  imagesc(rmse_smoothed_doa)
  h_colorbar = colorbar;
  set(get(h_colorbar,'YLabel'),'String','RMSE (deg.)');
  xlabel('Range line');
  ylabel('Range bin');
  title(sprintf('%s: Average RMSE with %s smoothing over %u runs',upper(est_method),param.smoothing_method,Ntrials))
  ylim([min_y_lim-5 max_y_lim+5])
  
  subplot(312)
  scatter(1:length(avg_rbin_rmse_smoothed_doa),avg_rbin_rmse_smoothed_doa,20,'fill', 'MarkerFaceColor','b');
  xlabel('Range-bin')
  ylabel('RMSE (deg)')
  title(sprintf('%s: Average RMSE (over range-lines) with %s smoothing',upper(est_method),param.smoothing_method))
  grid on
  
  subplot(313)
  scatter(1:length(avg_rline_rmse_smoothed_doa),avg_rline_rmse_smoothed_doa,20,'fill', 'MarkerFaceColor','b');
  xlabel('Range-line')
  ylabel('RMSE (deg)')
  title(sprintf('%s: Average RMSE (over range-bins) with %s smoothing',upper(est_method),param.smoothing_method))
  grid on
end
end
%% RMSE vs TEST_PARAM (Nsnaps, SNR, etc) averaged over all range-lines and range-bins
% -----------------------------------------------------------------------------------
if isfield(param,'test_param') && ~isempty(param.test_param)
  rmse_tests_mean = param.rmse_tests_mean;
  test_param = param.test_param;
  figure(9993);hold on
  if strcmp(param.test_type,'Number of snapshots')
    semilogx(test_param,rmse_tests_mean,'*','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',2);
  else
    plot(test_param,rmse_tests_mean,'*','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',2);
  end
  if test_param(end) > test_param(1)
  xlim([test_param(1)  test_param(end)])
  end
  if isfield(param,'test_type') && ~isempty(param.test_type)
    xlabel(param.test_type)
  end
  ylabel('RMSE (deg.)')
  title(sprintf('%s: %s vs average RMSE (over range-lines/bins)',upper(est_method),param.test_type))
  grid on
  legend('1000 particles: standard PF','1000 particles: MCMC PF','5000 particles: standard PF','5000 particles: MCMC PF')
end

return