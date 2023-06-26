% script sim.crosstrack_example_pf
%
% Example particle filter simulation. Calls crosstrack.m and
% crosstrack_rmse.m.
%
% Author: Mohanad Al-Ibadi, Sravya Athinarapu, Sean Holloway, John Paden,
% Theresa Stumpf
%
% See also: sim.crosstrack.m, sim.crosstrack_data.m,
% sim.crosstrack_example_*.m, sim.crosstrack_rmse.m

%% Particle filter parameters

physical_constants;

% Center frequency
param.fc                 = fc;
% Signal bandwidth
param.BW                 = BW;
% Sampling frequency
param.fs                 = BW;
% Maximum number of targets. Currently it is always 2
param.Nsrc               = array_param.Nsrc;
% Number of sensors
param.Nc                 = Nc;
% Maximum number of particles to start with. Set=Np_fixed for now untill it
% is fully supported in the right way.
param.Np_tmp             = 1000;
% After PF has converged, this will be the new Np
param.Np_fixed           = 1000;
% Switch from Np_tmp to Np_fixed after this many range-bins from the starting bin
bins_to_switch           = 0;
% Type of smoothing: 'Viterbi' or 'backward'. Backward smoothing requires
% larger apriori variance to give good results, but the large variance may
% ruin the PF results.
param.smoothing_method   = '';
% Type of particle filter: 'standard', RPF, or MCMC
param.pf_method          = 'standard';
% Location of the regularization step: 'before' or 'after' resampling
param.reg_loc            = 'before';
% Kernel used in RPF and MCMC PF: 'Epanechnikov', 'Gaussian without whitening',
% or 'Gaussian with whitening'
param.kernel_type        = 'Gaussian without whitening';
% Type of DoA estimation: map, mmse, quasi map, or hybrid
% TO DO: MAKE SURE THAT quasi MAP DOESN'T INCLUDE PARTICLES FROM THE WRONG
% MODE.
param.est_method         = 'hybrid';
% Choose resampling method: systematic, stratified, standard, or randsample
param.resamp_method      = 'systematic';
% Model of the change in DoA per range-bin: exact (sum previous changes)
% or approx (n*average DoA change over all bins)
param.doa_change_model   = 'exact';
% Initial prior distribution: Gaussian or uniform.
param.init_proposal_dist = 'uniform';
% Slice index that we need to plot
param.slice              = 1;
% Resample if the number of effectivs particles is < resampling_thre*Np
param.resampling_thr     = 1/2;
% In dB. Used only when the actual number of targets is unknown. If SNR of
% the DoA drops bellow this threshold, then this DoA is treated as noise
% or false target. This threshold should be around 3dB bellow the average
% SNR (estimated from data). If set to inf, then there is no thresholding.
param.snr_thr            = inf;
% Value >0 and <1. Used only when the actual number of targets is unknown.
% In case of more than one DoAs, any DoA has a likelihood less than this
% threshold will be treated as noise or false target.
param.likelihood_thr     = 0.0;

% Estimated (roughly) SNR in dB over all sensors. Required if CRLB is used
% to estimate the variance of the process model. CHANGE IF DATA HAS CHANGED.
% param.snr = 0.4+10*log10(param.Nc);

% phase_center = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
%  param.y_pc = phase_center(2,:).';
%  param.z_pc = phase_center(3,:).';

param.k = k;
param.steering_mtx_fh = @(DOA) ((1/sqrt(length(param.y_pc)))*exp(1i*(param.y_pc*k*sin(DOA) - param.z_pc*k*cos(DOA))));

% param.tx_weights    = [hanning(4).',0 0 0 0].';
param.tx_weights    = hanning(param.Nc);
param.steering_ang  = 0*pi/180;
steering_delay      = (param.y_pc*sin(param.steering_ang) - param.z_pc*cos(param.steering_ang))/c;
param.weight        = param.tx_weights.*exp(-1i*4*pi*param.fc*steering_delay); % Sensors complex weights
% Angular range for radiation patetrn, RP, calulation
theta_RP = linspace(-90,90,2048)*pi/180;
RP       = abs(param.weight.'*param.steering_mtx_fh(theta_RP)).^2;
RP       = RP./max(RP);
RP_dB    = 10*log10(RP);
% 3dB threshold. Here we take a little more than 3dB.
threshold_dB = -3.5; % in dB
idxs_3dB    = find(RP_dB>=threshold_dB);
RP_3dB      = RP_dB(idxs_3dB);
theta_3dB   = theta_RP(idxs_3dB);
% Start/end of the DoA limits
param.doa_lim_st  = min(ang) - 5*pi/180; %min(theta_3dB) ;
param.doa_lim_end =  max(ang) + 5*pi/180; %max(theta_3dB) ;

% Start/end of the DoA limits for plotting
param.plot_lim_st  = min(ang); %min(theta_3dB);
param.plot_lim_end = max(ang); %max(theta_3dB);

comp_tx_weight = param.weight .'*param.steering_mtx_fh(ang.');
param.src.tx_weights = comp_tx_weight;

% Number of fast-time snapshots
Nt_snaps = length(array_param.bin_rng);
% Number of azimuth snapshots
Nx_snaps = length(array_param.line_rng);
% Total number of snapshots
param.M  = Nt_snaps*Nx_snaps;

% Type of likelihood density: either taking the expectation/mean of
% the likelihood densities of all snapshots, or taking the
% product of the likelihood densities of all snapshots.
% Function handel to the (complex) likelihood density function
param.likelihood_dens_fh = @(x,sigma) (((pi*sigma)^(-param.Nc))*(1/param.M)*sum(exp(-((sum(abs(x).^2,1))./sigma))));% Mean
%   param.likelihood_dens_fh = @(x,sigma) (((pi*sigma)^(-param.Nc*param.M))*prod(exp(-((sum(abs(x).^2,1))./(sigma))))); Product

% Function handel to the (real) apriori density function
param.apriori_dens_fh    = @(x,mu,sigma,Ndoa) (((2*pi*sigma).^(-Ndoa/2)).*exp(-(sum(abs(x-mu).^2,2))./(2*sigma)));

% const is a scaling factor for the measurements noise variance.
% Setting this number too large (e.g. 20) may turn the Gaussian to uniform,
% unless you increase Np. Setting it to a too small number (e.g. 2) may make
% the Gaussian too narrowand and will not cover the support of posteropr, so some
% DoAs will not show up. A good value for const is usually 1 (high SNR),
% 2 (medium SNR), or 10 (low SNR).
param.const = 1;

%% Prepare the approximate change in DoA from range-bin to the next
% -----------------------------------------------------------------
% Change of DoA from rbin to rbin (max at nadir). Depends on the range
% resolution and the flight hight. The assumption here is that the surface
% is flat.
if 1
  % Simulation data
  % dt: Fast time sample spacing (sec)
  dt = 1/param.fs;
  % Nt: number of range bins (fast time samples)
  param.Nt = floor((param.src.t1-param.src.t0)/dt);
  time = param.src.t0 + dt*(0:param.Nt-1).';
  R = time*c/2;
  param.R = R;
  rng_res = c/(2*param.BW);
  param.rng_res = rng_res;
  h = flight_height;
  % Ignore range-bins above the surface
  good_R_idx = find(R>=h);
  good_R = R(good_R_idx);
  % Now calculate the delta theta
  delta_theta_tmp = diff(acos(h./good_R));
  delta_theta = [zeros(good_R_idx(1),1) ; delta_theta_tmp];
elseif 0
  % Real data: CHECK
  rng_res = c/(2*param.BW);
  param.pulse_width = 3e-6;
  h = param.pulse_width*c;
  R = data.fcs.surface(param.Nt_st:param.Nt_end)*c;
  good_R_idx = find(R>=h);
  theta_1 = acos(h./R(good_R_idx(1:end-1)));
  theta_2 = acos(h./(R(good_R_idx(2:end))+rng_res));
  delta_theta = zeros(param.Nt,1);
  delta_theta(good_R_idx(1)+1:end) = theta_2-theta_1;
  
  param.R = R;
  param.rng_res = rng_res;
end

if 0
  % Debug: Plot the change in DOA as a function of DOA
  figure(1000);clf
  plot_doa = acos(h./good_R(1:end-1))*180/pi;
  plot(plot_doa,delta_theta_tmp*180/pi,'b')
  xlabel('\theta^\circ')
  ylabel('\Delta\theta^\circ')
  title('Change in DOA as a function of DOA for a flat surface')
  grid on
end

% delta_theta = [zeros(good_R_idx(1)-1,1) ;(theta_2-theta_1)]; % (Nx-1)*1
if strcmp(param.doa_change_model,'exact')
  % Flat surface approximation..but uses the exact model (sum of previous
  % step sizes or DoA changes per range-bin)
  param.delta_theta = delta_theta;
  
elseif strcmp(param.doa_change_model,'approx')
  % Assumes the DoA change is constant from range-bin to the next. It
  % approximates the model as number of dead range-bins times the constant
  % step size.
  param.doa_change_per_rbin = mean(delta_theta);
end

% First/last range-bin index with targets.
param.Nt_st = good_R_idx(1)-max(array_param.bin_rng);
param.Nt_end = param.Nt - max(array_param.bin_rng);

% First/last range-line index with targets.
param.Nx_st = 1 + max(array_param.line_rng);
param.Nx_end = length(surf_param.x) - max(array_param.line_rng);

param.rbins  = [param.Nt_st:array_param.dbin:param.Nt_end];
param.rlines = [param.Nx_st:array_param.dline:param.Nx_end];

param.Nx = length(param.rlines);
param.Nt = length(param.rbins);

% clear array_param;
param.bin_rng = array_param.bin_rng;
param.line_rng = array_param.line_rng;

%% Recursion
param.Ntrials = 5; % =1 for real data
tic
for trial_idx = 1:param.Ntrials
  % Reset the RNG (Don not put it outside this loop).
  %   rng default
  param.monte.rng_seed = trial_idx+test_idx;
  
  % Generate input data and determine actual number of sources
  results = crosstrack(param);
  %   results = sim.crosstrack(param);
  
  sim_data = squeeze(results.sim_data{1}); % Nt*Nx*Nc
  if Nx_snaps == 1
    sim_data_tmp(:,1,:) = sim_data;
    sim_data = sim_data_tmp;
  end
  actual_doa_cell = results.actual_doa;    % Cell array
  
  param.data_in    =  sim_data;
  
  if 1
    % Restrict the number of DOAs per bin to 2 DOAs only
    actual_doa = cell(param.Nt,param.Nx);
    for lineIdx = 1:param.Nx
      for binIdx = 1:param.Nt
        if ~isempty(actual_doa_cell{param.rbins(binIdx),lineIdx})
          doa_tmp = sort(actual_doa_cell{param.rbins(binIdx),lineIdx},'ascend');
          if length(doa_tmp)<array_param.Nsrc
            doa = NaN(array_param.Nsrc,1);
            doa(1:length(doa_tmp)) = doa_tmp;
          else
            doa = doa_tmp(1:array_param.Nsrc);
          end
          actual_doa{binIdx,lineIdx} = doa; % Nt*Nsrc*Nx
        end
      end
    end
    param.actual_doa = actual_doa;
  else
    param.actual_doa = actual_doa_cell;
  end
  
  est_doa          = NaN(size(actual_doa_cell,1),size(actual_doa_cell,2),param.Nsrc);
  smoothed_est_doa = NaN(size(actual_doa_cell,1),size(actual_doa_cell,2),param.Nsrc);
  actual_doa_tmp   = NaN(size(actual_doa_cell,1),size(actual_doa_cell,2),param.Nsrc);
  est_doa_pwr      = NaN(size(actual_doa_cell,1),size(actual_doa_cell,2),param.Nsrc);
  
  rline_idx = 0;
  for rline =  param.rlines
    rline_idx = rline_idx + 1;
    fprintf('\n\n Test %u of %u ... Trial %u of %u ... processing range-line %u of %u\n', ...
      test_idx, Ntests,trial_idx,param.Ntrials,rline_idx,length(param.rlines))
    
    param.rline_idx           = rline_idx;
    param.Np                  = param.Np_tmp;
    
    % Process noise standard deviation
    param.doa_std             = 10*pi/180;
    % Process noise variance
    param.sigma_proposal      = param.doa_std^2;
    % Inital particles weight
    param.w                   = ones(param.Np,1)./param.Np;
    
    % Inital likelihood density
    %     param.likelihood_dens     = ones(param.Np,1)./param.Np;
    
    if 0
      % Initial left DoA (or left mode center)
      st_doa_L = (param.doa_lim_st+param.steering_ang)/2;
      % Initial right DoA (or right mode center)
      st_doa_R = (param.doa_lim_end+param.steering_ang)/2;
      % Initial DoA
      param.st_doa = [st_doa_L;st_doa_R];
    elseif 1
      param.st_doa  = 0*ones(param.Nsrc,1);
    end
    param.init_ref_doa = mean(param.st_doa);
    
    if strcmp(param.init_proposal_dist,'Gaussian')
      % Initial left mode std
      init_std_L = abs(param.steering_ang-param.doa_lim_st)/2;
      % Initial right mode std
      init_std_R = abs(param.doa_lim_end-param.steering_ang)/2;
      % Initial std (applied to both surface modes)
      init_std   = 2*max(init_std_L,init_std_R);
      param.doa_init_std  = init_std;
      %       param.doa_init_std  = 60*pi/180;
      % Initial Process noise variance. Used in the initial state only
      param.init_sigma_proposal = param.doa_init_std^2;
    end
    
    % May be needed in the multimodel algorithm.
    param.init_mu_proposal    = param.st_doa;
    
    % Range-bin index to switch from Np_tmp to Np_fixd
    param.switch_rbin = param.Nt_st +max(param.bin_rng) + bins_to_switch;
    
    % Call particle filter function
    pf_results = particle_filter(param);
    
    % Collect results
    est_doa(param.Nt_st:param.Nt_end,rline_idx,:)          = pf_results.est_doa;          % In rad
    smoothed_est_doa(param.Nt_st:param.Nt_end,rline_idx,:) = pf_results.smoothed_est_doa; % In rad
    est_doa_pwr(param.Nt_st:param.Nt_end,rline_idx,:)      = pf_results.est_doa_pwr;      % Linear units
    
    if isfield(param,'actual_doa') && ~isempty(param.actual_doa)
      actual_doa_tmp(param.Nt_st:param.Nt_end,rline_idx,:)   = pf_results.actual_doa_tmp;
      
      bad_actual_doa_idx_L = find(actual_doa_tmp(:,rline_idx,1) == param.init_ref_doa);
      bad_actual_doa_idx_R = find(actual_doa_tmp(:,rline_idx,2) == param.init_ref_doa);
      
      if isempty(bad_actual_doa_idx_L)
        if isnan(est_doa(bad_actual_doa_idx_R,rline_idx,2)) & ~isnan(est_doa(bad_actual_doa_idx_R,rline_idx,1))
          actual_doa_tmp(bad_actual_doa_idx_R,rline_idx,1) = actual_doa_tmp(bad_actual_doa_idx_R,rline_idx,2);
          actual_doa_tmp(bad_actual_doa_idx_R,rline_idx,2) = NaN;
        end
      else
        if isnan(est_doa(bad_actual_doa_idx_L,rline_idx,1)) & ~isnan(est_doa(bad_actual_doa_idx_L,rline_idx,2))
          actual_doa_tmp(bad_actual_doa_idx_L,rline_idx,2) = actual_doa_tmp(bad_actual_doa_idx_L,rline_idx,1);
          actual_doa_tmp(bad_actual_doa_idx_L,rline_idx,1) = NaN;
        end
      end
    end
    
  end
  
  % Ignore DOAs outside the 3dB beampattern
  est_doa(est_doa<param.doa_lim_st | est_doa>param.doa_lim_end) = NaN;
  est_doa = est_doa*180/pi;
  
  %% DoA error calculation
  if isfield(param,'actual_doa') && ~isempty(param.actual_doa)
    actual_doa_tmp(actual_doa_tmp<param.doa_lim_st | actual_doa_tmp>param.doa_lim_end) = NaN;
    actual_doa_tmp = actual_doa_tmp*180/pi;
    
    Error   = squeeze(nanmean(abs(actual_doa_tmp - est_doa),3)); % Nt*Nx, mean error over range-bins
    param.est_doa_error{trial_idx} = Error;
    
    if isfield(param,'smoothing_method') || isempty(param.smoothing_method)
      
      smoothed_est_doa(smoothed_est_doa<param.doa_lim_st | smoothed_est_doa>param.doa_lim_end) = NaN;
      smoothed_est_doa   = smoothed_est_doa*180/pi;
      
      Error_smoothed_doa = squeeze(nanmean(abs(actual_doa_tmp - smoothed_est_doa),3)); % Nt*Nx
      param.smoothed_doa_error{trial_idx} = Error_smoothed_doa;
    end
  end
end

%% RMSE calculation
if isfield(param,'actual_doa') && ~isempty(param.actual_doa)
  % Call the RMSE function
  rmse_results = pf_rmse(param);
  
  param.rmse           = rmse_results.rmse;
  param.avg_rline_rmse = rmse_results.avg_rline_rmse;
  param.avg_rbin_rmse  = rmse_results.avg_rbin_rmse;
  
  if isfield(param,'smoothing_method') && ~isempty(param.smoothing_method)
    param.rmse_smoothed_doa           = rmse_results.rmse_smoothed_doa;
    param.avg_rline_rmse_smoothed_doa = rmse_results.avg_rline_rmse_smoothed_doa;
    param.avg_rbin_rmse_smoothed_doa  = rmse_results.avg_rbin_rmse_smoothed_doa;
  end
  
  rmse_tests.rmse(:,:,test_idx)         = param.rmse;
  rmse_tests.avg_rline_rmse(:,test_idx) = param.avg_rline_rmse;
  rmse_tests.avg_rbin_rmse(:,test_idx)  = param.avg_rbin_rmse;
  test_param(test_idx) = param.M;
  %   test_param(test_idx) = Bf(test_idx);
end

end

%% Plot
slice = param.slice; % Surface of this slice will be plotted
param.est_doa    = est_doa;
if isfield(param,'actual_doa') && ~isempty(param.actual_doa)
  param.actual_doa = actual_doa_tmp;
end

if isfield(param,'smoothing_method') || isempty(param.smoothing_method)
  param.smoothed_est_doa = smoothed_est_doa;
end

rmse_tests_rmse = reshape(rmse_tests.rmse,[size(rmse_tests.rmse,1)*size(rmse_tests.rmse,2)  size(rmse_tests.rmse,3)]); % (Nt*Nx) by Ntests
rmse_tests_mean = nanmean(rmse_tests_rmse,1); % Ntests by 1
param.rmse_tests_mean = rmse_tests_mean;
param.test_param = test_param;
% Call the plotting functions
pf_plot(param)

if 0
  %% To plot Slice model
  
  slice = 11;
  surface_z = results.z_grid-results.param.monte.target_param{1}.z.mean;
  surface_z_groung_truth =results.surf_model.z-results.param.monte.target_param{1}.z.mean;
  figure(5); clf;  plot(results.surf_model.y, surface_z_groung_truth(:,slice),'b');
  hold on
  plot(results.surf_model.y, surface_z(:,slice),'r');
  xlim([-1500 1500])
  ylim([-200 250])
  title('Slice - surface model');
  xlabel('Cross-track (m)');
  ylabel('WGS84-Elevation (m)');
  hold off
  legend('Ground-truth surface','Actual surface');
  grid on
  
  % Slice - Range Bin v/s DOA
  figure(6),clf
  scatter(results.tomo.doa(:,1,slice)*(180/pi),results.array_param.bins, 20 , 10*log10(results.tomo.power(:,1,slice)),'fill');
  colorbar
  hold on
  scatter(results.tomo.doa(:,2,slice)*(180/pi),results.array_param.bins, 20 , 10*log10(results.tomo.power(:,2,slice)),'fill');
  colorbar
  set(gca,'Ydir','reverse')
  xlim([-60 60])
  ylim([1 100])
  title('Slice');
  xlabel('DOA (deg)');
  ylabel('Range bin');
  cc = caxis;
  h_cb = colorbar;
  set(get(h_cb,'YLabel'),'String','Relative power (dB)');
  grid on
end
toc

if param.debug_level >= 3
  return
end

%% Process and save the outputs
out_fn_dir = 'D:\tmp\RADAR_CONF\';
out_fn_name = 'example1';

RMSE = sim.crosstrack_rmse(param,results);

figure(1); clf
plot(param.monte.SNR(:,1)+10*log10(3),RMSE(:,:,1).','.','LineWidth',2);
grid on
xlabel('Source SNR (dB)')
ylabel('RMS error ( \circ )')
legend({'MUSIC','MLE','WB','WBMLE'})
title('Figure 2 from Wax and Ziskind 1988')

figure(2); clf
plot(param.monte.SNR(:,2)+10*log10(3),RMSE(:,:,2).','.','LineWidth',2);
grid on
xlabel('Source SNR (dB)')
ylabel('RMS error ( \circ )')
legend({'MUSIC','MLE','WB','WBMLE'})
title('Complement to figure 2 from Wax and Ziskind (source 2)')

% Save Outputs
if ~exist(out_fn_dir,'dir')
  mkdir(out_fn_dir);
end
out_fn = fullfile(out_fn_dir,out_fn_name);
savefig(1,[out_fn '_src1.fig']);
savefig(2,[out_fn '_src2.fig']);
save([out_fn '.mat'],'param','results')
