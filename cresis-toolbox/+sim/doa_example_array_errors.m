% script sim.doa_example_array_errors
%
% Example setup scripts for the sim.doa function. To run, enable one
% of the if {0|1} blocks below and then run the script.
%
% Author: John Paden, Theresa Stumpf
%
% See also: sim.doa.m

%% Array errors (uncalibrated array) on estimated DoAs
% =======================================================================

%% Simulation parameters
physical_constants;
param = [];

% Source parameters
fc     = 195e9;
BW     = 10e6;
param.src.f0 = fc-BW/2;
param.src.f1 = fc+BW/2;

Nc = 7;
lambda = c/fc;
d_y = lambda/2;
param.src.phase_center =  zeros(3,Nc);
param.src.phase_center(2,:) = d_y/2*(0:Nc-1); % -d_y/2*(0:Nc-1)
param.src.ft_wind                 = @(N) hanning(N);
%   param.src.lever_arm.fh            = @sim.lever_arm_example;
%   param.src.lever_arm.args          = {[],1,[1:Nc],[0; c/param.fc/2; 0]};
% DOA method parameters
param.method.list                   = [7];
param.method.Nsv                    = 32;
param.method.OneD_Nsv               = 128;

theta = [-1.2154, -1.0654, -0.9484, -0.8481, -0.7580, -0.6751, -0.5974, -0.5236, -0.4528, -0.3844, -0.3178, -0.2527, ...
  -0.1886, -0.1253, -0.0625, 0, 0.0625, 0.1253, 0.1886, 0.2527, 0.3178, 0.3844, 0.4528, 0.5236, 0.5974, 0.6751, ...
  0.7580, 0.8481, 0.9484, 1.0654, 1.2154]';

theta_left = flipud(theta(1:floor(length(theta)/2)));
theta_right = theta(floor(length(theta)/2)+2:end);

param.method.src_limits = [];
for idx = 1:length(theta_left)
  param.method.src_limits{idx}     = [min(theta) max(theta)];
end

param.method.theta_guard            = 1.5/180*pi;
param.method.nb_nd.init             = 'ap';
param.method.wb_td.init             = 'ap';
param.method.wb_td.widening_factor  = 1;
param.method.wb_fd.init             = 'ap';
param.method.wb_fd.filter_banks     = 1;

k = 2*pi/(lambda/2);

% Beam steering angle
steer_ang = 0*pi/180;

% Complex transmit weights
y_pc = param.src.phase_center(2,:).';
z_pc = param.src.phase_center(3,:).';

% If you steer left/right set the right/right sensors' weights to 0 because
% the left/right sensors pattern tend to be biased twords left/right.
%   tx_weights    = [hanning(4).',0 0 0 0].';
tx_weights    = hanning(Nc);    %% CHECK WHY HANNING WINDOW PRODUCES WRONG RESULTS
if ~exist('tx_weights','var')
  tx_weights = ones(Nc,1);
end

steering_delay = (y_pc*sin(steer_ang) - z_pc*cos(steer_ang))/c;
w = tx_weights.*exp(-1i*4*pi*fc*steering_delay);

steering_mtx = @(DoA) ((1/sqrt(length(y_pc)))*exp(1i*(y_pc*k*sin(DoA) - z_pc*k*cos(DoA))));
comp_tx_weight = w.'*steering_mtx(theta.');
param.src.tx_weights = comp_tx_weight;
param.src.theta = theta*180/pi;

% Set the limits of the DoAs to within 3dB of the antenna beampattern
theta_RP = linspace(-90,90,2048)*pi/180;
RP = abs(w.'*steering_mtx(theta_RP)).^2;
RP = RP./max(RP);
RP_dB = 10*log10(RP);

idxs_3dB = find(RP_dB>=-3);
RP_3dB = RP_dB(idxs_3dB);
theta_3dB = theta_RP(idxs_3dB);

plot_doa_lims = [min(theta_3dB) max(theta_3dB)]*180/pi;
%   plot_doa_lims = [-2 33];

% ------------------------------------------------------------------------
%                           Generate array errors
% ------------------------------------------------------------------------
error_ypc      = [0 0.001 -0.003 0 0.009 0.002 -0.001]'*lambda;
error_zpc      = [0 0.001 0.001 -0.002 0.001 0.002 0.001]'*lambda;
error_phase    = [0 0 0 0 0 0 0]';
error_g_s      = [0 0.8 1 0.9 1 0.8 1]';%[0 1 -0.1 0.5 0.1 -1 0.6]';
error_g_p      = [0 0 15 0 -5 0 10]'*pi/180;
error_g_offset = [0 -0.1 3 -2 0 0.1 0.01]';

Err(:,1) = error_ypc./lambda;
Err(:,2) = error_zpc./lambda;
Err(:,3) = error_phase;
Err(:,4) = error_g_s;
Err(:,5) = error_g_p;
Err(:,6) = error_g_offset;

% Array calibration errors
error_params.error_ypc      = error_ypc ;
error_params.error_zpc      = error_zpc ;
error_params.error_phase    = error_phase;
error_params.error_g_s      = error_g_s;
error_params.error_g_p      = error_g_p;
error_params.error_g_offset = error_g_offset;

param.error_params          = error_params; % With array errors
%   param.error_params          = [];           % No array errors

% ------------------------------------------------------------------------
%                           Run DoA estimator
% ------------------------------------------------------------------------
rmse_tmp = NaN(length(theta_left)+1,2);

SNR = 10;
Nsnaps = 3*21;
Nruns = 100;
% DoA = 0
param.monte.SNR   = SNR*ones(1,1);
param.monte.Nsnap = Nsnaps*ones(1,1);
param.monte.runs  = Nruns;
param.monte.random_seed_offset = 0;

param.monte.DOA   = [0]*180/pi;
comp_tx_weight = w.'*steering_mtx(param.monte.DOA*pi/180);
param.src.tx_weights = comp_tx_weight;

results = sim.doa(param);

err = abs(param.monte.DOA*pi/180 - results.theta_est{param.method.list});
%   sigma = std(err);
%   err(err > 3*sigma) = [];
RMSE = sqrt(nanmean(abs(err).^2));

%   RMSE = sim.doa_rmse(param,results);
rmse_tmp(1,1) = squeeze(RMSE);
rmse_tmp(1,2) = NaN;

est_doa0 = mean(results.theta_est{param.method.list}(:,1))*180/pi;

% All other DoAs
param.monte.SNR   = SNR*ones(1,2);
param.monte.Nsnap = Nsnaps*ones(1,2);
param.monte.runs  = Nruns;
param.monte.random_seed_offset = 0;

est_doa1 = [];
est_doa2 = [];
for doa_idx = 2:length(theta_left)+1
  param.monte.DOA   = [theta_left(doa_idx-1) theta_right(doa_idx-1)]*180/pi;
  
  comp_tx_weight = w.'*steering_mtx(param.monte.DOA*pi/180);
  param.src.tx_weights = comp_tx_weight;
  
  results = sim.doa(param);
  
  if 0
    % Ignore DoAs outside the 3dB beamwidth
    param.monte.DOA(param.monte.DOA < plot_doa_lims(1)) = NaN;
    param.monte.DOA(param.monte.DOA > plot_doa_lims(2)) = NaN;
  end
  
  %     RMSE = sim.doa_rmse(param,results);
  
  for src_idx = 1:length(param.monte.DOA)
    err = param.monte.DOA(src_idx)*pi/180 - results.theta_est{param.method.list}(:,src_idx);
    %       sigma = std(err);
    %       err(err > 3*sigma) = [];
    RMSE(src_idx) = sqrt(nanmean(abs(err).^2));
  end
  
  rmse_tmp(doa_idx,:) = squeeze(RMSE);
  est_doa1(doa_idx-1) = mean(results.theta_est{param.method.list}(:,1))*180/pi;
  est_doa2(doa_idx-1) = mean(results.theta_est{param.method.list}(:,2))*180/pi;
end
est_doa = [fliplr(est_doa1), est_doa0, est_doa2].';

if 0
  % Ignore DoAs outside the 3dB beamwidth
  est_doa(est_doa<plot_doa_lims(1)) = NaN;
  est_doa(est_doa>plot_doa_lims(2)) = NaN;
end

rmse_all = [flipud(rmse_tmp(:,1));rmse_tmp(2:end,2)];
rmse_all(isnan(est_doa)) = NaN;

true_doa = theta*180/pi;

if 0
  threshold = 2 *std(rmse_all(~isnan(rmse_all)));
  rmse_all(rmse_all>threshold) = NaN;
end
% ------------------------------------------------------------------------
%                           Plot
% ------------------------------------------------------------------------
figure(14);
scatter(theta*180/pi,rmse_all,20,'fill','r');
%   xlim([plot_doa_lims(1)-5 plot_doa_lims(2)+5])
xlim([-25 +25])
xlabel('True DoA (deg.)')
ylabel('RMSE (deg.)')
Title = sprintf('N_c = %d, SNR = %d dB, Nsnaps = %d, and Nruns = %d',Nc, param.monte.SNR(1), param.monte.Nsnap(1), Nruns);
title(Title)
grid on
%   legend('Without errors','With errors', 'Location','best')
%   hold on

figure(15);
scatter(true_doa,est_doa,20,'fill','r')
%   xlim([plot_doa_lims(1)-5 plot_doa_lims(2)+5])
%   ylim([plot_doa_lims(1)-5 plot_doa_lims(2)+5])
xlim([-25 +25])
ylim([-25 +25])
xlabel('True DoA (deg.)')
ylabel('Estimated DoA (deg.)')
Title = sprintf('N_c = %d, SNR = %d dB, Nsnaps = %d, and Nruns = %d',Nc, param.monte.SNR(1), param.monte.Nsnap(1), Nruns);
title(Title)
grid on
%   legend('Without errors','With errors', 'Location','best')
%   hold on

