% script sim.doa_example
%
% Example setup scripts for the sim.doa function. To run, enable one
% of the if {0|1} blocks below and then run the script.
%
% Author: John Paden, Theresa Stumpf
%
% See also: doa.m


%% Array calibration
% =========================================================================
param = [];
physical_constants

% -----------------------------------------------------------------------
% Source parameters
% -----------------------------------------------------------------------
Nc                = 7;
Nruns             = 1;
param.src.fc      = 195e9;%195e6;%320e6;
param.src.BW      = 10e6;%1e6;%32e6;
param.src.fs      = param.src.BW;
param.src.f0      = param.src.fc-param.src.BW/2;
param.src.f1      = param.src.fc+param.src.BW/2;
param.src.ft_wind = @(N) hanning(N);
param.src.SNR     = 30*ones(1,2);
param.src.Nsnap   = 5*21*ones(1,1);
%   param.src.lever_arm.fh     = @sim.lever_arm_example;
%   param.src.lever_arm.args   = {[],1,[1:Nc],[0; c/param.fc/2; 0]};

param.method.wb_td.widening_factor  = 1;
param.method.wb_fd.filter_banks     = 1;

% -----------------------------------------------------------------------
% Beam steering and 3dB beamwidth
% -----------------------------------------------------------------------
theta = [-1.2154, -1.0654, -0.9484, -0.8481, -0.7580, -0.6751, -0.5974, -0.5236, -0.4528, -0.3844, -0.3178, -0.2527, ...
  -0.1886, -0.1253, -0.0625, 0, 0.0625, 0.1253, 0.1886, 0.2527, 0.3178, 0.3844, 0.4528, 0.5236, 0.5974, 0.6751, ...
  0.7580, 0.8481, 0.9484, 1.0654, 1.2154]';
%   theta = [-75:3:75]'*pi/180;
param.src.theta = theta;
lambda = c/param.src.fc;
d_y    = lambda/2;

phase_center      =  zeros(3,Nc);
phase_center(2,:) = -d_y/2*(0:Nc-1);%-d_y/2*(0:Nc-1);

y_pc = phase_center(2,:).';
z_pc = phase_center(3,:).';

param.src.y_pc       = y_pc;
param.src.z_pc       = z_pc;
k = 2*pi/(lambda/2);

% Steering martix function handle
A = @(theta) sim.steering_mtx(theta,param);

% Beam steering angle
steer_ang = 0*pi/180;

% Complex transmit weights
tx_weights     = hanning(Nc);
%   tx_weights     = [hanning(4).',0 0 0 0].';
steering_delay = (y_pc*sin(steer_ang) - z_pc*cos(steer_ang))/c;
w = tx_weights.*exp(-1i*4*pi*param.src.fc*steering_delay);

% Set the limits of the DoAs to within 3dB of the antenna beampattern
theta_RP = linspace(-90,90,2048)*pi/180;
RP = abs(w.'*A(theta_RP)).^2;
RP = RP./max(RP);
RP_dB = 10*log10(RP);

idxs_3dB  = find(RP_dB>=-3);
RP_3dB    = RP_dB(idxs_3dB);
theta_3dB = theta_RP(idxs_3dB);

beam_doa_lims = [min(theta_3dB) max(theta_3dB)];

%   actual_doa = theta(theta>beam_doa_lims(1) & theta<beam_doa_lims(2));
actual_doa = theta;
comp_tx_weight = w.'*A(actual_doa.');
%   param.src.tx_weights = comp_tx_weight;

% -----------------------------------------------------------------------
% Array calibration: generate errors
% -----------------------------------------------------------------------
% Error bounds
error_bounds = [-0.5 0.5;-0.5 0.5;0 0;-1 1;-15*pi/180 15*pi/180;-sqrt(10) sqrt(10)];

if 0
  % Rando error generation (all errors must be within the bounds)
  ac_error_gen_param.Nc           = Nc;
  ac_error_gen_param.ref_sensor   = ceil(Nc/2);
  ac_error_gen_param.error_bounds = error_bounds;
  
  ac_errors      = ac_error_gen(ac_error_gen_param);
  error_ypc      = ac_errors.error_ypc * lambda;
  error_zpc      = ac_errors.error_zpc * lambda;
  error_phase    = ac_errors.error_phase;
  error_g_s      = ac_errors.error_g_s;
  error_g_p      = ac_errors.error_g_p;
  error_g_offset = ac_errors.error_g_offset;
elseif 1
  % Manually enter error values (all errors must be within the bounds)
  error_ypc      = [0 0.001 -0.003 0 0.009 0.002 -0.001]'*lambda;
  error_zpc      = [0 0.001 0.001 -0.002 0.001 0.002 0.001]'*lambda;
  error_phase    = [0 0 0 0 0 0 0]';
  error_g_s      = [0 0.8 1 0.9 1 0.8 1]';%[0 1 -0.1 0.5 0.1 -1 0.6]';
  error_g_p      = [0 0 15 0 -5 0 10]'*pi/180;
  error_g_offset = [0 -0.1 3 -2 0 0.1 0.01]';
end
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

param.error_params          = error_params;

ac_est_error = [];
%   for run_idx = 1:Nruns

% -----------------------------------------------------------------------
% Generate data covariance matrix
% -----------------------------------------------------------------------
% Call data generation function
for doa_idx = 0:ceil(length(actual_doa)/2)-1%2:floor(length(actual_doa)/2)+1
  param.src.DOAs       = actual_doa([1+doa_idx,end-doa_idx]).'*180/pi;
  param.src.tx_weights = comp_tx_weight([1+doa_idx,end-doa_idx]);
  if param.src.DOAs(1)==param.src.DOAs(2)
    param.src.DOAs       = param.src.DOAs(1);
    param.src.tx_weights = param.src.tx_weights(1);
  end
  [~,DCM{doa_idx+1}] = sim.doa_wideband_data(param);
  
end

param.error_params = [];

% -----------------------------------------------------------------------
% Array calibration: run the optimizer
% -----------------------------------------------------------------------

extra_error_params.error_ypc      = error_ypc ;
extra_error_params.error_zpc      = error_zpc ;
extra_error_params.error_phase    = error_phase;
extra_error_params.error_g_s      = error_g_s;
extra_error_params.error_g_p      = error_g_p;
extra_error_params.error_g_offset = error_g_offset;

param.extra_error_params = extra_error_params;

ac_cost_params.actual_doa = actual_doa;
ac_cost_params.y_pc       = y_pc;
ac_cost_params.z_pc       = z_pc;
ac_cost_params.fc         = param.src.fc;
ac_cost_params.DCM        = DCM;

LB = repmat(error_bounds(:,1).',[Nc 1]);
UB = repmat(error_bounds(:,2).',[Nc 1]);

LB(1,:) = 0;
UB(1,:) = 0;

LB = LB(:);
UB = UB(:);

initial_ac = 0*ones(size(LB)); % Nc*6 matrix

tic
% Run the optimizer
% -------------------------------------------------------------------
if 1
  % Local solver (fmincon solver)
  options =  optimoptions(@fmincon,'TolX',1e-6,'TolFun',1e-6,'MaxIter',10^8,'MaxFunEvals',10^8, ...
    'PlotFcn',{@optimplotfval,@optimplotstepsize},'Algorithm','sqp');
  [ac_est_error,~,exitflag] = fmincon(@(est_errors) sim.array_calibration_cost_1D(est_errors,ac_cost_params),...
    initial_ac, [],[],[],[],LB,UB,[],options);
  exitflag = exitflag
elseif 0
  % Local solver (patternsearch solver) .. Too slow for TolMesh>1e-2
  options =  psoptimset('TolX',1e-6,'TolFun',1e-6,'TolMesh',1e-4,'InitialMeshSize',0.001,'MaxIter',10^8,'MaxFunEvals',10^8,...
    'PlotFcn',@psplotbestf);
  % options =  psoptimset('TolX',1e-15,'InitialMeshSize',0.001,'MaxIter',10^100,'MaxFunEvals',10^100,'TolMesh',1e-15);
  ac_est_error = patternsearch(@(est_errors) sim.array_calibration_cost_1D(est_errors,ac_cost_params),...
    initial_ac, [],[],[],[],LB,UB,[],options);
  
elseif 0
  % Global solver
  % Sometimes it fails and recover .. Wast of time
  gs = GlobalSearch;
  objFun = @(est_errors) sim.array_calibration_cost_1D(est_errors,ac_cost_params);
  options = optimoptions(@fmincon,'TolX',1e-6,'TolFun',1e-3,'MaxIter',10^6,'MaxFunEvals',10^6, ...
    'PlotFcn',{@optimplotfval,@optimplotstepsize});
  problem = createOptimProblem('fmincon','x0',initial_ac,'objective',objFun,'lb',LB,'ub',UB,'options',options);
  ac_est_error = run(gs,problem);
elseif 0
  
  % Alternating Projection (AP)-like implementation
  % -----------------------------------------------
  % TO BE DONE LATER ...
  
end

ac_est_error = reshape(ac_est_error,[Nc 6]);
% Measure the accuracy of the result
% ----------------------------------
% 1) RMSE
sprintf('\n')
rmse           = sqrt(mean(abs(ac_est_error-Err).^2,1))
% 2) Mean of the absolute error
sprintf('\n')
mean_abs_error = mean(abs(ac_est_error-Err),1)
sprintf('\n')
% 3) Mean abs error relative to the maximum error
relative_max   = mean_abs_error./max(abs(Err))

%   end
toc
%     ac_est_error        = mean(ac_est_error,2);
%     ac_est_error        = reshape(ac_est_error,[Nc 6]);

est_error_ypc       = ac_est_error(:,1);
est_error_zpc       = ac_est_error(:,2);
est_error_phase     = ac_est_error(:,3);
est_error_g_s       = ac_est_error(:,4);
est_error_g_p       = ac_est_error(:,5);
est_error_g_offset  = ac_est_error(:,6);

% Determine the difference between actual and estimated errors
%     ac_error_diff(trial_idx,:,1) = abs(error_ypc./lambda-est_error_ypc);
%     ac_error_diff(trial_idx,:,2) = abs(error_zpc./lambda-est_error_zpc);
%     ac_error_diff(trial_idx,:,3) = abs(error_phase-est_error_phase);
%     ac_error_diff(trial_idx,:,4) = abs(error_g_s-est_error_g_s);
%     ac_error_diff(trial_idx,:,5) = abs(error_g_p-est_error_g_p);
%     ac_error_diff(trial_idx,:,6) = abs(error_g_offset-est_error_g_offset);

% -----------------------------------------------------------------------
%           Plot the gain and phase deviation patterns for each sensor
% -----------------------------------------------------------------------

if 1
  %  theta_RP = linspace(-90,90,2048)*pi/180;
  sv_params.src.y_pc = y_pc;
  sv_params.src.z_pc = z_pc;
  sv_params.src.fc   = param.src.fc;
  extra_error_params.error_ypc      = error_ypc;
  extra_error_params.error_zpc      = error_zpc;
  extra_error_params.error_phase    = error_phase;
  extra_error_params.error_g_s      = error_g_s;
  extra_error_params.error_g_p      = error_g_p;
  extra_error_params.error_g_offset = error_g_offset;
  
  % Calculate the 3dB limits of the beampattern (to plot)
  sv_params.extra_error_params = [];
  A = @(theta,param)sim.steering_mtx(theta,param);
  SV = A(theta_RP,sv_params);
  rad_pattern = abs(w.'*SV).^2;
  rad_pattern = rad_pattern./max(rad_pattern);
  gain_pattern_dB = 10*log10(rad_pattern);
  
  figure(100);clf
  plot(theta_RP*180/pi,gain_pattern_dB,'r')
  xlabel('\theta^\circ')
  ylabel('Power (dB)')
  title('Array radiation pattern')
  xlim([beam_doa_lims(1) beam_doa_lims(2)]*180/pi)
  grid on
  
  sv_params.extra_error_params = extra_error_params;
  A = @(theta,param)sim.steering_mtx(theta,sv_params);
  SV = A(theta_RP,sv_params);
  
  % Plot sensors gain pattern
  figure(101);clf
  hold on;
  for chan_idx = 1:Nc
    chan_resp = SV(chan_idx,:);
    chan_gain = abs(chan_resp).^2;
    chan_gain = chan_gain./max(chan_gain);
    chan_gain_dB = 10*log10(chan_gain);
    plot(theta_RP*180/pi,chan_gain_dB)
  end
  xlabel('\theta^\circ')
  ylabel('Power (dB)')
  title('Sensors gain pattern (sensor 1 is the reference)')
  xlim([-30 +30])
  grid on
  legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
  
  % Plot sensors phase patten (with phase errors)
  figure(102);clf
  hold on
  for chan_idx = 1:Nc
    chan_resp = SV(chan_idx,:);
    chan_phase = angle(chan_resp)*180/pi;
    plot(theta_RP*180/pi,chan_phase)
  end
  xlabel('\theta^\circ')
  ylabel('Phase (deg.)')
  title('Sensors phase pattern (sensor 1 is the reference)')
  xlim([-9.5 +9.5])
  grid on
  legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
  
  % Plot phase deiviation (ideal+error-ideal) -- Only phase errors
  sv_params.extra_error_params = [];
  A = @(theta,param)sim.steering_mtx(theta,sv_params);
  SV_ideal = A(theta_RP,sv_params);
  SV_phase_error_only = conj(SV_ideal) .* SV;
  
  figure(103);clf
  hold on
  for chan_idx = 1:Nc
    chan_resp = SV_phase_error_only(chan_idx,:);
    chan_phase = angle(chan_resp)*180/pi;
    plot(theta_RP*180/pi,chan_phase)
  end
  xlabel('\theta^\circ')
  ylabel('Phase (deg.)')
  title('Sensors phase deviation (sensor 1 is the reference)')
  xlim([-30 +30])
  grid on
  legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
  
end

% This is to test for a good gain error model
% DoA = [-pi/2:0.1:pi/2]*180/pi;
% g_phase = pi/10*180/pi;
%
% f=@(DoA,g_phase) 10.^((((sind(DoA)-sind(g_phase)).^2+2)./20));
% % f=@(DoA,alpha_phase) ((sind(DoA)-sind(g_phase)).^2+2);
% % f=@(DoA,alpha_phase) (1+(sind(DoA)-sind(g_phase)).^2+2);
% % f=@(DoA,alpha_phase) (1+(sind(DoA-g_phase)).^2+2);
%
%
%  figure(47);plot(DoA,f(DoA,g_phase))
%
%  xlabel('theta (deg)')
%  ylabel('gain')
%
%  title('(10.^((((sind(DoA)-sind(alpha_phase)).^2+2)./20))','Interpreter','Latex')
% grid on
%
% % ------------------------
% DoA = [-pi/2:0.1:pi/2];
% g_phase = pi/10;
% %
% % f=@(DoA,g_phase) 10.^((((sinc(DoA)-sinc(g_phase)).^2+2)./20));
% % f=@(DoA,g_phase) ((sinc(DoA)-sinc(g_phase)).^2+2);
% % f=@(DoA,g_phase) (1+(sinc(DoA)-sinc(g_phase)).^2+2);
% f=@(DoA,g_phase) (1+(sinc(DoA-g_phase)).^2+2);
%
%
%  figure(13);plot(DoA*180/pi,f(DoA,g_phase))
%
%  xlabel('theta (deg)')
%  ylabel('gain')
%
%  title('(1+(sinc(DoA-g_phase)).^2+2)','Interpreter','Latex')
% grid on
