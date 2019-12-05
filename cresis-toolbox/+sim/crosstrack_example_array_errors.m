% script sim.crosstrack_example_array_errors
%
% Example array calibration errors simulation. Calls crosstrack.m and
% crosstrack_rmse.m.
%
% Author: Mohanad Al-Ibadi, Sravya Athinarapu, Sean Holloway, John Paden,
% Theresa Stumpf
%
% See also: sim.crosstrack.m, sim.crosstrack_data.m,
% sim.crosstrack_example_*.m, sim.crosstrack_rmse.m

%%    Array Calibration: Effect of array errors on estimated DoAs
% =======================================================================
physical_constants;
param = [];

% Debug level of 3 causes this function to stop early and output
% simulated data, simulation parameters, and array processing results
% in a "results" structure.
param.debug_level = 3;

% -----------------------------------------------------------------------
%                          Source parameters
% -----------------------------------------------------------------------
fc = 195e9;
BW = 32e6;
fs = BW;
lambda = c/fc;
%   targets_doa = [-85:5:85].'*pi/180;
%   targets_doa = [-1.2154, -1.0654, -0.9484, -0.8481, -0.7580, -0.6751, -0.5974, -0.5236, -0.4528, -0.3844, -0.3178, -0.2527, ...
%     -0.1886, -0.1253, -0.0625, 0, 0.0625, 0.1253, 0.1886, 0.2527, 0.3178, 0.3844, 0.4528, 0.5236, 0.5974, 0.6751, ...
%     0.7580, 0.8481, 0.9484, 1.0654, 1.2154]';
targets_doa = [-1.0291 -0.9742 -0.9075 -0.8366 -0.7684 -0.6867 -0.6303 -0.5606 -0.4861 -0.4190 -0.3480 -0.2773 -0.2127 -0.1398 ...
  -0.0679 0.0050 0.0671 0.1381 0.2086 0.2804 0.3497 0.4198 0.4894 0.5502 0.6294 0.7036 0.7491 0.8367 0.9089 0.9836 1.0389 1.0913].';
% targets_doa = [-0.4407   -0.3534   -0.2651   -0.1840   -0.0995    0.0034    0.0581    0.1609    0.2494    0.3394    0.4244].';
%   targets_doa(abs(targets_doa)>31*pi/180) = [];
flight_height  = 1500;
range_vec      = flight_height./cos(targets_doa);
max_range      = max(range_vec);
min_range      = min(range_vec);

param.src.f0                      = fc-BW/2;
param.src.f1                      = fc+BW/2;
param.src.t0                      = 2*(min_range-500)/c;
param.src.t1                      = 2*(max_range+500)/c;
param.src.ft_func                 = @(t) tukeywin_cont(t * BW);
param.src.ft_wind                 = @(N) hanning(N);
param.src.lever_arm.fh            = @sim.lever_arm_example;
% Nsep: number of lambda/4 steps between sensors
Nsep = 1;
% Nc: number of sensors
Nc = 7;
% Arguments for a linear array in y-dimension
% param.src.lever_arm.fh_args: Arguments for a linear array in y-dimension
param.src.lever_arm.fh_args       = {[], 1, 1:Nc, [0; Nsep*c/fc/2; 0]};
%   param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
param.src.noise_power             = zeros(1,Nc);


% DOA method parameters
param.method.list                   = [7];

[phase_center] = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
y_pc = phase_center(2,:).';
z_pc =0.5*lambda  + phase_center(3,:).'; % Make it a multiple of lambda/2 for best results

%   y_pc = [1.4126    1.0409    0.6891    0.3111   -0.0669   -0.4190   -0.7909].';
%   z_pc = [0.0612    0.0368    0.0125   -0.0031    0.0093    0.0306    0.0518].';

% -----------------------------------------------------------------------
%                        Simulation Runs Setup
% -----------------------------------------------------------------------
% Cross track monte carlo setup
param.monte.target_func = @sim.surface_gen;
param.monte.runs = 1;
param.monte.random_seed_offset = 0;

% Target surface parameters
surf_param = [];
surf_param.z.mean = -flight_height;
surf_param.z.rms_height = 0;
surf_param.z.corr_length_x = 400;
surf_param.z.corr_length_y = 400;
surf_param.rcs_in.mean = 0;
surf_param.rcs_in.var = 1e4;
surf_param.dy = 10;
%   surf_param.y_range = [-2500 2500];
surf_param.dx = 10;
surf_param.x_range = [-2500 2500];
%   surf_param.x = [-500:1:500]; % Defined later
%   surf_param.y = [-1500:20:1500].';
surf_param.y = range_vec .* sin(targets_doa);

% y_range should >= maximum y (i.e. targets shoulb be inside the imaged
% swath)
surf_param.y_range = [-1.5*max(surf_param.y) 1.5*max(surf_param.y)];
if max(surf_param.y_range) == 0
  surf_param.y_range = [-2500 2500];
end
%   param.monte.target_param{1} = surf_param;
% -----------------------------------------------------------------------
%                   Array Processing parameters
% -----------------------------------------------------------------------
array_param = [];

%% NN: total length of the sensor array (in lambda/4 units)
NN = Nc*Nsep;
% lambda: wavelength
lambda = c/fc;
% k: wavenumber
k = 2*pi/(lambda/2);
% My: over sampling factor
My = 4;
% dy: phase center spacing
dy = Nsep*lambda/4;
% dky and ky: y-component of wavenumber (spacing and axis)
dky = 2*pi / (Nc*dy) / My;
ky = dky * ifftshift(-floor(My*Nc/2) : floor((My*Nc-1)/2));
% theta: theta values associated with ky axis
theta = fftshift(asin(ky/k));
array_param.Nsv = {'theta', asin(ky/k)};

array_param.sv_fh = @array_proc_sv;

array_param.dbin = 1;
array_param.dline = 1;

array_param.bin_rng = 0;
array_param.line_rng = -50:50;

Nsrc = 2;
array_param.Nsrc = Nsrc;

array_param.init = 'ap';
array_param.doa_theta_guard = (max(theta)-min(theta))/(4*Nc);

array_param.Nsubband = 1;
dt = 1/fs;
array_param.imp_resp.time_vec = -3*array_param.Nsubband*dt : dt/8 : 3*array_param.Nsubband*dt;
BW = abs(param.src.f1 - param.src.f0);
array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / BW);

for idx = 1:array_param.Nsrc
  array_param.doa_constraints(idx).method = 'fixed';
  array_param.doa_constraints(idx).init_src_limits = [-90 90];%[min(theta) max(theta)]*180/pi;
  array_param.doa_constraints(idx).src_limits = [-90 90];%[min(theta) max(theta)]*180/pi;
end

param.array_param = array_param;

N_skipped_rlines = 1;
N_reqd_rlines    = 1;
surf_param.x = [-750:N_skipped_rlines:-750+N_skipped_rlines*(N_reqd_rlines+2*max(array_param.line_rng))-1];
param.monte.target_param{1} = surf_param;

%   clear array_param;
if 0
  % Debug/test code
  surf_model = param.monte.target_func(param.monte.target_param{1});
  var(surf_model.rcs(:))
  std(surf_model.z(:))
  imagesc(surf_model.dem)
  surf(surf_model.x, surf_model.y, surf_model.z)
  xlabel('x (m), along-track snapshots')
  ylabel('y (m), cross-track')
  return
end

% ------------------------------------------------------------------------
%                           Generate array errors
% ------------------------------------------------------------------------
% Location errors and error bounds are have length units (i.e. meters)
%     error_bounds = [-0.5 0.5;-0.5 0.5;-15*pi/180 15*pi/180;-1 1;-15*pi/180 15*pi/180;-sqrt(10) sqrt(10)];
error_bounds =  [-0.2 0.2; -0.2 0.2; -1*pi 1*pi; -10 10; -1*pi 1*pi; -10 10];
%    error_bounds =  [-0.2 0.2; -0.2 0.2; -1*pi 1*pi; 0 0; 0 0; 0 0];
%   error_bounds = [-0.2    0.2;  % Location errors are in meters
%                   -0.2    0.2;
%                   -1*pi   1*pi;
%                   -10      10;
%                   -1*pi   1*pi;
%                   -10     10];
%
error_ypc      = [-0.1943    0.0452   -0.1972         0   -0.1958   -0.1943    0.0875]';%* lambda; % In meters units
error_zpc      = [-0.1924    0.0452    0.1944         0    0.1965   -0.0128   -0.1859]';%* lambda; % In meters units
error_phase    = [-1.1228    0.6236    1.6716         0    1.6941    0.1720   -1.0669]';
error_g_s      = [-3.7225    4.5863    4.0187         0    9.8245    4.0963    3.7311]';
error_g_p      = [-1.4270   -1.6135   -1.7050         0   -0.0157   -1.9984   -1.5337]';
error_g_offset = [-9.8911    6.6751    8.8213         0   -9.7071    9.6572    9.7167]';

% y_pc and z_pc errors are generated as percentage fo the actual y_pc
% and z_pc values. This guarantees that the errors are sensible.
y_pc_err_percentage = [1 5 0 4 9 2 8].' ./100;
z_pc_err_percentage = [1  7  1  1  3  5  2].' ./100;
%     error_ypc      = y_pc .* y_pc_err_percentage;%  [0.02 0.01 -0.03 0 0.009 0.002 -0.001]';% * lambda; % In meters units
%     error_zpc      = z_pc .* z_pc_err_percentage;% [0 0.001 0.001 -0.02 0.005 0.01 0.001]';% *lambda ; % In meters units
%     error_phase    = [5 1 5 0 -15 0.2 10]'*pi/180;%[0 0 0 0 0 0 0]';
%     error_g_s      = [0.2 0.8 1 0.9 1 0.8 1]';
%     error_g_p      = [0 0 15 0 -5 0 10]'*pi/180;
%     error_g_offset = [-1 -0.1 3 -2 0 0.1 8]';

Err(:,1) = error_ypc;%./lambda;
Err(:,2) = error_zpc;%./lambda;
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

% Sensors gain model
% ------------------
param.array_gain_model_fh = @(x) (exp(-x./40));
%       sigma = 1;
%       param.array_gain_model_fh = @(x) (1/sqrt(2*pi*sigma^2)*exp(-((x-1).^2)./(2*sigma^2)));
%     param.array_gain_model_fh = @(x) (10^(x./20));
%     param.array_gain_model_fh = @(x) (1+x);

% -----------------------------------------------------------------------
%                       Run the simulation
% -----------------------------------------------------------------------
results = crosstrack(param);
%   results = sim.crosstrack(param);

sim_data        = squeeze(results.sim_data{1}); % Nt*Nx*Nc data matrix
est_doa         = results.tomo.doa;
actual_doa_cell = results.actual_doa;
array_param     = results.array_param;

Nt = size(sim_data,1); %length(array_param.bins);
Nx = size(sim_data,2); %length(array_param.lines);

% Convert actual_doa from cell array into a matrix (maximum Nsrc targets
% per range-bin).
actual_doa = NaN(Nt,Nsrc,Nx);
lines = array_param.lines(1):array_param.dline:array_param.lines(end);
bins = array_param.bins(1):array_param.dbin:array_param.bins(end);
for lineIdx = 1:length(lines) %length(array_param.lines)
  lineIdx_idx = lines(lineIdx);
  for binIdx = 1:length(bins)  %length(array_param.bins)
    binIdx_idx = bins(binIdx);
    if ~isempty(actual_doa_cell{binIdx,lineIdx})
      doa_tmp = actual_doa_cell{binIdx,lineIdx};
      if length(doa_tmp)<Nsrc
        doa = NaN(Nsrc,1);
        if sign(doa_tmp) <= 0
          doa(1:length(doa_tmp)) = doa_tmp;
        else
          doa(length(doa_tmp)+1:end) = doa_tmp;
        end
      else
        doa = doa_tmp(1:Nsrc);
      end
      actual_doa(binIdx_idx,1:Nsrc,lineIdx_idx) = doa; % Nt*Nsrc*Nx
    end
  end
end

%   binIdx_vec = [];
%   for lineIdx = 1:length(lines)
%     for binIdx = 1:length(bins)
%       if ~isempty(actual_doa_cell{binIdx,lineIdx})
%         binIdx_vec(end+1,lineIdx) = binIdx;
%       end
%     end
%   end
%

if 0
  % Debug: plot the actual vs the estimated surfaces
  figure(99);clf
  hold on
  h1 = plot(squeeze(actual_doa(:,1,1))*180/pi,1:Nt,'*r');
  h2 = plot(squeeze(actual_doa(:,2,1))*180/pi,1:Nt,'*r');
  
  h3 = plot(squeeze(est_doa(:,1,1))*180/pi,1:Nt,'*k');
  h4 = plot(squeeze(est_doa(:,2,1))*180/pi,1:Nt,'*k');
  
  set(gca,'YDir','reverse')
  xlim([-30 30])
  ylim([100 180])
  %     ylim([1 size(est_doa,1)])
  grid on
  grid minor
  
  xlabel('\theta^\circ')
  ylabel('Range-bin index')
  title('True vs estimated surfaces')
  legend([h1 h3],'Actual DoA','Estimated DoA','Location','northeast')
  %     legend([h1 h3 h5],'Actual DoA','Estimated DoA: without array errors','Estimated DoA: with array errors','Location','best')
end

if 0
  % Plot true vs estimated DoA with and without errors
  figure(50);clf
  hold on
  h1 = plot(squeeze(actual_doa(:,1,1))*180/pi,squeeze(est_doa(:,1,1))*180/pi,'*r');
  h2 = plot(squeeze(actual_doa(:,2,1))*180/pi,squeeze(est_doa(:,2,1))*180/pi,'*r');
  xlim([-25 25])
  ylim([-25 25])
  grid on
  
  xlabel('True DoA (deg.)')
  ylabel('Estimated DoA (deg.)')
  Title = sprintf('N_c = %1.0f, rcs.var = %3.0f, Nx = %3.0f, Nt = %3.0f',Nc,surf_param.rcs_in.var,Nx, Nt);
  title(Title)
  %     legend([h1 h3],'Without array errors','With array errors','Location','northwest')
end

if 0
  % Plot RMSE (average over range-lines)
  doa_error = abs(actual_doa-est_doa);
  sigma_doa = sqrt(nanvar(doa_error(:)));
  doa_error(doa_error>=sigma_doa) = NaN;
  
  rmse = sqrt(nanmean(doa_error.^2,3));
  doa = squeeze(actual_doa(:,:,1));
  
  figure(51);clf
  hold on
  h2 = plot(doa2*180/pi,rmse*180/pi,'*r');
  xlim([-25 25])
  grid on
  
  xlabel('True DoA (deg.)')
  ylabel('RMSE (deg.)')
  Title = sprintf('N_c = %1.0f, rcs.var = %3.0f, Nx = %3.0f, Nt = %3.0f',Nc,surf_param.rcs_in.var,Nx, Nt);
  title(Title)
  %     legend('Without array errors','With array errors','Location','best')
end
% -----------------------------------------------------------------------
%                         Cost function parameters
% -----------------------------------------------------------------------
% extra_error_params is used in steering_mtx function, which is called
% from array_calibration_cost_2D
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
ac_cost_params.fc         = fc;
ac_cost_params.sim_data   = sim_data;
ac_cost_params.array_param = array_param;
ac_cost_params.array_gain_model_fh = param.array_gain_model_fh;

tic
% -----------------------------------------------------------------------
%                          Run the optimizer
% -----------------------------------------------------------------------
LB = repmat(error_bounds(:,1).',[Nc 1]);
UB = repmat(error_bounds(:,2).',[Nc 1]);

%   LB(1,:) = 0;
%   UB(1,:) = 0;
LB(ceil(Nc/2),:) = 0;
UB(ceil(Nc/2),:) = 0;

LB = LB(:);
UB = UB(:);
initial_ac = 0*ones(size(LB)); % Nc*6 matrix

[nadir_doa,nadir_doa_idx] = min(abs(targets_doa));

param.fc = fc;
param.nadir_y_pc = y_pc;
param.nadir_z_pc = z_pc;
param.nadir_doa  = nadir_doa;
param.doa = targets_doa.';
array_calib_nonlcon_fh = @(est_errors) array_calib_nonlcon(est_errors,param);

if 1
  % Local solver (fmincon solver)
  options =  optimoptions(@fmincon,'TolX',1e-6,'TolFun',1e-6,'MaxIter',10^8,'MaxFunEvals',10^8, ...
    'PlotFcn',{@optimplotfval,@optimplotstepsize});%,'Algorithm','sqp');
  [ac_est_error, ~, exitflag] = fmincon(@(est_errors) array_calibration_cost_2D(est_errors,ac_cost_params),...
    initial_ac, [],[],[],[],LB,UB,array_calib_nonlcon_fh,options);
  
  exitflag = exitflag
elseif 0
  % Local solver (patternsearch solver) .. Too slow for TolMesh>1e-2
  options =  psoptimset('TolX',1e-6,'TolFun',1e-6,'TolMesh',1e-4,'InitialMeshSize',0.001,'MaxIter',10^8,'MaxFunEvals',10^8,...
    'PlotFcn',@psplotbestf);
  % options =  psoptimset('TolX',1e-15,'InitialMeshSize',0.001,'MaxIter',10^100,'MaxFunEvals',10^100,'TolMesh',1e-15);
  ac_est_error = patternsearch(@(est_errors) array_calibration_cost_2D(est_errors,ac_cost_params),...
    initial_ac, [],[],[],[],LB,UB,[],options);
elseif 0
  % Global solver
  % Sometimes it fails and recover .. Wast of time
  gs = GlobalSearch;
  objFun = @(est_errors) array_calibration_cost_2D(est_errors,ac_cost_params);
  options = optimoptions(@fmincon,'TolX',1e-6,'TolFun',1e-3,'MaxIter',10^6,'MaxFunEvals',10^6, ...
    'PlotFcn',{@optimplotfval,@optimplotstepsize});
  problem = createOptimProblem('fmincon','x0',initial_ac,'objective',objFun,'lb',LB,'ub',UB,'options',options);
  ac_est_error = run(gs,problem);
elseif 0
  % Alternating Projection (AP)-like implementation:
  % It doesn't give better results than fmincon, which is much faster.
  %
  %     options =  optimoptions(@fmincon,'TolX',1e-6,'TolFun',1e-6,'MaxIter',10^8,'MaxFunEvals',10^8, ...
  %       'PlotFcn',{@optimplotfval,@optimplotstepsize});%,'Algorithm','sqp');
  options =  optimoptions(@fmincon,'TolX',1e-6,'TolFun',1e-6,'MaxIter',10^8,'MaxFunEvals',10^8);%,'Algorithm','sqp');
  
  LB_mtx = repmat(error_bounds(:,1).',[Nc 1]); % Nc by 6 matrix
  UB_mtx = repmat(error_bounds(:,2).',[Nc 1]);
  ac_est_error_chan = reshape(initial_ac,[Nc 6]);
  ac_est_error_loop = zeros(Nc,6);
  exitflag_all = [];
  stop_metric = 0;
  Tol = 1e-4*ones(6,1);
  Ntrials = 1;
  while (all(stop_metric > Tol) || Ntrials <= 1e2)
    for chan_idx = 1:Nc
      % Update bounds: trun off the bounds on all channels escept the
      % current.
      LB_chan = LB_mtx;
      LB_chan([1:chan_idx-1 chan_idx+1:end],:) = 0;
      LB_chan = LB_chan(:);
      
      UB_chan = UB_mtx;
      UB_chan([1:chan_idx-1 chan_idx+1:end],:) = 0;
      UB_chan = UB_chan(:);
      
      [ac_est_error, ~, exitflag] = fmincon(@(est_errors) array_calibration_cost_2D(est_errors,ac_cost_params),...
        initial_ac, [],[],[],[],LB_chan,UB_chan,array_calib_nonlcon_fh,options);
      
      exitflag_all(end+1) = exitflag;
      
      % Update initialization
      ac_est_error_mtx = reshape(ac_est_error,[Nc 6]);
      ac_est_error_chan(chan_idx,:) = ac_est_error_mtx(chan_idx,:);
      initial_ac = ac_est_error_chan(:);
      sprintf('\n Working on channel %u ... Loop %u \n', chan_idx, Ntrials)
      %       pause(1)
    end
    
    stop_metric = mean(abs(ac_est_error_loop - ac_est_error_chan),1).';
    ac_est_error_loop = ac_est_error_chan;
    Ntrials = Ntrials+1;
  end
  ac_est_error = ac_est_error_loop(:);
end
toc

ac_est_error = reshape(ac_est_error,[Nc 6]);

% Measure the accuracy of the result
% ----------------------------------
% 1) RMSE over all sensors for each error type
sprintf('\n')
rmse           = sqrt(mean(abs(ac_est_error-Err).^2,1));
% 2) Mean of the absolute error over all sensors for each error type
sprintf('\n')
mean_abs_error = mean(abs(ac_est_error-Err),1);
sprintf('\n')
% 3) Mean abs error relative to the maximum error for each error type
relative_max   = mean_abs_error./max(abs(Err));
if isnan(relative_max(3))
  relative_max(3) = 0;
end

sprintf('\nRMSE = ')
disp(rmse)
sprintf('\nMean absolute error = ')
disp(mean_abs_error)
sprintf('\nRelative max = ')
disp(relative_max)

% -----------------------------------------------------------------------
%           Plot the gain and phase deviation patterns for each sensor
% -----------------------------------------------------------------------

if 1
  %  theta_RP = linspace(-90,90,2048)*pi/180;
  sv_params.array_gain_model_fh = param.array_gain_model_fh;
  sv_params.src.y_pc = y_pc;
  sv_params.src.z_pc = z_pc;
  sv_params.src.fc   = fc;
  % Array error parameters
  extra_error_params = [];
  extra_error_params.error_ypc      = error_ypc;
  extra_error_params.error_zpc      = error_zpc;
  extra_error_params.error_phase    = error_phase;
  extra_error_params.error_g_s      = error_g_s;
  extra_error_params.error_g_p      = error_g_p;
  extra_error_params.error_g_offset = error_g_offset;
  
  % Estimated array errors (i.e. calibration parameters)
  % You don't need these calib_params in real data case. Here we use them
  % to see how close the estimated errors to the actual errors. Ideally,
  % if we subtract the phase terms of the parameters from the actual
  % phase errors we would get 0 deg. Also, if we divide the actual gain
  % errors by the estimated gain errors we would get 1, which is the
  % ideal case.
  calib_params = [];
  calib_params.calib_ypc      = ac_est_error(:,1);% * lambda;
  calib_params.calib_zpc      = ac_est_error(:,2);% * lambda;
  calib_params.calib_phase    = ac_est_error(:,3);
  calib_params.calib_g_s      = ac_est_error(:,4);
  calib_params.calib_g_p      = ac_est_error(:,5);
  calib_params.calib_g_offset = ac_est_error(:,6);
  if 0
    % Debug
    calib_params.calib_ypc      = error_ypc;
    calib_params.calib_zpc      = error_zpc;
    calib_params.calib_phase    = error_phase;
    calib_params.calib_g_s      = error_g_s;
    calib_params.calib_g_p      = error_g_p;
    calib_params.calib_g_offset = error_g_offset;
  end
  %     sv_params.calib_params      = calib_params;
  
  DOA = 0*pi/180; % Plot radiation patterns for this DOA
  theta_RP = linspace(-90,90,2048)*pi/180;
  [~, DOA_idx] = min(abs(DOA-theta_RP));
  w = hanning(Nc); %ones(Nc,1);
  
  % Plot 3dB radiation pattern for 3 cases
  % ---------------------------------------------------------------------
  % Case 1: ideal radiation pattern
  sv_params.extra_error_params = [];
  sv_params.calib_params       = [];
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs = A(theta_RP,sv_params);
  A0 = SVs(:,DOA_idx);
  
  RP = abs((w.*A0)'*SVs).^2;
  RP = RP./max(RP);
  RP_dB = 10*log10(RP);
  
  idxs_3dB  = find(RP_dB>=-3);
  RP_3dB    = RP_dB(idxs_3dB);
  theta_3dB = theta_RP(idxs_3dB);
  
  beam_doa_lims = [min(theta_3dB) max(theta_3dB)];
  
  RP_dB_ideal = RP_dB;
  beam_doa_lims_vec(:,1) = beam_doa_lims;
  
  figure(100);clf
  subplot(131)
  plot(theta_RP*180/pi,RP_dB,'b')
  xlabel('\theta^\circ')
  ylabel('Power (dB)')
  title('No array errors')
  xlim([beam_doa_lims(1) beam_doa_lims(2)]*180/pi)
  grid on
  
  % Case 2: radiation pattern of a calibrated array (i.e. with errors)
  sv_params.extra_error_params = extra_error_params;
  sv_params.calib_params       = calib_params;
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs = A(theta_RP,sv_params);
  A0 = SVs(:,DOA_idx);
  
  RP = abs((w.*A0)'*SVs).^2;
  RP = RP./max(RP);
  RP_dB = 10*log10(RP);
  
  idxs_3dB  = find(RP_dB>=-3);
  RP_3dB    = RP_dB(idxs_3dB);
  theta_3dB = theta_RP(idxs_3dB);
  
  beam_doa_lims = [min(theta_3dB) max(theta_3dB)];
  
  RP_dB_calib = RP_dB;
  beam_doa_lims_vec(:,2) = beam_doa_lims;
  
  figure(100);
  subplot(132)
  plot(theta_RP*180/pi,RP_dB,'b')
  xlabel('\theta^\circ')
  %     ylabel('Power (dB)')
  title('Calibrated array')
  xlim([beam_doa_lims(1) beam_doa_lims(2)]*180/pi)
  grid on
  
  % Case 3: radiation pattern of an uncalibrated array
  sv_params.extra_error_params = extra_error_params;
  sv_params.calib_params       = [];
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs = A(theta_RP,sv_params);
  A0 = SVs(:,DOA_idx);
  
  RP = abs((w.*A0).'*SVs).^2;
  RP = RP./max(RP);
  RP_dB = 10*log10(RP);
  
  idxs_3dB  = find(RP_dB>=-3);
  RP_3dB    = RP_dB(idxs_3dB);
  theta_3dB = theta_RP(idxs_3dB);
  
  beam_doa_lims = [min(theta_3dB) max(theta_3dB)];
  
  RP_dB_no_calib = RP_dB;
  beam_doa_lims_vec(:,3) = beam_doa_lims;
  
  figure(100);
  subplot(133)
  plot(theta_RP*180/pi,RP_dB,'b')
  xlabel('\theta^\circ')
  %     ylabel('Power (dB)')
  title('Uncalibrated array')
  xlim([beam_doa_lims(1) beam_doa_lims(2)]*180/pi)
  grid on
  
  suptitle(sprintf('Array 3dB radiation pattern of %u deg. target',DOA*180/pi))
  
  % Plot the three radiation patterns above in one plot
  figure(101);clf
  hold on
  plot(theta_RP*180/pi,RP_dB_ideal,'b')
  plot(theta_RP*180/pi,RP_dB_calib,'r')
  plot(theta_RP*180/pi,RP_dB_no_calib,'k')
  xlim([beam_doa_lims_vec(1,1) beam_doa_lims_vec(2,1)]*180/pi)
  %     xlim([min(beam_doa_lims_vec(1,:)) max(beam_doa_lims_vec(2,:))]*180/pi)
  xlabel('\theta^\circ')
  ylabel('Power (dB)')
  title(sprintf('Array 3dB radiation pattern of %u deg. target',DOA*180/pi))
  grid on
  legend('No array errors','Calibrated array','Uncalibrated array','Location','best')
  
  % Plot sensors gain pattern before and after array calibration
  % ---------------------------------------------------------------------
  % Case 1: before array calibration (i.e. with array errors)
  sv_params.extra_error_params = extra_error_params;
  sv_params.calib_params       = [];
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs = A(theta_RP,sv_params);
  
  figure(102);clf
  subplot(121)
  hold on;
  for chan_idx = 1:Nc
    chan_resp = SVs(chan_idx,:);
    chan_gain = abs(chan_resp).^2;
    chan_gain = chan_gain./max(chan_gain);
    chan_gain_dB = 10*log10(chan_gain);
    plot(theta_RP*180/pi,chan_gain_dB)
  end
  xlabel('\theta^\circ')
  ylabel('Power (dB)')
  title('Uncalibrated array')
  xlim(DOA*180/pi+[-30 +30])
  grid on
  legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
  
  % Case 2: after array calibration (i.e. with array errors)
  sv_params.extra_error_params = extra_error_params;
  sv_params.calib_params       = calib_params;
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs = A(theta_RP,sv_params);
  
  %     figure(103);clf
  subplot(122)
  hold on;
  for chan_idx = 1:Nc
    chan_resp = SVs(chan_idx,:);
    chan_gain = abs(chan_resp).^2;
    chan_gain = chan_gain./max(chan_gain);
    chan_gain_dB = 10*log10(chan_gain);
    plot(theta_RP*180/pi,chan_gain_dB)
  end
  xlabel('\theta^\circ')
  %     ylabel('Power (dB)')
  title('Calibrated array')
  xlim(DOA*180/pi+[-30 +30])
  grid on
  %     legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
  suptitle('Sensors gain pattern')
  
  % Plot sensors phase pattern before and after array calibration
  % ---------------------------------------------------------------------
  % Case 1: before array calibration (i.e. with array errors)
  sv_params.extra_error_params = extra_error_params;
  sv_params.calib_params       = [];
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs = A(theta_RP,sv_params);
  
  figure(103);clf
  subplot(211)
  hold on
  for chan_idx = 1:Nc
    chan_resp = SVs(chan_idx,:);
    chan_phase = angle(chan_resp)*180/pi;
    plot(theta_RP*180/pi,chan_phase)
  end
  %     xlabel('\theta^\circ')
  ylabel('Phase (deg.)')
  title('Uncalibrated array')
  xlim(DOA*180/pi+[-10 +10])
  grid on
  legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
  
  % Case 2: after array calibration (i.e. with array errors)
  sv_params.extra_error_params = extra_error_params;
  sv_params.calib_params       = calib_params;
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs = A(theta_RP,sv_params);
  
  figure(103);
  subplot(212)
  hold on
  for chan_idx = 1:Nc
    chan_resp = SVs(chan_idx,:);
    chan_phase = angle(chan_resp)*180/pi;
    plot(theta_RP*180/pi,chan_phase)
  end
  xlabel('\theta^\circ')
  ylabel('Phase (deg.)')
  title('Calibrated array')
  xlim(DOA*180/pi+[-10 +10])
  grid on
  %     legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
  
  suptitle('Sensors phase pattern')
  
  % Plot sensors phase deviation before and after array calibration
  % ---------------------------------------------------------------------
  % Case 1: before array calibration (i.e. with array errors)
  sv_params.extra_error_params = [];
  sv_params.calib_params       = [];
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs_ideal = A(theta_RP,sv_params);
  
  sv_params.extra_error_params = extra_error_params;
  sv_params.calib_params       = [];
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs = A(theta_RP,sv_params);
  
  SV_phase_dev = conj(SVs_ideal) .* SVs; % (ideal+error)-ideal
  
  figure(104);clf
  subplot(211)
  hold on
  for chan_idx = 1:Nc
    chan_resp = SV_phase_dev(chan_idx,:);
    chan_phase = angle(chan_resp)*180/pi;
    plot(theta_RP*180/pi,chan_phase)
  end
  %     xlabel('\theta^\circ')
  ylabel('Phase (deg.)')
  title('Uncalibrated array')
  xlim(DOA*180/pi+[-10 +10])
  grid on
  legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
  
  % Case 2: after array calibration (i.e. with array errors)
  sv_params.extra_error_params = [];
  sv_params.calib_params       = [];
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs_ideal = A(theta_RP,sv_params);
  
  sv_params.extra_error_params = extra_error_params;
  sv_params.calib_params       = calib_params;
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs = A(theta_RP,sv_params);
  
  SV_phase_dev = conj(SVs_ideal) .* SVs; % (ideal+error-calib)-ideal
  
  figure(104);
  subplot(212)
  hold on
  for chan_idx = 1:Nc
    chan_resp = SV_phase_dev(chan_idx,:);
    chan_phase = angle(chan_resp)*180/pi;
    plot(theta_RP*180/pi,chan_phase)
  end
  xlabel('\theta^\circ')
  ylabel('Phase (deg.)')
  title('Calibrated array')
  xlim(DOA*180/pi+[-10 +10])
  grid on
  %     legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
  
  suptitle('Sensors phase deviation pattern')
  
end

if param.debug_level >= 3
  return
end
