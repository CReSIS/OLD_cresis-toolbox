%%    Array Calibration: The new implementation
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
  fc = 195e6;
  BW = 30e6;
  fs = BW;
  lambda = c/fc;
%   targets_doa = [-85:5:85].'*pi/180;
%   targets_doa = [-1.2154, -1.0654, -0.9484, -0.8481, -0.7580, -0.6751, -0.5974, -0.5236, -0.4528, -0.3844, -0.3178, -0.2527, ...
%     -0.1886, -0.1253, -0.0625, 0, 0.0625, 0.1253, 0.1886, 0.2527, 0.3178, 0.3844, 0.4528, 0.5236, 0.5974, 0.6751, ...
%     0.7580, 0.8481, 0.9484, 1.0654, 1.2154]';
% targets_doa = [-1.0291 -0.9742 -0.9075 -0.8366 -0.7684 -0.6867 -0.6303 -0.5606 -0.4861 -0.4190 -0.3480 -0.2773 -0.2127 -0.1398 ...
%     -0.0679 0.0050 0.0671 0.1381 0.2086 0.2804 0.3497 0.4198 0.4894 0.5502 0.6294 0.7036 0.7491 0.8367 0.9089 0.9836 1.0389 1.0913].';
% targets_doa = [-0.4407   -0.3534   -0.2651   -0.1840   -0.0995    0.0034    0.0581    0.1609    0.2494    0.3394    0.4244].';
% targets_doa(abs(targets_doa)>40*pi/180) = [];

delta_theta = 0.75; % In degrees
roll_bins_centers = [-42:2*delta_theta:42].'* pi/180;
roll_bins_centers(roll_bins_centers==0) = [];
roll_angles = roll_bins_centers(roll_bins_centers>-pi/2 & roll_bins_centers<pi/2);
roll_angles_deg = roll_angles.'*180/pi;

flight_height  = 1500;
  range_vec      = flight_height./cos(roll_angles);
  max_range      = max(range_vec);
  min_range      = min(range_vec);
  
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.t0                      = 2*(min_range-500)/c;
  param.src.t1                      = 2*(max_range+500)/c;
  param.src.ft_func                 = @(t) tukeywin_cont(t * BW);
  param.src.ft_wind                 = @(N) hanning(N);
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  % Nc: number of sensors
  Nc = 7;
  % DOA method parameters
  param.method.list                   = [7];
  %   param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
  param.src.noise_power             = zeros(1,Nc);
  
%   % Nsep: number of lambda/4 steps between sensors
%   Nsep = 1;
%   % Arguments for a linear array in y-dimension
%   % param.src.lever_arm.fh_args: Arguments for a linear array in y-dimension
%   param.src.lever_arm.fh_args       = {[], 1, 1:Nc, [0; Nsep*c/fc/2; 0]};
%   [phase_center] = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
%   y_pc = phase_center(2,:).';
%   y_pc = y_pc-y_pc(1);
%   z_pc = 0*lambda  + phase_center(3,:).'; % Make it a multiple of lambda/2 for best results
%   z_pc = z_pc-z_pc(1);
%   
%   phase_center(2,:) = y_pc;
%   phase_center(3,:) = z_pc;
%   param.src.phase_center = -phase_center;
  
 % Nsep: number of lambda/4 steps between sensors
  Nsep = 1;
  param.src.lever_arm.fh_args       = {[], 1, 1:Nc, [0; Nsep*c/fc/2; 0]};
  phase_center = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
  
  % Use this if you want a non-level array (i.e. there z component).
%   phase_center(3,:) = phase_center(3,:) + cumsum(0.5*lambda*ones(Nc,1)).';
%     phase_center(3,:) = phase_center(3,:) + cumsum(0.5*lambda*ones(Nc,1)).';

  % Change reference sensor if you want to
  ref_chan = 1;
  phase_center(2,:) = phase_center(2,:) - phase_center(2,ref_chan);
  phase_center(3,:) = phase_center(3,:) - phase_center(3,ref_chan);
  
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
  surf_param.rcs_in.var = 1e2;
  param.snr_db = 20;%-10*log10(Nc);
  surf_param.dy = 10;
  %   surf_param.y_range = [-2500 2500];
  surf_param.dx = 10;
  surf_param.x_range = [-2500 2500];
  %   surf_param.x = [-500:1:500]; % Defined later
  %   surf_param.y = [-1500:20:1500].';
  surf_param_y = range_vec .* sin(roll_angles);
%   surf_param.y = 0;
  % y_range should >= maximum y (i.e. targets shoulb be inside the imaged
  % swath)
  surf_param.y_range = [-1.5*max(surf_param_y) 1.5*max(surf_param_y)];
  if max(surf_param.y_range) == 0
    surf_param.y_range = [-2500 2500];
  end
%     param.monte.target_param{1} = surf_param;
  % -----------------------------------------------------------------------
  %%                   Array Processing parameters
  % -----------------------------------------------------------------------
  array_param = [];
  
  % NN: total length of the sensor array (in lambda/4 units)
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
  array_param.rline_rng = -10:10;
  
  % Nsig: maximum number of targets (not the actual number of targets)
  Nsig = 1;
  array_param.Nsig = Nsig;
  
  array_param.init = 'ap';
  array_param.theta_guard = 1.5;%(max(theta)-min(theta))/(4*Nc);
  
  array_param.W = 1;
  dt = 1/fs;
  array_param.imp_resp.time_vec = -3*array_param.W*dt : dt/8 : 3*array_param.W*dt;
  BW = abs(param.src.f1 - param.src.f0);
  array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / BW);
  
  for idx = 1:array_param.Nsig
    array_param.doa_constraints(idx).method = 'fixed';
    array_param.doa_constraints(idx).init_src_limits = [-90 90];%[min(theta) max(theta)]*180/pi;
    array_param.doa_constraints(idx).src_limits = [-90 90];%[min(theta) max(theta)]*180/pi;
  end
  
  param.array_param = array_param;
  
  N_skipped_rlines = 1;
  N_reqd_rlines    = 1;
  surf_param.x = [-750:N_skipped_rlines:-750+N_skipped_rlines*(N_reqd_rlines+2*max(array_param.rline_rng))-1];
%   param.monte.target_param{1} = surf_param;
  
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
%   error_bounds =  [-0.2 0.2; -0.2 0.2; -1*pi 1*pi; -10 10; -1*pi 1*pi; -10 10];
%    error_bounds =  [-0.2 0.2; -0.2 0.2; -1*pi 1*pi; 0 0; 0 0; 0 0];
error_bounds = [-inf  inf; -inf  inf; 0  0 ...
               ;0  0; 0  0; 0  0];
  
% error_ypc      = zeros(Nc,1);
% error_zpc      = zeros(Nc,1);
% error_phase    = zeros(Nc,1);
% error_g_s      = zeros(Nc,1);
% error_g_p      = zeros(Nc,1);
% error_g_offset = zeros(Nc,1);

% y_pc and z_pc errors are generated as percentage fo the actual y_pc
% and z_pc values. This guarantees that the errors are sensible.
% y_pc_err_percentage = [1  20 -10   4   -9   2   8].' ./100;
z_pc_err_percentage = [1  -7  1   -5    20  5   10].' ./100;
% y_pc_err_percentage = [1 20 10].' ./100;
% z_pc_err_percentage = [1 7 10].' ./100;

error_ypc      = phase_center(2,:).' .* y_pc_err_percentage;%  In meters units
error_zpc      = phase_center(2,:).' .* z_pc_err_percentage;%  In meters units
% error_zpc      = phase_center(3,:).' .* z_pc_err_percentage;%  In meters units
error_phase    = zeros(Nc,1);%
error_g_s      = zeros(Nc,1);%[0 0.8 2 0.9 1 0.8 5]';
error_g_p      = zeros(Nc,1);
error_g_offset = zeros(Nc,1);

clear Err
Err(:,1) = error_ypc;%./lambda;
Err(:,2) = error_zpc;%./lambda;
Err(:,3) = error_phase;
Err(:,4) = error_g_s;
Err(:,5) = error_g_p;
Err(:,6) = error_g_offset;

% Array calibration errors
% error_params.error_ypc      = error_ypc ;
% error_params.error_zpc      = error_zpc ;
error_params.error_phase    = error_phase;
error_params.error_g_s      = error_g_s;
error_params.error_g_p      = error_g_p;
error_params.error_g_offset = error_g_offset;

% error_params_tmp = error_params;
param.error_params = error_params; % With array errors
    
  % -----------------------------------------------------------------------
  %                       Run the simulation
  % -----------------------------------------------------------------------  
  if 1
    est_doa_before_calib = [];
    Rxx_all = [];
    est_doa = [];
    clear ypc_all zpc_all
    for target_i = 1:length(roll_angles)
      roll_angle_tmp = roll_angles(target_i);
      
%       surf_param.y = surf_param_y(target_i);
      surf_param.y = 0;
      param.monte.target_param{1} = surf_param;
     
      % Account for the roll angle, if any. Default is theta_roll=0. This
      % can be done by rotating the level array (i.e. no roll) by theta_roll.
%       theta_roll = 0*pi/180;
      theta_roll = roll_angle_tmp;
      if 0
        % Wrong way .. Remove later
      for chan_idx = 1:Nc
        d_chan_ref = sqrt(abs(diff(phase_center(2,[ref_chan  chan_idx])))^2 + abs(diff(phase_center(3,[ref_chan  chan_idx])))^2);
        ypc_all(chan_idx,target_i) = d_chan_ref .* cos(theta_roll);
        zpc_all(chan_idx,target_i) = d_chan_ref .* sin(theta_roll);
%         d_chan_ref_error = sqrt(abs(diff(error_ypc([ref_chan  chan_idx])))^2 + abs(diff(error_zpc([ref_chan  chan_idx])))^2);
%         ypc_error(chan_idx) = d_chan_ref_error .* cos(theta_roll);
%         zpc_error(chan_idx) = d_chan_ref_error .* sin(theta_roll);
      end
      elseif 1
        % This is the right way of rotating the array
        % rotation_mtx: rotates any point (sensor) counter-clockwise by 
        % theta_roll. To rotate back, change the sign of theta_roll. 
      rotation_mtx = [cos(theta_roll), -sin(theta_roll);sin(theta_roll), cos(theta_roll)];
      rotated_phase_center = rotation_mtx * phase_center([2 3],:) *(-1);
      
       % Change reference sensor to ref_chan
%       rotated_phase_center(1,:) = rotated_phase_center(1,:) - rotated_phase_center(1,ref_chan);
%       rotated_phase_center(2,:) = rotated_phase_center(2,:) - rotated_phase_center(2,ref_chan);
  
      ypc_all(:,target_i) = rotated_phase_center(1,:).';
      zpc_all(:,target_i) = rotated_phase_center(2,:).';
      
      rotated_error = rotation_mtx*[error_ypc error_zpc].';
      ypc_error = rotated_error(1,:).';
      zpc_error = rotated_error(2,:).';
%       ypc_error = (rotated_error(1,:) - rotated_error(1,ref_chan)).';
%       zpc_error = (rotated_error(2,:) - rotated_error(2,ref_chan)).';
      end
      ypc_error_all(:,target_i) = ypc_error;
      zpc_error_all(:,target_i) = zpc_error;
      
      param.src.phase_center(2,:) = ypc_all(:,target_i);
      param.src.phase_center(3,:) = zpc_all(:,target_i);
      
      error_params.error_ypc      = ypc_error; 
      error_params.error_zpc      = zpc_error; 
      param.error_params          = error_params;
      
      results = crosstrack(param);
      
      actual_doa = NaN(length(results.actual_doa),1);
      for actual_doa_i = 1:length(results.actual_doa)
        if ~isnan(results.actual_doa{actual_doa_i})
        actual_doa(actual_doa_i) = results.actual_doa{actual_doa_i};
        rbin_i = actual_doa_i;
        end
      end
      
      est_doa(:,target_i) = results.tomo.doa(rbin_i,:);
      sim_data = squeeze(results.sim_data{1}(rbin_i,:,:,:,:)).'; % Nt*Nx*Nc data matrix
      Data{target_i} = sim_data; % Used for calibration plots
      Rxx_all{end+1} = 1/size(sim_data,2) * sim_data * sim_data';
            
%       good_doa_i = find(~isnan(results.tomo.doa(:)));
%       if ~isempty(good_doa_i)
%         est_doa_before_calib(end+1) = results.tomo.doa(good_doa_i);
%         sim_data = squeeze(results.sim_data{1}(good_doa_i,:,:,:,:)).'; % Nt*Nx*Nc data matrix
%         Data{target_i} = sim_data;
%         Rxx_all{end+1} = 1/size(sim_data,2) * sim_data * sim_data';
%         %       actual_doa(end+1) = targets_doa(target_i);
%         if ~isempty(results.actual_doa{good_doa_i})
%           est_doa(end+1) = results.actual_doa{good_doa_i};
%         else
%           est_doa(end+1) = 0;%targets_doa(target_i);
%         end
%       end
    end
%     est_doa_deg = est_doa*180/pi;
  end  
  
  % The ideal DCM
  if 1
    Rxx_ideal_all = [];
    param.error_params.error_ypc      = zeros(Nc,1);
    param.error_params.error_zpc      = zeros(Nc,1);
    param.error_params.error_phase    = zeros(Nc,1);
    param.error_params.error_g_s      = zeros(Nc,1);
    param.error_params.error_g_p      = zeros(Nc,1);
    param.error_params.error_g_offset = zeros(Nc,1);
    for target_i = 1:length(roll_angles)
      roll_angle_tmp = roll_angles(target_i);
      
%       surf_param.y = surf_param_y(target_i);
      surf_param.y = 0;
      param.monte.target_param{1} = surf_param;
      
      param.src.phase_center(2,:) = ypc_all(:,target_i);
      param.src.phase_center(3,:) = zpc_all(:,target_i);
      
      results = crosstrack(param);
      
      actual_doa = NaN(length(results.actual_doa),1);
      for actual_doa_i = 1:length(results.actual_doa)
        if ~isnan(results.actual_doa{actual_doa_i})
        actual_doa(actual_doa_i) = results.actual_doa{actual_doa_i};
        rbin_i = actual_doa_i;
        end
      end
      
%       rbin_i = find(~isnan(results.tomo.doa(:)));
      sim_data = squeeze(results.sim_data{1}(rbin_i,:,:,:,:)).'; % Nt*Nx*Nc data matrix
      Rxx_ideal_all{end+1} = 1/size(sim_data,2) * sim_data * sim_data';
        
%        good_doa_i = find(~isnan(results.tomo.doa(:)));
%       if ~isempty(good_doa_i)
%         sim_data = squeeze(results.sim_data{1}(good_doa_i,:,:,:,:)).'; % Nt*Nx*Nc data matrix
% %         Data{target_i} = sim_data;
%         Rxx_ideal_all{end+1} = 1/size(sim_data,2) * sim_data * sim_data';
%       end
    end
   
  end
   
% param.error_params  = error_params_tmp; 

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

  ac_cost_params.roll_angles = roll_angles;
  ac_cost_params.y_pc       = ypc_all;% + ypc_error_all;
  ac_cost_params.z_pc       = zpc_all;% + zpc_error_all; 
  ac_cost_params.fc         = fc;
  ac_cost_params.BW = BW;
  ac_cost_params.Rxx_all = Rxx_all;
%   ac_cost_params.sim_data   = sim_data;
  ac_cost_params.array_param = array_param;
  ac_cost_params.ref_chan = ref_chan;
%   ac_cost_params.array_gain_model_fh = param.array_gain_model_fh;
  
  tic
  % -----------------------------------------------------------------------
  %                          Run the optimizer
  % -----------------------------------------------------------------------
  LB = repmat(error_bounds(:,1).',[Nc 1]);
  UB = repmat(error_bounds(:,2).',[Nc 1]);
  
%   ref_chan = 1;
  LB(ref_chan,:) = 0;
  UB(ref_chan,:) = 0;
%   LB(ceil(Nc/2),:) = 0;
%   UB(ceil(Nc/2),:) = 0;
  
  LB = LB(:);
  UB = UB(:);
  initial_ac = 0*ones(size(LB)); % Nc*6 matrix
%   initial_ac = -2*rand(size(LB))+1;
  
  [nadir_doa,nadir_doa_idx] = min(abs(roll_angles));
  
%   param.fc = fc;
%   param.nadir_y_pc = y_pc;
%   param.nadir_z_pc = z_pc;
%   param.nadir_doa  = nadir_doa;
%   param.doa = targets_doa.';
  array_calib_nonlcon_fh = []; %@(est_errors) array_calib_nonlcon(est_errors,param);
  
  if 1
    % Local solver (fmincon solver)
    % Global solver then local steeoest descent
%     options =  psoptimset('TolX',1e-3,'TolFun',1e-3,'TolMesh',1e-3,'InitialMeshSize',0.001,'MaxIter',10^6,'MaxFunEvals',10^6,...
%       'PlotFcn',@psplotbestf);
%     % options =  psoptimset('TolX',1e-15,'InitialMeshSize',0.001,'MaxIter',10^100,'MaxFunEvals',10^100,'TolMesh',1e-15);
%     [ac_est_error, fval, exitflag] = patternsearch(@(est_errors) array_calibration_cost_2D(est_errors,ac_cost_params),...
%       initial_ac, [],[],[],[],LB,UB,array_calib_nonlcon_fh,options);    
%     
%     initial_ac = ac_est_error(:);
    options =  optimoptions(@fmincon,'TolX',1e-6,'TolFun',1e-6,'MaxIter',10^7,'MaxFunEvals',10^7, ...
      'PlotFcn',{@optimplotfval,@optimplotstepsize});%,'Algorithm','sqp');
    [ac_est_error, fval, exitflag] = fmincon(@(est_errors) array_calibration_cost_2D(est_errors,ac_cost_params),...
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
  
  ac_est_error = reshape(ac_est_error,[Nc length(ac_est_error)/Nc]);
  
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
  relative_max(isnan(relative_max)) = 0;
  % 4) Absolute error
  abs_error = abs(ac_est_error-Err);
  
  sprintf('\nRMSE = ')
  disp(rmse)
  sprintf('\nMean absolute error = ')
  disp(mean_abs_error)
  sprintf('\nRelative max = ')
  disp(relative_max)
  sprintf('\nAbsolute error = ')
  disp(abs_error)
  
  return
  % -----------------------------------------------------------------------
  %%           Plot the gain and phase deviation patterns for each sensor
  % -----------------------------------------------------------------------
  
  if 1
    DOA = 12* pi/180;
    N = 1;
    [~,doa_idx] = min(abs(DOA-actual_doa));
    
    phase_unwrapping = 1;
    k = 4*pi/(c/fc);
    calib_params = [];
    calib_params.calib_ypc      = ac_est_error(:,1) * sin(actual_doa(doa_idx));
    calib_params.calib_zpc      = ac_est_error(:,2) * cos(actual_doa(doa_idx));
    calib_params.calib_phase    = ac_est_error(:,3);% + k*ac_est_error(:,2);
    calib_params.calib_g_s      = ac_est_error(:,4);
    calib_params.calib_g_p      = ac_est_error(:,5);
    calib_params.calib_g_offset = ac_est_error(:,6);
    
    % actual array phase centers and length
    uncalib_ypc = ac_cost_params.y_pc(:,doa_idx);
    uncalib_zpc = ac_cost_params.z_pc(:,doa_idx);

    calib_ypc = uncalib_ypc + calib_params.calib_ypc ;
    calib_zpc = uncalib_zpc + calib_params.calib_zpc - ac_est_error(:,2);
    
    if 0
      % Debug: calculate the angle between y and z coordinates from geometry.
      % This angle should match the roll angle (or DOA defined above).
      DOA_hat_uncalib = nanmean(atan(uncalib_zpc./uncalib_ypc));
      DOA_hat_calib = nanmean(atan(calib_zpc./calib_ypc));
      
      h1 = sprintf('\nActual roll angle is %2.2f\n',actual_doa(doa_idx)*180/pi);
      h2 = sprintf('Uncalib. array (average) roll angle from geometry is %2.2f deg. \n',DOA_hat_uncalib*180/pi);
      h3 = sprintf('Calib. array (average) roll angle from geometry is %2.2f deg. \n',DOA_hat_calib*180/pi);
      disp(h1)
      disp(h2)
      disp(h3)
    end
    
    %% Plot phase centers of the sensors
    if 1
      figure(9998);clf
      hold on
      % Actual phase centers
      plot(uncalib_ypc,uncalib_zpc ,'b*','LineWidth',2)
      % Calibrated phase centers
      plot(calib_ypc,calib_zpc ,'r*','LineWidth',1.5)
      % Calibrated phase centers projected onto y axis
      %     plot(y_pos{fn_idx}{doa_idx}+ac_est_error(:,1)*sin(actual_doa(doa_idx)),z_pos{fn_idx}{doa_idx},'k*')
      
      xlim([min(min(uncalib_ypc),min(calib_ypc))  max(max(uncalib_ypc),max(calib_ypc))])
      ylim([min(min(uncalib_zpc),min(calib_zpc))  max(max(uncalib_zpc),max(calib_zpc))])
      
      grid on
      xlabel('y-axis')
      ylabel('z-axis')
      title('Phase centers, PCs, of the sensors')
      legend('Uncalib. PCs','Calib. PCs','Location','southeast')
      %   legend('Actual PCs','Calib. PCs','calb. PCs projected onto y-axis','Location','best')
    end
    %% Plot DCM magnitude and angle
    if 1
      % 1- Before calibration
      % --------------------------
      if 1
        % L0:array length, L_doa:array length projected in the direction of
        % DOA, L_extra:extra distance the signal travels from sensor 1 to
        % sensor Nc (from which we determine the maximum phase difference
        % accross the array for a target at angle DOA). The DCM should show the
        % same phase between the first and last sensors (Sanity check).
        % Remember that the phase centers for a roll=DOA are projected onto the
        % y-axis (the nominal array axis) and now the array is no longer orthogomal
        % to the range vector in the direction of DOA. So, the new array length
        % is shorter and the maximum phase accross the array from a target at DOA
        % is the two way phase (i.e. 4*pi/lambda * L_extra).
        
        % 1) Maximum phase of the array before calibration
        %   L0 = lambda/4*(Nc-1);
        L0 = sqrt(abs(diff(uncalib_ypc([1 end])))^2 + abs(diff(uncalib_zpc([1 end])))^2);
        L_doa  = L0*cos(actual_doa(doa_idx));
        L_extra = L0*sin(actual_doa(doa_idx));
        max_actual_expected_phase = 4*pi/lambda * L_extra * 180/pi;
        h1 = sprintf('\nMaximum phase before calib. is %4.2f deg. (%4.2f after wraping to pi) \n',max_actual_expected_phase,wrapToPi(max_actual_expected_phase*pi/180)*180/pi);
        
        % 2) Maximum phase of the array after calibration
        L0 = sqrt(abs(diff(calib_ypc([1 end])))^2 + abs(diff(calib_zpc([1 end])))^2);
        L_doa  = L0*cos(actual_doa(doa_idx));
        L_extra = L0*sin(actual_doa(doa_idx));
        max_new_expected_phase = 4*pi/lambda * L_extra * 180/pi;
        h2 = sprintf('\nMaximum phase after calib. is %4.2f deg. (%4.2f after wraping to pi) \n',max_new_expected_phase,wrapToPi(max_new_expected_phase*pi/180)*180/pi);
        
        disp(h1)
        disp(h2)
      end
      
      DataSample = squeeze(Data{doa_idx});      
      Rxx_uncalib = Rxx_all{doa_idx};
      
      figure(9990);clf
      Title = sprintf('DOA = %2.2f deg.',actual_doa(doa_idx)*180/pi);      
      t = suptitle(strcat('Before calib.: ',Title));
      
      subplot(211)
      abs_Rxx_uncalib = 10*log10(abs(Rxx_uncalib)./max(abs(Rxx_uncalib(:))));
      imagesc(abs_Rxx_uncalib)
      h = colorbar;
      ylabel(h,'Normalized |R| (dB)')
      %   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      
      subplot(212)
      phase_Rxx_uncalib = angle(Rxx_uncalib);
      if phase_unwrapping
        phase_Rxx_uncalib = unwrap(phase_Rxx_uncalib,[],1);
        phase_Rxx_uncalib = phase_Rxx_uncalib - diag(diag(phase_Rxx_uncalib));
        
        phase_Rxx_uncalib = unwrap(phase_Rxx_uncalib,[],2);
        phase_Rxx_uncalib = phase_Rxx_uncalib - diag(diag(phase_Rxx_uncalib));
      end
      imagesc(phase_Rxx_uncalib*180/pi)
      h = colorbar;
      ylabel(h,'\angle{R} (\circ)')
      h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      if ~isempty(Ticks)
        h.Ticks = Ticks;
      end
      
      % 2- After calibration
      % -------------------------
      gain_dB = calib_params.calib_g_s.*(sin(actual_doa(doa_idx)) - sin(calib_params.calib_g_p)).^2 + calib_params.calib_g_offset;
      phase_deg = (calib_params.calib_ypc*k*sin(actual_doa(doa_idx)) - calib_params.calib_zpc*k*cos(actual_doa(doa_idx)) + ...
        calib_params.calib_phase + + ac_est_error(:,2)*k)*180/pi;
      
      % Correct the DCM using the estimated calibration parameters.
      % The following two options are mathematically identical
      if 1
        calib_vec = 10.^(-gain_dB./20) .*exp(-1i*phase_deg*pi/180);
        DataSample_calib = repmat(calib_vec,[1 size(DataSample,2)]) .* DataSample;
        Rxx_calib =  1/size(DataSample_calib,2) * DataSample_calib*DataSample_calib';
      else
        calib_vec = 10.^(gain_dB./20) .*exp(1i*phase_deg*pi/180);
        H = diag(calib_vec);
        DataSample_calib = inv(H) * DataSample; % flipud(DataSample);
        Rxx_calib =  1/size(DataSample_calib,2) * DataSample_calib*DataSample_calib';
        % OR, do this
        %     Rxx_calib = inv(H) * Rxx_uncalib *inv(H');
      end
      
      figure(9991);clf
      suptitle(strcat('After calib.: ',Title))
      
      subplot(211)
      abs_Rxx_calib = 10*log10(abs(Rxx_calib)./max(abs(Rxx_calib(:))));
      imagesc(abs_Rxx_calib)
      h = colorbar;
      ylabel(h,'Normalized |R| (dB)')
      %   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      
      subplot(212)
      phase_Rxx_calib = angle(Rxx_calib);
      if phase_unwrapping
        phase_Rxx_calib = unwrap(phase_Rxx_calib,[],1);
        phase_Rxx_calib = phase_Rxx_calib - diag(diag(phase_Rxx_calib));
        
        phase_Rxx_calib = unwrap(phase_Rxx_calib,[],2);
        phase_Rxx_calib = phase_Rxx_calib - diag(diag(phase_Rxx_calib));
      end
      imagesc(phase_Rxx_calib*180/pi)
      h = colorbar;
      ylabel(h,'\angle{R} (\circ)')
      Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      if ~isempty(Ticks)
        h.Ticks = Ticks;
      end
      
      % 3) Before calibration with forward-backward averaging
      % ------------------------------------------------------
      if 1
        % FB averaging
        reflect_mtx = flipud(eye(size(Rxx_uncalib)));
        Rxx_uncalib_fb = (1/2)*(Rxx_uncalib + reflect_mtx * conj(Rxx_uncalib) * reflect_mtx);
      elseif 0
        % FB averaging AND spatial smoothing
        Rxx1 = 1/(size(DataSample(1:end-1,:),2))*DataSample(1:end-1,:) * DataSample(1:end-1,:)';
        Rxx2 = 1/(size(DataSample(2:end,:),2))*DataSample(2:end,:) * DataSample(2:end,:)';
        
        % Apply FB averaging for each subarray
        reflect_mtx = flipud(eye(size(Rxx1)));
        Rxx1_fb = (1/2)*(Rxx1 + reflect_mtx * conj(Rxx1) * reflect_mtx);
        Rxx2_fb = (1/2)*(Rxx2 + reflect_mtx * conj(Rxx2) * reflect_mtx);
        
        % Average the two DCMs
        Rxx_ss = 1/2*(Rxx1_fb+Rxx2_fb); % (Nc-1)-by-(Nc-1)
        
        % Handle the lost sensor (or DOF) such that the final DCM is Nc-by-Nc
        Rxx_tmp = zeros(Nc);
        Rxx_tmp(1:end-1,1:end-1) = Rxx_ss;
        Rxx_tmp(end,:) = Rxx_uncalib(end,:);
        Rxx_tmp(:,end) = Rxx_uncalib(:,end);
        
        Rxx_uncalib_fb = Rxx_tmp; % Nc-by-Nc matrix
      end
      
      figure(9992);clf
      suptitle(strcat('Before calib. and FB avging:',',', Title))
      %      suptitle(strcat('After array calib. and FB avging: ',frame_name(1:8),'\_',frame_name(10:11),'\_',frame_name(13:end),',',Title))
      
      subplot(211)
      imagesc(10*log10(abs(Rxx_uncalib_fb)./max(abs(Rxx_uncalib_fb(:)))))
      h = colorbar;
      ylabel(h,'Normalized |R| (dB)')
      %   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      
      subplot(212)
      phase_Rxx_uncalib_fb = angle(Rxx_uncalib_fb);
      if phase_unwrapping
        phase_Rxx_uncalib_fb = unwrap(phase_Rxx_uncalib_fb,[],1);
        phase_Rxx_uncalib_fb = phase_Rxx_uncalib_fb - diag(diag(phase_Rxx_uncalib_fb));
        
        phase_Rxx_uncalib_fb = unwrap(phase_Rxx_uncalib_fb,[],2);
        phase_Rxx_uncalib_fb = phase_Rxx_uncalib_fb - diag(diag(phase_Rxx_uncalib_fb));
      end
      imagesc(phase_Rxx_uncalib_fb*180/pi)
      h = colorbar;
      ylabel(h,'\angle{R} (\circ)')
      Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      if ~isempty(Ticks)
        h.Ticks = Ticks;
      end
      
      % 4) After calibration with forward-backward averaging
      % ----------------------------------------------------
      if 1
        % FB averaging
        reflect_mtx = flipud(eye(size(Rxx_calib)));
        Rxx_calib_fb = (1/2)*(Rxx_calib + reflect_mtx * conj(Rxx_calib) * reflect_mtx);
      elseif 0
        % FB averaging AND spatial smoothing
        Rxx1 = 1/(size(DataSample_calib(1:end-1,:),2))*DataSample_calib(1:end-1,:) * DataSample_calib(1:end-1,:)';
        Rxx2 = 1/(size(DataSample_calib(2:end,:),2))*DataSample_calib(2:end,:) * DataSample_calib(2:end,:)';
        
        % Apply FB averaging for each subarray
        reflect_mtx = flipud(eye(size(Rxx1)));
        Rxx1_fb = (1/2)*(Rxx1 + reflect_mtx * conj(Rxx1) * reflect_mtx);
        Rxx2_fb = (1/2)*(Rxx2 + reflect_mtx * conj(Rxx2) * reflect_mtx);
        
        % Average the two DCMs
        Rxx_ss = 1/2*(Rxx1_fb+Rxx2_fb); % (Nc-1)-by-(Nc-1)
        
        % Handle the lost sensor (or DOF) such that the final DCM is Nc-by-Nc
        Rxx_tmp = zeros(Nc);
        Rxx_tmp(1:end-1,1:end-1) = Rxx_ss;
        Rxx_tmp(end,:) = Rxx_calib(end,:);
        Rxx_tmp(:,end) = Rxx_calib(:,end);
        
        Rxx_calib_fb = Rxx_tmp; % Nc-by-Nc matrix
      end
      
      figure(9993);clf
      suptitle(strcat('After calib. and FB avging:',',', Title))
      %      suptitle(strcat('After array calib. and FB avging: ',frame_name(1:8),'\_',frame_name(10:11),'\_',frame_name(13:end),',',Title))
      
      subplot(211)
      imagesc(10*log10(abs(Rxx_calib_fb)./max(abs(Rxx_calib_fb(:)))))
      h = colorbar;
      ylabel(h,'Normalized |R| (dB)')
      %   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      
      subplot(212)
      phase_Rxx_calib_fb = angle(Rxx_calib_fb);
      if phase_unwrapping
        phase_Rxx_calib_fb = unwrap(phase_Rxx_calib_fb,[],1);
        phase_Rxx_calib_fb = phase_Rxx_calib_fb - diag(diag(phase_Rxx_calib_fb));
        
        phase_Rxx_calib_fb = unwrap(phase_Rxx_calib_fb,[],2);
        phase_Rxx_calib_fb = phase_Rxx_calib_fb - diag(diag(phase_Rxx_calib_fb));
      end
      imagesc(phase_Rxx_calib_fb*180/pi)
      h = colorbar;
      ylabel(h,'\angle{R} (\circ)')
      Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      if ~isempty(Ticks)
        h.Ticks = Ticks;
      end
    end
    
    % return
    
    
    %% Eigenvalues plots (in dB)
    if 1
      % 1) Before calibrarion
      % ----------------------
      [V_uncalib, D_uncalib] = eig(Rxx_uncalib);
      [D_uncalib, D_uncalib_idxs] = sort(real(diag(D_uncalib)),'descend');
      V_uncalib = V_uncalib(:,D_uncalib_idxs);
      
      figure(9994);clf
      subplot(121)
      stem(10*log10(D_uncalib./max(D_uncalib)),'fill','b','LineWidth',2)
      xlabel('Eigenvalue index')
      ylabel('Eigenvalue (dB)')
      title('Before calib.')
      xlim([1 Nc])
      grid on
      grid minor
      
      % 2) After calibrarion
      % ----------------------
      [V_calib, D_calib] = eig(Rxx_calib);
      [D_calib, D_calib_idxs] = sort(real(diag(D_calib)),'descend');
      V_calib = V_calib(:,D_calib_idxs);
      
      figure(9994);
      subplot(122)
      stem(10*log10(D_calib./max(D_calib)),'fill','b','LineWidth',2)
      xlabel('Eigenvalue index')
      %    ylabel('Eigenvalue (dB)')
      title('After calib.')
      xlim([1 Nc])
      grid on
      grid minor
      
      % 3) Before calibration with forward-backward averaging
      % ------------------------------------------------------
      [V_uncalib_fb, D_uncalib_fb] = eig(Rxx_uncalib_fb);
      [D_uncalib_fb, D_uncalib_fb_idxs] = sort(real(diag(D_uncalib_fb)),'descend');
      V_uncalib_fb = V_uncalib_fb(:,D_uncalib_fb_idxs);
      
      figure(9995);clf
      subplot(121)
      stem(10*log10(D_uncalib_fb./max(D_uncalib_fb)),'fill','b','LineWidth',2)
      xlabel('Eigenvalue index')
      ylabel('Eigenvalue (dB)')
      title('Before calib. & FB')
      xlim([1 Nc])
      grid on
      grid minor
      
      % 4) After calibration with forward-backward averaging
      % ------------------------------------------------------
      [V_calib_fb, D_calib_fb] = eig(Rxx_calib_fb);
      [D_calib_fb, D_calib_fb_idxs] = sort(real(diag(D_calib_fb)),'descend');
      V_calib_fb = V_calib_fb(:,D_calib_fb_idxs);
      
      figure(9995);
      subplot(122)
      stem(10*log10(D_calib_fb./max(D_calib_fb)),'fill','b','LineWidth',2)
      xlabel('Eigenvalue index')
      %    ylabel('Eigenvalue (dB)')
      title('After calib. & FB')
      xlim([1 Nc])
      grid on
      grid minor
    end
    %% Plot eigenvectors spectrum (or MUSIC pseudospectrum)
    if 1
      theta = linspace(-90,90,2048)*pi/180;
      
      % Steering vectors before calibration
      sv_params = [];
      sv_params.src.y_pc = uncalib_ypc;%y_pos{fn_idx}{doa_idx};
      sv_params.src.z_pc = uncalib_zpc;%z_pos{fn_idx}{doa_idx};
      sv_params.src.fc   = fc;
      
      sv_params.extra_error_params = [];
      sv_params.calib_params       = [];
      
      A = @(theta,param)steering_mtx(theta,sv_params);
      SVs_uncalib = A(theta,sv_params)./sqrt(Nc);
      
      % Steering vectors after calibration
      sv_params = [];
      sv_params.src.y_pc = calib_ypc;%y_pos{8}{32};
      sv_params.src.z_pc = calib_zpc;% - calib_params.calib_zpc;%z_pos{8}{32};
      sv_params.src.fc   = fc;
      
      sv_params.extra_error_params = [];
      sv_params.calib_params       = [];
      
      A = @(theta,param)steering_mtx(theta,sv_params);
      SVs_calib = A(theta,sv_params)./sqrt(Nc);
      
      % Calculate the spectrum. FB averaging doesn't change the eigenvectors
      % of the DCM, but it does change the eigenvalues.
      N            = N;
      N_calib      = N;%1;
      N_uncalib_fb = N;%1;
      N_calib_fb   = N;%1;
      
      % 1) Before calibration
      % ----------------------
      f_uncalib = 1./sum(abs(V_uncalib(:,N+1:end)' * SVs_uncalib).^2,1);
      f_uncalib = 10*log10(f_uncalib./max(abs(f_uncalib)));
      
      % 2) After calibration
      % ----------------------
      % Either calibrate SVs or measurements, BUT NOT BOTH. Calibrating the
      % measurements correct the eigenvectors, and thus the DOA estimation, and
      % this is what you need to do here.
      f_calib = 1./sum(abs(V_calib(:,N_calib+1:end)' * SVs_uncalib).^2,1);
      %   f_calib = 1./sum(abs(V_uncalib(:,N_calib+1:end)' * SVs_calib).^2,1);
      f_calib = 10*log10(f_calib./max(abs(f_calib)));
      
      % 3) Before calibration with forward-backward averaging
      % ------------------------------------------------------
      f_uncalib_fb = 1./sum(abs(V_uncalib_fb(:,N_uncalib_fb+1:end)' * SVs_uncalib).^2,1);
      f_uncalib_fb = 10*log10(f_uncalib_fb./max(abs(f_uncalib_fb)));
      
      % 4) After calibration with forward-backward averaging
      % ------------------------------------------------------
      f_calib_fb = 1./sum(abs(V_calib_fb(:,N_calib_fb+1:end)' * SVs_uncalib).^2,1);
      f_calib_fb = 10*log10(f_calib_fb./max(abs(f_calib_fb)));
      
      % DOA at maximum point in the spectrum
      [~ ,theta_uncalib_idx] = max(f_uncalib);
      theta_max_uncalib = theta(theta_uncalib_idx)*180/pi;
      
      [~, theta_max_f_calib_idx] = max(f_calib);
      theta_max_calib = theta(theta_max_f_calib_idx)*180/pi;
      
      [~, theta_uncalib_fb_idx] = max(f_uncalib_fb);
      theta_max_uncalib_fb = theta(theta_uncalib_fb_idx)*180/pi;
      
      [~, theta_calib_fb_idx] = max(f_calib_fb);
      theta_max_calib_fb = theta(theta_calib_fb_idx)*180/pi;
      
      sprintf('\nUncalib.: %2.2f deg. \n Calib.: %2.2f deg. \n Uncalib. with FB: %2.2f deg. \n Calib. with FB: %2.2f deg. \n',...
        theta_max_uncalib,theta_max_calib,theta_max_uncalib_fb,theta_max_calib_fb)
      
      figure(9996);clf
      subplot(121)
      plot(theta*180/pi,f_uncalib,'b')
      xlim(DOA*180/pi+[-60  60])
      xlabel('\theta ^\circ')
      ylabel('Pseudo-power (dB)')
      title('Before calib.')
      grid on
      grid minor
      
      subplot(122)
      plot(theta*180/pi,f_calib,'b')
      xlim(DOA*180/pi+[-60  60])
      xlabel('\theta ^\circ')
      %     ylabel('Pseudo-power (dB)')
      title('After calib.')
      grid on
      grid minor
      suptitle(sprintf('Eigenvectors pattern: Nsig=%2d',N))

      figure(9997);clf
      subplot(121)
      plot(theta*180/pi,f_uncalib_fb,'b')
      xlim(DOA*180/pi+[-60  60])
      xlabel('\theta ^\circ')
      ylabel('Pseudo-power (dB)')
      title('Before calib.&FB')
      grid on
      grid minor
      
      subplot(122)
      plot(theta*180/pi,f_calib_fb,'b')
      xlim(DOA*180/pi+[-60  60])
      xlabel('\theta ^\circ')
      %   ylabel('Pseudo-power (dB)')
      title('After calib.&FB')
      grid on
      grid minor
      suptitle(sprintf('Eigenvectors pattern: Nsig=%2d',N))
    end
    
    % return
    
    %%
    Color  = {'b','r','k','c','m','g',[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560]};
    Marker = {'-*','-o','-+','-.','-s','-^','-<','->'};
    
    if 1     
      theta_RP = linspace(-90,90,2048)*pi/180;
      [~, DOA_idx] = min(abs(DOA-theta_RP)); % For main beam direction
      
      w = hanning(Nc);
      %     w = ones(Nc,1);
      
      % Plot 3dB radiation pattern for 2 cases
      % ---------------------------------------------------------------------
      % Case 1: Before array calibration
      sv_params.src.y_pc = uncalib_ypc;%y_pos{fn_idx}{doa_idx};
      sv_params.src.z_pc = uncalib_zpc;%z_pos{fn_idx}{doa_idx};   % z_pc(:,doa_idx);
      sv_params.src.fc   = fc;
      
      sv_params.extra_error_params = [];
      sv_params.calib_params       = [];
      
      A = @(theta,param)steering_mtx(theta,sv_params);
      SVs_uncalib = A(theta_RP,sv_params)./sqrt(Nc);
      A0 = SVs_uncalib(:,DOA_idx); % Main-beam direction
      
      RP = abs((w.*A0)'*SVs_uncalib).^2;
      RP = RP./max(RP);
      RP_dB = 10*log10(RP);
      
      idxs_3dB  = find(RP_dB>=-3);
      RP_3dB_uncalib    = RP_dB(idxs_3dB);
      theta_3dB_uncalib = theta_RP(idxs_3dB);
      
      beam_doa_lims = [min(theta_3dB_uncalib) max(theta_3dB_uncalib)];
      
      RP_dB_uncalib = RP_dB;
      beam_doa_lims_vec(:,1) = beam_doa_lims;
      
      figure(100);clf
      subplot(121)
      plot(theta_3dB_uncalib*180/pi,RP_3dB_uncalib,'b')
      xlabel('\theta^\circ')
      ylabel('Power (dB)')
      title('Before calibration')
      xlim([min(theta_3dB_uncalib*180/pi)  max(theta_3dB_uncalib*180/pi)])
      %     xlim([beam_doa_lims(1) beam_doa_lims(2)]*180/pi)
      grid on
      
      % Case 2: After array calibration
      sv_params.src.y_pc = calib_ypc;
      sv_params.src.z_pc = calib_zpc;
      sv_params.src.fc   = fc;
      
      sv_params.extra_error_params = [];
      sv_params.calib_params       = [];
      
      A = @(theta,param)steering_mtx(theta,sv_params);
      SVs_calib = A(theta_RP,sv_params)./sqrt(Nc);
      A0 = SVs_calib(:,DOA_idx);
      
      RP = abs((w.*A0)'*SVs_calib).^2;
      RP = RP./max(RP);
      RP_dB = 10*log10(RP);
      
      idxs_3dB  = find(RP_dB>=-3);
      RP_3dB_calib    = RP_dB(idxs_3dB);
      theta_3dB_calib = theta_RP(idxs_3dB);
      
      beam_doa_lims = [min(theta_3dB_calib) max(theta_3dB_calib)];
      
      RP_dB_calib = RP_dB;
      beam_doa_lims_vec(:,2) = beam_doa_lims;
      
      figure(100);
      subplot(122)
      plot(theta_3dB_calib*180/pi,RP_3dB_calib,'b')
      xlabel('\theta^\circ')
      ylabel('Power (dB)')
      title('After calibration')
      xlim([min(theta_3dB_calib*180/pi)  max(theta_3dB_calib*180/pi)])
      %     xlim([min(theta_3dB*180/pi)  max(theta_3dB*180/pi)])
      %     xlim([beam_doa_lims(1) beam_doa_lims(2)]*180/pi)
      grid on
      suptitle(sprintf('Array 3dB radiation pattern of %2.2f deg. target',DOA*180/pi))
      
      close(figure(100))
      
      % Plot the three radiation patterns above in one plot
      figure(101);clf
      hold on
      plot(theta_RP*180/pi,RP_dB_calib,'b')
      plot(theta_RP*180/pi,RP_dB_uncalib,'r')
      %     xlim([beam_doa_lims_vec(1,1) beam_doa_lims_vec(2,1)]*180/pi)
      %     xlim([min(beam_doa_lims_vec(1,:)) max(beam_doa_lims_vec(2,:))]*180/pi)
      ylim([-3 0])
      xlabel('\theta^\circ')
      ylabel('Power (dB)')
      title(sprintf('Array 3dB radiation pattern of %2.2f deg. target',DOA*180/pi))
      grid on
      legend('After calibration','Before calibration','Location','best')
      
      % Plot sensors gain pattern
      % ---------------------------------------------------------------------
      % After array calibration (before calibration it was 1 or uniform)
      theta_mtx = repmat(theta_RP,[Nc 1]);
      error_g_p_mtx = repmat(calib_params.calib_g_p,[1 length(theta_RP)]);
      error_g_s_mtx = repmat(calib_params.calib_g_s,[1 length(theta_RP)]);
      error_g_offset_mtx = repmat(calib_params.calib_g_offset,[1 length(theta_RP)]);
      gain_dB = error_g_s_mtx.* (sin(theta_mtx) - sin(error_g_p_mtx)).^2 + error_g_offset_mtx; % In dB
      
      figure(102);clf
      %     subplot(121)
      hold on;
      for chan_idx = 1:Nc
        chan_resp = 10.^(gain_dB(chan_idx,:));
        %     chan_resp = SVs_calib(chan_idx,:);
        chan_gain = abs(chan_resp).^2;
        %     chan_gain = chan_gain./max(chan_gain);
        chan_gain_dB = 10*log10(chan_gain);
        if all(abs(chan_gain_dB)<1e-4)
          chan_gain_dB = zeros(size(chan_gain_dB));
        end
        plot(theta_RP*180/pi,chan_gain_dB,'Color',Color{chan_idx})
      end
      xlabel('\theta^\circ')
      ylabel('Power (dB)')
      title(sprintf('Estimated gain relative to sensor 1 (after calib.)'))
      xlim([-20  20])
      grid on
      legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
      
      % Plot sensors phase deviation
      % ---------------------------------------------------------------------
      % After array calibration (before calibration it was 0)
      phase_deg = (calib_params.calib_ypc*k*sin(theta_RP) - calib_params.calib_zpc*k*cos(theta_RP) + ...
        repmat(calib_params.calib_phase + ac_est_error(:,2)*k,[1 length(theta_RP)]))*180/pi;
      
      figure(103);clf
      hold on
      for chan_idx = 1:Nc
        % chan_phase = SV_phase_dev(chan_idx,:);
        %     chan_phase = angle(chan_resp)*180/pi;
        chan_phase = phase_deg(chan_idx,:);
        if all(abs(chan_phase)<1e-4)
          chan_phase = zeros(size(chan_phase));
        end
        plot(theta_RP*180/pi,chan_phase,'Color',Color{chan_idx})
      end
      xlabel('\theta^\circ')
      ylabel('Phase (deg.)')
      title(sprintf('Estimated phase deviation relative to sensor 1 (after calib.)'))
      xlim([-20  20])
      %   xlim([ceil(min(theta_3dB_calib)*180/pi)   floor(max(theta_3dB_calib)*180/pi)])
      grid on
      legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
      
      % Plot sensors phase pattern
      % ---------------------------------------------------------------------
      % Case 1: Before  calibration
      figure(104);clf
      subplot(211)
      hold on
      for chan_idx = 1:Nc
        chan_resp = SVs_uncalib(chan_idx,:);
        chan_phase = angle(chan_resp)*180/pi;
        plot(theta_RP*180/pi,chan_phase,'Color',Color{chan_idx})
      end
      %     xlabel('\theta^\circ')
      ylabel('Phase (deg.)')
      title('Before calibration')
      xlim([ceil(min(theta_3dB_uncalib)*180/pi)   floor(max(theta_3dB_uncalib)*180/pi)])
      %     xlim(DOA*180/pi+[-10 +10])
      grid on
      legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
      
      % Case 2: After calibration
      figure(104);
      subplot(212)
      hold on
      for chan_idx = 1:Nc
        chan_resp = SVs_calib(chan_idx,:) ;
        chan_phase = angle(chan_resp)*180/pi;
        plot(theta_RP*180/pi,chan_phase,'Color',Color{chan_idx})
      end
      xlabel('\theta^\circ')
      ylabel('Phase (deg.)')
      title('After calibration')
      xlim([ceil(min(theta_3dB_calib)*180/pi)   floor(max(theta_3dB_calib)*180/pi)])
      %     xlim(DOA*180/pi+[-10 +10])
      grid on
      %     legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
      suptitle('Phase pattern relative to sensor 1')
      
    end
    
  end
  
  if param.debug_level >= 3
    return
  end
  
  return