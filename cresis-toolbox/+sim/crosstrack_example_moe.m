% script sim.crosstrack_example_moe
%
% Example model order estimation simulation. Calls crosstrack.m and
% crosstrack_rmse.m.
%
% Author: Mohanad Al-Ibadi, Sravya Athinarapu, Sean Holloway, John Paden,
% Theresa Stumpf
%
% See also: sim.crosstrack.m, sim.crosstrack_data.m,
% sim.crosstrack_example_*.m, sim.crosstrack_rmse.m

%% Model order estimation
% =========================================================================
param = [];
physical_constants;
tic

%DECIMATION USED IN OPTIMIZER FOR NT TO AVOID CORRELATED BINS FOR TRAINED
%BUT ALL DATA IS SAVED. WHILE PLOTTING IF DECIMATED PLOTS ARE REQUIRED THEN
%IT SHOULD BE DONE IN WHERE WE ARE PLOTTING
%%
% Set norm_saved/penalty_saved to 1 if they are already generated and saved
norm_saved = 0;
penalty_saved = 0;
save_eigenvalues = 1;
% opt_norm=1 if  normalization is used (normalize optimal methods). So, if
% the normalization coefficients are not already saved, generate them.
opt_norm = ~ norm_saved;

% If suboptimal/optimal_test=1, then NT will be activated (i.e.
% optimizer_NT will be called).
suboptimal_test = 1;
optimal_test    = 1;  %running optimal methods therefore calculating normalization term for optimal methods
% 1 Set to 1 to run numerical tuning (i.e. generate the penalty coefficients)
optimizer       = 1;

% norm_allign_zero=1 is the best case.It means the reference for normalizing
% loglikelihoods is q=0 case. It is done separately for suboptimal and
% optimal.
param.norm_allign_zero = 1;

% Set layers = 1 if multiple surfaces/layers were used.
% Sravya generated all her results with layers = 0 (i.e place targets on
% range-bins, not on surfaces).
layers = 0;%1;

% For decimated results. Decimation is to eliminate frequency leakage
% This parameter is used inside otimizer_fit_NT_2D. I don't it is
% implemented correctly..so, don't use till make sure it works (i.e. set to 0 for now).
SS_decimation      = 1;
% decimation_factor: number of neighboring range-bins to be skipped to avoid
% correlation
decimation_factor = 2;
% Number of targets in a cluster. Default is 1 for sparse surface
dist_target_factor = 1;
%   dist_target_factor = 3;

% SS stands for sparse surface. So, SS=1 enables sparse surface scenario.
if dist_target_factor>1
  SS = 0;
else
  SS = 1;
end
% Angular distance, in deg, between the distributed targets in a given cluster
dist_deg           = 0.2;

% Plots control parameters
likelihood_plots = 1;
stat_test        = 1;
IMAGESC_plots    = 1;

AICc_both        = 0;

if AICc_both == 0  % NOT USED IN THIS SCRIPTS>>> CHECK WHERE TO USE
  methods = 0:6;
elseif AICc_both == 1
  methods = 0:7;
end
surf_param_all.z.mean =[];
surf_param_all.y = [];

%%
%array_proc
param_debug_all_testing = [];

%(in dB)
SNR_training_Q_0 = -inf;
SNR_training     = [10 20 30];
SNR_testing      = [10 20 30];

% ---------------------------------------------------------------------
%%         Setup simulation parameters  (TRAINING)
% ---------------------------------------------------------------------
param.suboptimal_test    = suboptimal_test ;
param.optimal_test       = optimal_test;
param.dist_target_factor = dist_target_factor;
param.SS                 = SS;

% Debug level of 3 causes this function to stop early and output
% simulated data, simulation parameters, and array processing results
% in a "results" structure.
param.debug_level = 3;

% ---------------------------------------------------------------------
%% Source parameters
% ---------------------------------------------------------------------
%narrowband
fc = 195e9;
BW = 32e6;

% % wideband
%     fc = 195e6;
%     BW = 32e6;

param.src.f0                      = fc-BW/2;
param.src.f1                      = fc+BW/2;
fs = BW;

% VARY tpd TO 3 us FOR REAL DATA GREENLAND 2014 P3
flight_h = 1500;
if 0
  % If a specific number of range-bins is required
  N_reqd_rbins = 825;
  dt = 1/BW;
  param.src.t0                      = 2*(flight_h-500)/c;
  param.src.t1                      = param.src.t0 + N_reqd_rbins*dt;
elseif 1
  param.src.t0                      = 2*(flight_h-500)/c;
  param.src.t1                      = 2*(flight_h+1000)/c;
end
param.src.ft_func                 = @(t) tukeywin_cont(t * BW);
param.src.ft_wind                 = @(N) hanning(N);
% param.src.ft_wind1                = @(N) blackman(N);
param.src.lever_arm.fh            = @sim.lever_arm_example;
% Nsep: number of lambda/4 steps between sensors
Nsep = 1;
% Nc: number of sensors
Nc = 7;
% M: max number of sources we want to estimate (it can go max upto Nc-1)
M = 6; %Nc-1

% Arguments for a linear array in y-dimension
% param.src.lever_arm.fh_args: Arguments for a linear array in y-dimension
param.src.lever_arm.fh_args       = {[], 1, 1:Nc, [0; Nsep*c/fc/2; 0]};

%NOISE POWER (dB)
% param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
param.src.noise_power             = 0*ones(1,Nc);

% DOA method parameters
param.method.list                   = [7];

% ---------------------------------------------------------------------
%% Simulation Runs Setup [For training phase]
% ---------------------------------------------------------------------
% Cross track monte carlo setup (surface generation function handle)
param.monte.target_func = @sim.surface_gen;

%RUNS
param.monte.runs = 5;

% ALTER HERE TO Q= 0 WITH DIFFERENT DATA
param.monte.random_seed_offset = 0;  %(0 same as testing data)

% Target surface parameters
surf_param                 = [];
surf_param.z.mean          = -1500;
surf_param.z.rms_height    = 0;
surf_param.z.corr_length_x = 400;
surf_param.z.corr_length_y = 400;
surf_param.rcs_in.mean     = 0;
surf_param.dy      = 10;
surf_param.y_range = [-2500 2500];
surf_param.dx      = 10;
surf_param.x_range = [-2500 2500];

% This is for Q=0 case. [] is not allowed so one target with -inf SNR used.
surf_param.y       = [0].' ;

% param.monte.target_param{1} = surf_param;

% surf_param = [];
% surf_param.z.mean = 0.7;
% surf_param.z.rms_height = 0.5;
% surf_param.z.corr_length_x = 1000;
% surf_param.z.corr_length_y = 40;
% surf_param.rcs_in.mean = 0;
% surf_param.rcs_in.var = 100;
% surf_param.dy = 1;
% surf_param.y_range = [-1500 1500];
% surf_param.dx = 1;
% surf_param.x_range = [-1500 1500];
% surf_param.x = [-1:0.1:1];
% surf_param.y = [-50:1:50].';
% param.monte.target_param{2} = surf_param;

% ---------------------------------------------------------------------
%% Array Processing parameters
% ---------------------------------------------------------------------
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

array_param.dbin  = 1;
array_param.dline = 1;

array_param.bin_rng   = 0;
array_param.line_rng = [-10:10];
Nsnap = length(array_param.bin_rng)*length(array_param.line_rng);

N_skipped_rlines = 1;
N_reqd_rlines    = 1;
surf_param.x = [-750:N_skipped_rlines:-750+N_skipped_rlines*(N_reqd_rlines+2*max(array_param.line_rng))-1];

param.monte.target_param{1} = surf_param;

%SRAVYA
array_param.Nsrc  = 1:M;
array_param.init = 'ap';
array_param.doa_theta_guard = 2/180*pi;%(max(theta)-min(theta))/(4*Nc);
array_param.Nsubband = 1;
dt = 1/fs;
array_param.imp_resp.time_vec = -3*array_param.Nsubband*dt : dt/8 : 3*array_param.Nsubband*dt;
array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / BW);

for idx = 1:max(array_param.Nsrc)
  array_param.doa_constraints(idx).method          = 'fixed';
  array_param.doa_constraints(idx).init_src_limits = [min(theta) max(theta)]*180/pi;
  array_param.doa_constraints(idx).src_limits      = [min(theta) max(theta)]*180/pi;
end

param.array_param = array_param;
clear array_param;

% ----------------------------------------------------------------------------------
%% Training (Q=0): generated normalizing term of the loglikelihood function
% ----------------------------------------------------------------------------------
% Here, Q=0 doesn't mean no targets. We still have 0:M targets, but the
% SNR=-inf, so as if there is only noise (no signal, no targets).
param.testing   = 0;
% Training SNR for the case q=0 (usually -inf)
param.SNR_db    = SNR_training_Q_0;
param.opt_norm = opt_norm;

if ~norm_saved
  % Normalization terms are not yet generated. So, generate them here.
  
  % MOHANAD: Generate log-likelihoods for all SNRs and runs.
  param.moe_methods = methods(1); % For training, only method=0 (NT) is needed
  LL_results = sim.crosstrack(param);
  LL_subopt = LL_results.LL_subopt;
  LL_opt    = LL_results.LL_opt;
  
  % Mohanad: Determine the log-likelihood normilizatin coefficiens
  LL_subopt_tmp = [];
  LL_opt_tmp    = [];
  for run_idx = 1:length(LL_subopt)
    LL_subopt_tmp = [LL_subopt_tmp;LL_subopt{run_idx}{1}];
    LL_opt_tmp    = [LL_opt_tmp;LL_opt{run_idx}{1}];
  end
  LL_subopt_mean = nanmean(LL_subopt_tmp,1);
  LL_opt_mean = nanmean(LL_opt_tmp,1);
  
  norm_coeff.opt_norm_term        = LL_subopt_mean - LL_opt_mean ;
  norm_coeff.norm_term_suboptimal = - LL_subopt_mean;
  norm_coeff.norm_term_optimal    = - LL_opt_mean ;
  
  if param.opt_norm ==1
    opt_norm_term = norm_coeff.opt_norm_term;
  else
    opt_norm_term = zeros(1,max(param.array_param.Nsrc )+1);
  end
  
  if param.norm_allign_zero ==1
    norm_term_suboptimal = norm_coeff.norm_term_suboptimal;
    norm_term_optimal    = norm_coeff.norm_term_optimal;
  else
    norm_term_suboptimal =  zeros(1,max(param.array_param.Nsrc )+1);
    norm_term_optimal    = zeros(1,max(param.array_param.Nsrc )+1);
  end
else
  % Load the normalization parameters if you already have them.
end

param.opt_norm_term        = opt_norm_term;
param.norm_term_suboptimal = norm_term_suboptimal;
param.norm_term_optimal    = norm_term_optimal;

param.opt_norm = ~opt_norm;  % since section is done.

% ********* Generating normalization coefficints is done here *************
% *************************************************************************

%% Place targets on the surface at specific locations/angles
% -------------------------------------------------------------------------
param.optimizer = optimizer;  % 1 FOR NUMERICAL TUNING

% ALTER HERE TO TRAIN WITH DIFFERENT DATA
param.monte.random_seed_offset = 1;

if layers ==0
  %% Targets in range cylinders (targets are spaced along range-bins, not surfaces)
  dt = 1/fs;
  Nt = floor((param.src.t1-param.src.t0)/dt);
  time = param.src.t0 + dt*(0:Nt-1).';
  R_bins_values = time*c/2;
  
  R_shell_values = [R_bins_values;R_bins_values(end)+(c/(2*BW))];
  
  R_shell_mid_pt = R_bins_values + (1/2)*(c/(2*BW));
  
  clear N_bins_Q
  
  for source_idx = 0:M
    N_bins_Q(:,source_idx+1) = source_idx*ones(floor(Nt/(M+1)),1) ;
  end
  N_bins_Q = N_bins_Q(:);
  N_bins_Q =[N_bins_Q; zeros((Nt-numel(N_bins_Q)),1)];
  
  doa_bins = NaN*ones(Nt,max(param.array_param.Nsrc)*dist_target_factor);
  y_bin    = NaN*ones(Nt,max(param.array_param.Nsrc)*dist_target_factor);
  z_bin    = NaN*ones(Nt,max(param.array_param.Nsrc)*dist_target_factor);
  
  % ORTHOGONAL SV
  p = Nc;
  d = lambda/2;
  L = p*d;
  
  k_spacing = [-(p-1)/2:1:(p-1)/2]*(lambda/L); % as implemented in paper for orthogonal sv
  %k_spacing = [-p/2:1:(p/2)-1]*(lambda/L); % John
  DOA_orthogonal = asin(k_spacing)*180/pi;
  
  if 0
    % Debug: Make sure SVs are orhogonal
    for k_spacing_idx = 1:length(k_spacing)
      %keyboard;
      sv(:,k_spacing_idx) = exp(1i*pi*sind(DOA_orthogonal(k_spacing_idx))*(0:p-1)).';
    end
    for k_spacing_idx0= 1:length(k_spacing)
      for k_spacing_idx = 1:length(k_spacing)
        dot_result(k_spacing_idx0,k_spacing_idx) = dot(sv(:,k_spacing_idx0),sv(:,k_spacing_idx));
      end
    end
    real(dot_result);
    imag(dot_result);
  end
  
  for N_bins_idx = 1:length(N_bins_Q)
    source_idx = N_bins_Q(N_bins_idx);
    if 1
      index_ref = ceil(p/2);
      % +DOA and -DOA are getting switched (now fixed)
      if source_idx == 0
        % angs_elecs: 1 x N matrix, angle of arrival for each source (deg)
        q = [];
      else
        if mod(source_idx,2)==1 %odd
          index =((index_ref)-(source_idx-1)/2):1:(index_ref)+(source_idx-1)/2;
        else
          index =((index_ref)-(source_idx-2)/2):1:(index_ref)+(source_idx-2)/2;
          index = [index index(end)+1];
        end
        q = DOA_orthogonal(index);
        clear index
      end
    else  % using both +DOA and -DOA (even if they are switched no problem)
      if source_idx == 0
        % angs_elecs: 1 x N matrix, angle of arrival for each source (deg)
        q = [];
      else
        if mod(source_idx,2)==1 %odd
          
          index =((index_ref)-(source_idx-1)/2):1:(index_ref)+(source_idx-1)/2;
        else
          index =((index_ref)-(source_idx)/2):1:(index_ref)+(source_idx)/2;
          index(find(index==index_ref))= [];
        end
        q = DOA_orthogonal(index);
        clear index
      end
      
    end
    
    % w.r.t chan = 4;
    if ~isempty(q)
      if SS==0
        q = [q-dist_deg q q+dist_deg];
        %                 if dist_target_doa ==1
        %                     q = [q-dist_deg q q+dist_deg];%FILLED SURFACE WITH CLOSE DOA
        %                 end
      end
      
      z_bin(N_bins_idx,1:length(q)) = -1*sqrt((R_shell_mid_pt(N_bins_idx,:)^2)./(tand(q).^2 +1) );% -1 multiplied as sqrt always return only positive
      surf_param_all.z.mean = [surf_param_all.z.mean    z_bin(N_bins_idx,1:length(q))];
      
      y_bin(N_bins_idx,1:length(q)) = z_bin(N_bins_idx,1:length(q)).*tand(q);
      
      surf_param_all.y = [surf_param_all.y    y_bin(N_bins_idx,1:length(q))];
    end
    doa_bins(N_bins_idx,1:length(q)) = q;
  end
  
  if 0 && SS==0
    % This is used if we would like to add more targets. So, the total
    % number of targets becomes:
    % length(surf_param_all.y)*dist_target_factor*3.
    temp_y = surf_param_all.y;
    clear surf_param_all.y
    surf_param_all.y = [temp_y+4  temp_y temp_y-4];
    
    temp_z = surf_param_all.z.mean;
    clear surf_param_all.z.mean
    surf_param_all.z.mean = [surf_param_all.z.mean surf_param_all.z.mean surf_param_all.z.mean];
  end
  
  for y_idx = 1: length(surf_param_all.y)
    y_temp{y_idx} = surf_param_all.y(y_idx);
  end
  
  surf_param_all.y = y_temp;
  clear y_temp
else
  if 0
    surf_param_all.z.mean = [-1500 -1650 -1950];
    surf_param_all.y{1} = [-1430 -1250 -987 -860 -695 -420  0 420 695 860 987 1250 1430].';
    surf_param_all.y{2} = [-1505 -1250 -1045 -860 -710 -510 0 510 710 860 1045 1250 1505].';
    surf_param_all.y{3} = [-1240 -700 -510 0 510 700 1240 ].'  ;
  elseif 0
    surf_param_all.z.mean = [-1500];
    surf_param_all.y{1} = [-1430 -1250 -987 -860 -695 -420  0 420 695 860 987 1250 1430].';
  elseif 1
    surf_param_all.z.mean = -flight_h;
    DOA = [-[17 15 13 11 9 7 5 3], [0 3 5 7 9 11 13 15 17]]'*pi/180;
    surf_param_all.y{1} = flight_h*tan(DOA);
  end
end

%%         Training: generated penalty coefficients (NT)
% -------------------------------------------------------------------------
% P_md_coeff is a scalar to control Probability of missed detection, P_md, (underestimation)
%If P_md is very low, then false alarms will be higher?. i.e. we?ll overestimate a lot.
% Alternatively, the probability of false alarm, P_fa, can be set (overestimation),
% but not both P_md and P_fa together since they depend on one another.
% If both are set to 1, that means no overestimation and underestimation.
P_md_coeff = 1;
P_fa_coeff = 1;

if penalty_saved ==0
  param.moe_methods = methods;
  LL_results = [];
  LL_subopt  = [];
  LL_opt     = [];
  if param.optimizer ==1
    param.SNR_db = SNR_training;
    for surf_idx = 1:length(surf_param_all.y)
      surf_param.y = surf_param_all.y{surf_idx};
      surf_param.z.mean =  surf_param_all.z.mean(surf_idx);
      param.monte.target_param{surf_idx} = surf_param;
    end
    
    % Generate (normalized) log-likelihoods for all {SNRs} and {runs}.
    LL_results = sim.crosstrack(param);
    LL_subopt = LL_results.LL_subopt;
    LL_opt    = LL_results.LL_opt;
    eigenvalues_all = LL_results.eigenvalues_all;
    actual_num_targets = LL_results.actual_num_targets;
    
    % Mohanad: Pass data and parameters to the optimizer
    optimizer_NT_param.actual_num_targets = actual_num_targets;
    
    optimizer_NT_param.Nsnap           = Nsnap;
    optimizer_NT_param.SNR_db          = param.SNR_db;
    optimizer_NT_param.Nsrc            = M;
    optimizer_NT_param.Nc              = Nc;
    optimizer_NT_param.P_md_coeff      = P_md_coeff;
    optimizer_NT_param.P_fa_coeff      = P_fa_coeff;
    optimizer_NT_param.SS_decimation = SS_decimation;
    optimizer_NT_param.decimation_factor = decimation_factor;
    
    if suboptimal_test
      optimizer_NT_param.log_func_all = LL_subopt;
      param.NT = optimizer_NT_2D(optimizer_NT_param);
      %         param.NT = sim.optimizer_NT_2D(optimizer_NT_param);
    end
    
    if optimal_test
      optimizer_NT_param.log_func_all = LL_opt;
      param.NT_opt   = optimizer_NT_2D(optimizer_NT_param);
      %         param.NT_opt   = sim.optimizer_NT_2D(optimizer_NT_param);
    end
  else
    param.NT = zeros(1,M+1);
    param.NT_opt = zeros(1,M+1);
  end
else
  % If penalty coefficients already saved, then load them here.
end

param.optimize = 0;   % since NT section is done.

% ********* Generating penalty coefficints is done here *************
% *************************************************************************

%% Store DCMs and eigenvalues from all runs and SNRs, as well as simulation parameters
if save_eigenvalues
  sim_param.fc           = fc;
  sim_param.BW           = BW;
  sim_param.Nc           = Nc;
  sim_param.Nsnap        = Nsnap;
  sim_param.M            = M;
  sim_param.SNR_training = param.SNR_db;
  %     sim_param.phase_center = phase_center;
  %     sim_param.Nb = Nb;
  sim_param.Nruns = param.monte.runs;
  sim_param.notes{1} ='LL_opt, LL_subopt,eigenvalues_all, and actual_num_targets have the following forms: {run_idx}{snr_idx}. The first 3 have dimension of NtNx-by-Nc, while the last one has a dimension of NtNx-by-1.';
  
  penalty.opt    = param.NT_opt;
  penalty.subopt = param.NT;
  
  normalization.norm_allign_zero     = param.norm_allign_zero;
  normalization.norm_term_optimal    =  norm_coeff.norm_term_optimal;
  normalization.norm_term_suboptimal = norm_coeff.norm_term_suboptimal;
  normalization.opt_norm_term        = norm_coeff.opt_norm_term;
  
  %     out_fn_dir = '/users/mohanad/IceSheetProject/MOE work/DCMandEigenvalues/';
  out_fn_dir = 'H:\IceSheetProject\MOE work\DCMandEigenvalues\';
  out_fn_name = '2D';
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  save([out_fn '.mat'],'eigenvalues_all','actual_num_targets','LL_opt','LL_subopt','sim_param','penalty','normalization')
end

%%                      TESTING
% -------------------------------------------------------------------------
param.testing = 1;  % to run array_proc from crosstrack

if 1
  fprintf('=======================================================================\n');
  fprintf('Running sim.crosstrack_example example #2\n');
  fprintf('  Narrowband radar depth sounder simulation\n');
  fprintf('  Linear array in y-dimension\n');
  fprintf('=======================================================================\n');
  
  param.monte.random_seed_offset = 2;
  
  % Target surface parameters
  surf_param = [];
  surf_param.z.rms_height = 0;
  surf_param.z.corr_length_x = 400;
  surf_param.z.corr_length_y = 400;
  surf_param.rcs_in.mean = 0;
  
  surf_param.dy = 10;
  surf_param.y_range = [-2500 2500];
  surf_param.dx = 10;
  surf_param.x_range = [-2500 2500];
  %     N_skipped_rlines = 5; % Same as for training
  %     N_reqd_rlines    = 1; % Same as for training
  surf_param.x = [-750:N_skipped_rlines:-750+N_skipped_rlines*(N_reqd_rlines+2*max(param.array_param.line_rng))-1];
  
  param.SNR_db = SNR_testing; % SNR TRAINING cases
  
  if 0  %CHECK AGAIN
    noise_power = 10.^(param.src.noise_power(1)/10)
    SNR_db = SNR_testing  % SNR TESTING%
    surf_param.rcs_in.var = noise_power * 10.^(SNR_db./10)   % CHECK
  end
  
  for surf_idx = 1:length(surf_param_all.y)
    surf_param.y = surf_param_all.y{surf_idx};
    surf_param.z.mean =  surf_param_all.z.mean(surf_idx);
    param.monte.target_param{surf_idx} = surf_param;
  end
  
  %% Run the simulation
  param.suboptimal_test = 1;
  param.optimal_test    = 1;
  %     param.moe_methods = 0;
  param.moe_methods = methods;
  %     param.NT = [];
  %     param.NT_opt    = [2.3979  164.8862  251.2631];
  
  % fcs.surface is used in array_proc to skip range-bins above the surface.
  % Do not use when training MOE.
  if 0 && (isfield(param,'testing') && ~isempty(param.testing) && param.testing==1) ...
      || (~isfield(param,'testing'))
    param.array_param.surface = flight_h*ones(length(surf_param.x),1)./c;
  end
  
  results = sim.crosstrack(param);
  actual_num_targets = results.actual_num_targets;
  %     actual_doa_targets = results.actual_doa_targets;
  
  if 0
    % Plot sample slice. Used to test MLE vs S-MLE
    dout = results.tomo{1}{1};
    doa_NT = dout.model_order_results_optimal.NT.doa*180/pi;
    actual_doa = actual_doa_targets{1}{3}*180/pi;
    Nt = size(doa_NT,1);
    
    figure(1116);clf
    hold on
    p1 = plot(doa_NT(:,1),1:Nt,'*b');
    plot(doa_NT(:,2),1:1:Nt,'*b')
    
    p2 = plot(actual_doa(:,1),1:Nt,'*r');
    plot(actual_doa(:,2),1:Nt,'*r')
    
    set(gca,'Ydir','reverse')
    xlabel('\theta^\circ')
    ylabel('Range-bin index')
    ylim([30 45])
    grid on
    title('MLE -- Numerical tuning -- Optimal MOE -- 10dB SNR')
    legend([p1 p2],{'Estimated surface','Actual surface'},'Location','northeast')
  end
  
  %%                       Decimation
  % ----------------------------------------------------------------------
  if 0& SS_decimation
    for snr_idx = 1:length(SNR_testing)
      for run_idx = 1:length(actual_num_targets)
        log_func_all_tmp_new    = [];
        q_actual_rline_all      = [];
        LL_opt_tmp_rline_all    = [];
        LL_subopt_tmp_rline_all = [];
        
        q_actual      = actual_num_targets{run_idx}{snr_idx};
        dout          = results.tomo{run_idx}{snr_idx};
        LL_subopt_tmp = LL_subopt{run_idx}{snr_idx};
        LL_opt_tmp    = LL_opt{run_idx}{snr_idx};
        
        Nt = size(q_actual,1);
        Nx = size(q_actual,2);
        for line_idx = 1:Nx
          q_actual_rline = q_actual(:,line_idx);
          % Find the indices of the range-bins where Nsrc>0.
          idx = find(q_actual_rline>0);
          idx = idx(idx>decimation_factor & idx <= (Nt-decimation_factor));
          % Set the numbe of targets in decimation_factor neighboring
          % range-bins to NaN temporarily
          for i = 1:length(idx)
            n = idx(i);
            if ~isnan(q_actual_rline(n))
              q_subset = q_actual_rline([n-decimation_factor:n-1, n+1:n+decimation_factor]);
              bad_q_idx = find(q_subset>0);
              q_subset(bad_q_idx) = NaN;
              q_actual_rline([n-decimation_factor:n-1, n+1:n+decimation_factor]) = q_subset;
            end
          end
          
          % Ignore log-likelihoods of affected range-bins
          LL_subopt_tmp_rline = LL_subopt_tmp(1+(line_idx-1)*Nt:line_idx*Nt,:);
          LL_subopt_tmp_rline(isnan(q_actual_rline),:) = [];
          LL_subopt_tmp_rline_all = [LL_subopt_tmp_rline_all;LL_subopt_tmp_rline];
          
          LL_opt_tmp_rline = LL_opt_tmp(1+(line_idx-1)*Nt:line_idx*Nt,:);
          LL_opt_tmp_rline(isnan(q_actual_rline),:) = [];
          LL_opt_tmp_rline_all = [LL_opt_tmp_rline_all;LL_opt_tmp_rline];
          
          % Set estimated Nsrc of affected range-bins to NaN (will
          % be ignored later)
          for method_idx = 0:6
            switch method_idx
              case 0
                dout.model_order_results_suboptimal.NT.Nest(isnan(q_actual_rline),line_idx)    = NaN;
                dout.model_order_results_optimal.NT.Nest(isnan(q_actual_rline),line_idx)       = NaN;
              case 1
                dout.model_order_results_suboptimal.AIC.Nest(isnan(q_actual_rline),line_idx)   = NaN;
                dout.model_order_results_optimal.AIC.Nest(isnan(q_actual_rline),line_idx)      = NaN;
              case 2
                dout.model_order_results_suboptimal.HQ.Nest(isnan(q_actual_rline),line_idx)    = NaN;
                dout.model_order_results_optimal.HQ.Nest(isnan(q_actual_rline),line_idx)       = NaN;
              case 3
                dout.model_order_results_suboptimal.MDL.Nest(isnan(q_actual_rline),line_idx)   = NaN;
                dout.model_order_results_optimal.MDL.Nest(isnan(q_actual_rline),line_idx)      = NaN;
              case 4
                dout.model_order_results_suboptimal.AICc.Nest(isnan(q_actual_rline),line_idx)  = NaN;
                dout.model_order_results_optimal.AICc.Nest(isnan(q_actual_rline),line_idx)     = NaN;
              case 5
                dout.model_order_results_suboptimal.KICvc.Nest(isnan(q_actual_rline),line_idx) = NaN;
                dout.model_order_results_optimal.KICvc.Nest(isnan(q_actual_rline),line_idx)    = NaN;
              case 6
                dout.model_order_results_suboptimal.WIC.Nest(isnan(q_actual_rline),line_idx)   = NaN;
                dout.model_order_results_optimal.WIC.Nest(isnan(q_actual_rline),line_idx)      = NaN;
              otherwise
                error('Not supported')
            end
          end
          
          % Set actual Nsrc of affected range-bins to NaN (will
          % be ignored later)
          actual_num_targets{run_idx}{snr_idx}(isnan(q_actual_rline),line_idx) = NaN;
        end
        
        LL_subopt{run_idx}{snr_idx} = LL_subopt_tmp_rline_all;
        LL_opt{run_idx}{snr_idx}    = LL_opt_tmp_rline_all;
        
        results.tomo{run_idx}{snr_idx}.model_order_results_suboptimal = dout.model_order_results_suboptimal;
        results.tomo{run_idx}{snr_idx}.model_order_results_optimal    = dout.model_order_results_optimal;
      end
    end
  end
  
  %%                      Plot results
  % ---------------------------------------------------------------------
  % LOG_LIKELIHOOD PLOTS
  % --------------------
  if likelihood_plots == 1
    warning('off','MATLAB:legend:IgnoringExtraEntries')
    
    if opt_norm == 1
      figure(3);clf;
      for SNR_idx = 1:length(SNR_training_Q_0)
        subplot(length(SNR_training_Q_0),1,SNR_idx)
        plot(0:M,LL_subopt_mean,'b-*')%CHECK HOW MANY RLINES DO WE WANTED TO CONSIDER
        hold on
        plot(0:M,LL_opt_mean,'r--*' )
        
        xlabel('k','interpreter','none')
        ylabel('-2L','interpreter','none')
        %                 set(gca,'TickLabelInterpreter','Latex')
        h_legend = legend('Suboptimal','Optimal','Location','best');
        set(h_legend,'Interpreter','Latex')
        title([  num2str(SNR_training_Q_0(SNR_idx))'' ' dB,    Q=0 '],'interpreter','Latex'  )
        grid on
      end
    end
    
    if optimizer==1
      %TRAINING DATA
      figure(1);clf
      figure(2);clf
      
      Color         = {'b','r','k','c','m','g',[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560]};
      Marker_opt    = {'-*','-o','-+','-.','-s','-^','-<','->'};
      Marker_subopt = {'--*','--o','--+','--.','--s','--^','--<','-->'};
      
      for SNR_idx = 1: length(SNR_training)
        
        % PLOT LOG-LIKELIHOOD (NOT COST) WITH NORMALIZATION:
        % TRAINING
        %--------------------------------------------------
        LL_opt_tmp    = [];
        LL_subopt_tmp = [];
        actual_num_targets_tmp = [];
        for run_idx = 1:length(LL_opt)
          LL_opt_tmp    = [LL_opt_tmp;LL_opt{run_idx}{SNR_idx}]; % Nrows*Nruns-by-(M+1)
          LL_subopt_tmp = [LL_subopt_tmp;LL_subopt{run_idx}{SNR_idx}]; % Nrows*Nruns-by-(M+1)
          actual_num_targets_tmp = [actual_num_targets_tmp;actual_num_targets{run_idx}{SNR_idx}];
        end
        actual_num_targets_tmp(isnan(actual_num_targets_tmp)) = [];
        
        % Plot optimal methods
        figure(1);
        subplot(length(SNR_training),2,2*SNR_idx-1)
        hold on
        for idx = 0:M
          plot(0:M,nanmean(LL_opt_tmp(actual_num_targets_tmp==idx,:),1),Marker_opt{idx+1},'Color',Color{idx+1})
        end
        grid on
        
        % Plot suboptimal methods
        figure(2);
        subplot(length(SNR_training),2,2*SNR_idx-1)
        hold on
        for idx = 0:M
          plot(0:M,nanmean(LL_subopt_tmp(actual_num_targets_tmp==idx,:),1),Marker_subopt{idx+1},'Color',Color{idx+1})
        end
        grid on
        
        if SNR_idx ==1
          figure(1);
          title([  num2str(SNR_training(SNR_idx))'' ' dB      --- Optimal  Training (NT)'],'interpreter','Latex'  )
          
          figure(2);
          title([  num2str(SNR_training(SNR_idx))'' ' dB    -- -- Suboptimal      Training (NT)'],'interpreter','Latex'  )
        else
          figure(1);
          title([  num2str(SNR_training(SNR_idx))'' ' dB '],'interpreter','Latex'  )
          
          figure(2);
          title([  num2str(SNR_training(SNR_idx))'' ' dB '],'interpreter','Latex'  )
        end
        
        figure(1);
        ylabel('-2L','interpreter','none')
        
        figure(2);
        ylabel('-2L','interpreter','none')
        
        if SNR_idx == length(SNR_training)
          figure(1);
          xlabel('k','interpreter','none')
          h_legend = legend('q = 0','q = 1','q = 2','q = 3','q = 4','q = 5','q = 6','Location','best');
          
          figure(2);
          xlabel('k','interpreter','none')
          h_legend = legend('q = 0','q = 1','q = 2','q = 3','q = 4','q = 5','q = 6','Location','best');
        end
        hold off
        
        % PLOT LOG-LIKELIHOOD (NOT COST) WITHOUT NORMALIZATION:
        % TRAINING
        %--------------------------------------------------
        % Plot optimal methods
        figure(1);
        subplot(length(SNR_training),2,2*SNR_idx)
        hold on
        for idx = 0:M
          if param.norm_allign_zero ==1
            plot(0:M,nanmean(LL_opt_tmp(actual_num_targets_tmp==idx,:)-repmat(norm_term_optimal,[size(LL_opt_tmp(actual_num_targets_tmp==idx,:),1),1]),1),Marker_opt{idx+1},'Color',Color{idx+1})
          else
            plot(0:M,nanmean(LL_opt_tmp(actual_num_targets_tmp==idx,:)-repmat(opt_norm_term,[size(LL_opt_tmp(actual_num_targets_tmp==idx,:),1),1]),1),Marker_opt{idx+1},'Color',Color{idx+1})
          end
        end
        grid on
        
        % Plot suboptimal methods
        figure(2);
        subplot(length(SNR_training),2,2*SNR_idx)
        hold on
        for idx = 0:M
          if param.norm_allign_zero ==1
            plot(0:M,nanmean(LL_subopt_tmp(actual_num_targets_tmp==idx,:)-repmat(norm_term_suboptimal,[size(LL_opt_tmp(actual_num_targets_tmp==idx,:),1),1]),1),Marker_subopt{idx+1},'Color',Color{idx+1})
          else
            plot(0:M,nanmean(LL_subopt_tmp(actual_num_targets_tmp==idx,:),1),Marker_subopt{idx+1},'Color',Color{idx+1})
          end
        end
        grid on
        
        if SNR_idx == length(SNR_training)
          figure(1);
          xlabel('k','interpreter','none')
          
          figure(2);
          xlabel('k','interpreter','none')
        end
        
        if SNR_idx ==1
          figure(1);
          title([  num2str(SNR_training(SNR_idx))'' ' dB    Without Normalization'],'interpreter','Latex'  )
          
          figure(2);
          title([  num2str(SNR_training(SNR_idx))'' ' dB    Without Normalization'],'interpreter','Latex'  )
        else
          figure(1);
          title([  num2str(SNR_training(SNR_idx))'' ' dB'],'interpreter','Latex'  )
          
          figure(2);
          title([  num2str(SNR_training(SNR_idx))'' ' dB'],'interpreter','Latex'  )
        end
        hold off
      end
    end
  end
  
  
  % HISTOGRAM FOR A METHOD W.R.T MODEL ORDER
  % ----------------------------------------
  figure(4),clf;
  figure(5),clf;
  
  AVG_subopt = 0;
  AVG_opt    = 0;
  Nest_2D_subopt_all = [];
  Nest_2D_opt_all    = [];
  for SNR_idx = 1:length(SNR_testing)
    for method_idx = 0:length(methods)-1
      Nest_subopt = [];
      Nest_opt    = [];
      switch method_idx
        case 0
          % NT
          for run_idx = 1:length(results.tomo)
            dout = results.tomo{run_idx}{SNR_idx};
            Nest_subopt = [Nest_subopt;dout.model_order_results_suboptimal.NT.Nest];
            Nest_opt    = [Nest_opt;dout.model_order_results_optimal.NT.Nest];
          end
        case 1
          % AIC
          for run_idx = 1:length(results.tomo)
            dout = results.tomo{run_idx}{SNR_idx};
            Nest_subopt = [Nest_subopt;dout.model_order_results_suboptimal.AIC.Nest];
            Nest_opt    = [Nest_opt;dout.model_order_results_optimal.AIC.Nest];
          end
        case 2
          % HQ
          for run_idx = 1:length(results.tomo)
            dout = results.tomo{run_idx}{SNR_idx};
            Nest_subopt = [Nest_subopt;dout.model_order_results_suboptimal.HQ.Nest];
            Nest_opt    = [Nest_opt;dout.model_order_results_optimal.HQ.Nest];
          end
        case 3
          % MDL
          for run_idx = 1:length(results.tomo)
            dout = results.tomo{run_idx}{SNR_idx};
            Nest_subopt = [Nest_subopt;dout.model_order_results_suboptimal.MDL.Nest];
            Nest_opt    = [Nest_opt;dout.model_order_results_optimal.MDL.Nest];
          end
        case 4
          % AICc
          for run_idx = 1:length(results.tomo)
            dout = results.tomo{run_idx}{SNR_idx};
            Nest_subopt = [Nest_subopt;dout.model_order_results_suboptimal.AICc.Nest];
            Nest_opt    = [Nest_opt;dout.model_order_results_optimal.AICc.Nest];
          end
        case 5
          % KICvc
          for run_idx = 1:length(results.tomo)
            dout = results.tomo{run_idx}{SNR_idx};
            Nest_subopt = [Nest_subopt;dout.model_order_results_suboptimal.KICvc.Nest];
            Nest_opt    = [Nest_opt;dout.model_order_results_optimal.KICvc.Nest];
          end
        case 6
          % WIC
          for run_idx = 1:length(results.tomo)
            dout = results.tomo{run_idx}{SNR_idx};
            Nest_subopt = [Nest_subopt;dout.model_order_results_suboptimal.WIC.Nest];
            Nest_opt    = [Nest_opt;dout.model_order_results_optimal.WIC.Nest];
          end
        otherwise
          error('Not supported')
      end
      
      Nest_subopt = Nest_subopt(:);
      Nest_subopt(isnan(Nest_subopt)) = [];
      Nest_opt    = Nest_opt(:);
      Nest_opt(isnan(Nest_opt)) = [];
      
      actual_num_targets_tmp = [];
      for run_idx = 1:length(results.tomo)
        actual_num_targets_tmp = [actual_num_targets_tmp;actual_num_targets{run_idx}{SNR_idx}];
      end
      actual_num_targets_tmp = actual_num_targets_tmp(:);
      actual_num_targets_tmp(isnan(actual_num_targets_tmp)) = [];
      
      percentage_correct_subopt = [];
      percentage_correct_opt    = [];
      for k_idx = 0:M
        if isempty(actual_num_targets_tmp == k_idx)
          percentage_correct_subopt(k_idx+1) = NaN;
          percentage_correct_opt(k_idx+1)    = NaN;
        else
          actual_Nsrc_idxs = find(actual_num_targets_tmp == k_idx);
          
          Nest_subopt_subset = Nest_subopt(actual_Nsrc_idxs);
          Nest_opt_subset    = Nest_opt(actual_Nsrc_idxs);
          
          Nest_subopt_good = Nest_subopt_subset(Nest_subopt_subset==k_idx);
          Nest_opt_good    = Nest_opt_subset(Nest_opt_subset==k_idx);
          
          percentage_correct_subopt(k_idx+1) =  numel(Nest_subopt_good)/numel(actual_num_targets_tmp(actual_Nsrc_idxs)) * 100;
          percentage_correct_opt(k_idx+1) =  numel(Nest_opt_good)/numel(actual_num_targets_tmp(actual_Nsrc_idxs)) * 100;
        end
      end
      
      %performance can't be tested for the below case.
      %since no number of sorces equal the number we are testing in true sources
      %percentage(find(~isnan(percentage) == 1)) = ? ;
      
      percentage_method_subopt(:,method_idx+1) = percentage_correct_subopt;
      percentage_method_opt(:,method_idx+1)    = percentage_correct_opt;
      %
      percentage_method_subopt_all{SNR_idx} = percentage_method_subopt;
      percentage_method_opt_all{SNR_idx}    = percentage_method_opt;
      
      Nest_2D_subopt_all(:,method_idx+1) = Nest_subopt;
      Nest_2D_opt_all(:,method_idx+1)    = Nest_opt;
    end
    
    if stat_test ==1
      % Plot suboptimal methods MOE results: testing
      figure(4);
      subplot(2,2,SNR_idx)
      
      plot(0:M , percentage_method_subopt(:,2),'m'),hold on
      plot(0:M , percentage_method_subopt(:,3),'-.r')
      plot(0:M , percentage_method_subopt(:,4),'-.g')
      plot(0:M , percentage_method_subopt(:,5),'--b')
      if   AICc_both ==1
        plot(0:M , percentage_method_subopt(:,8),'-ob')
      end
      plot(0:M ,percentage_method_subopt(:,6),'--c.')
      plot(0:M , percentage_method_subopt(:,7),'-*k')
      
      plot(0:M , percentage_method_subopt(:,1),'-*g')
      
      xlim([0 M])
      ylim([0 100])
      grid on
      %             set(gca,'TickLabelInterpreter','Latex')
      if SNR_idx == length(SNR_testing)
        xlabel('Number of sources, q','interpreter','none')
        ylabel('% Correct','interpreter','none')
        %                 set(gca,'TickLabelInterpreter','Latex')
        
        if AICc_both ==1
          h_legend = legend('AIC', 'HQ', 'MDL', 'AICc', 'AICc 19', 'KICvc', 'WIC','NT','Location','best');
        else
          h_legend = legend('AIC', 'HQ', 'MDL', 'AICc', 'KICvc', 'WIC','NT','Location','best');
          
        end
        set(h_legend,'Interpreter','Latex')
      end
      
      if  SNR_idx ==1
        title([  num2str(SNR_testing(SNR_idx))'' ' dB '],'interpreter','Latex'  )
      else
        title([  num2str(SNR_testing(SNR_idx))'' ' dB'],'interpreter','Latex'  )
      end
      
      figure(5);
      subplot(2,2,SNR_idx)
      
      plot(0:M , percentage_method_opt(:,2),'m'), hold on
      plot(0:M , percentage_method_opt(:,3),'-.r')
      plot(0:M , percentage_method_opt(:,4),'-.g')
      plot(0:M , percentage_method_opt(:,5),'--b')
      if   AICc_both ==1
        plot(0:M , percentage_method_opt(:,8),'-ob')
      end
      plot(0:M ,percentage_method_opt(:,6),'--c.')
      plot(0:M , percentage_method_opt(:,7),'-*k')
      plot(0:M , percentage_method_opt(:,1),'-*g')
      
      xlim([0 M])
      ylim([0 100])
      grid on
      if SNR_idx == length(SNR_testing)
        xlabel('Number of sources, q','interpreter','none')
        ylabel('% Correct','interpreter','none')
        %                 set(gca,'TickLabelInterpreter','Latex')
        
        if AICc_both ==1
          h_legend = legend('AIC', 'HQ', 'MDL', 'AICc', 'AICc 19', 'KICvc', 'WIC','NT','Location','best');
        else
          h_legend = legend('AIC', 'HQ', 'MDL', 'AICc', 'KICvc', 'WIC','NT','Location','best');
        end
        set(h_legend,'Interpreter','Latex')
      end
      
      if  SNR_idx ==1
        title([  num2str(SNR_testing(SNR_idx))'' ' dB '],'interpreter','Latex'  )
      else
        title([  num2str(SNR_testing(SNR_idx))'' ' dB '],'interpreter','Latex'  )
      end
    end
    
    % add the actual number of targets to the Nest matrix in the last
    % column.
    method_idx = method_idx+1;
    actual_num_targets_tmp(actual_num_targets_tmp>M) = NaN;
    Nest_2D_subopt_all(:,method_idx+1) = actual_num_targets_tmp;
    Nest_2D_opt_all(:,method_idx+1)   = actual_num_targets_tmp;
    
    
    %        IMAGESC PLOTS
    % ---------------------------
    if IMAGESC_plots ==1
      figure(SNR_idx+210);clf
      if AICc_both==1
        imagesc(Nest_2D_subopt_all(:,[2:5 end-1 6:end-2 1 end]));
        methods_name = {'AIC'  , 'HQ'  , 'MDL'  , 'AICc'  ,'AICc 19' , 'KICvc'  , 'WIC'  , 'NT'  , 'Q'};
      else
        imagesc([Nest_2D_subopt_all(:,2:end-1),Nest_2D_subopt_all(:,1),Nest_2D_subopt_all(:,end)]);
        methods_name = {'AIC'  , 'HQ'  , 'MDL'  , 'AICc'  , 'KICvc'  , 'WIC'  , 'NT'  , 'Q'};
      end
      
      %             set(gca,'TickLabelInterpreter','Latex')
      set(gca,'XtickLabel',methods_name)
      ylabel('Range bins','interpreter','Latex')
      %title('comaparision for a range line')
      cbr = colorbar;
      set(cbr,'YTick',0:1:M)
      title([  num2str(SNR_testing(SNR_idx))'' ' dB   Suboptimal'],'interpreter','Latex'  )
      
      if AICc_both==1
        figure(SNR_idx+110);clf;
        imagesc(Nest_2D_opt_all(:,[2:5 end-1 6:end-2 1 end]));
        methods_name = {'AIC'  , 'HQ'  , 'MDL'  , 'AICc'  ,'AICc 19' , 'KICvc'  , 'WIC'  , 'NT'  , 'Q'};
        
      else
        figure(SNR_idx+110);clf;
        imagesc([Nest_2D_opt_all(:,2:end-1),Nest_2D_opt_all(:,1),Nest_2D_opt_all(:,end)]);
        methods_name = {'AIC'  , 'HQ'  , 'MDL'  , 'AICc'  , 'KICvc'  , 'WIC'  , 'NT'  , 'Q'};
      end
      
      %             set(gca,'TickLabelInterpreter','Latex')
      set(gca,'XtickLabel',methods_name)
      ylabel('Range bins','interpreter','Latex')
      %title('comaparision for a range line')
      cbr = colorbar;
      set(cbr,'YTick',0:1:M)
      title([  num2str(SNR_testing(SNR_idx))'' ' dB   Optimal'],'interpreter','Latex'  )
    end
    
    % Plot log-likelihoods from the testing data
    % ------------------------------------------
    % NOTE (MOHANAD): Sravya had these plots in here script, but I
    % really don't see any benefit/reason in these plots. We already
    % had log-likelihood plots in the training phase. I left ir as is.
    %         if likelihood_plots == 1
    %             % TESTING DATA
    %
    %             Color         = {'b','r','k','c','m','g',[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560]};
    %             Marker_opt    = {'-*','-o','-+','-.','-s','-^','-<','->'};
    %             Marker_subopt = {'--*','--o','--+','--.','--s','--^','--<','-->'};
    %
    %             % PLOT LOG-LIKELIHOOD (NOT COST) WITH NORMALIZATION:
    %             % TESTING
    %             %--------------------------------------------------
    %             figure(6);
    %             subplot(length(SNR_testing),2,2*SNR_idx-1)
    %             for idx = 0:M
    %                 plot(0:M,mean(param_debug_testing.opt{SNR_idx,1}(find(param_debug_testing.sources_true{SNR_idx,1}==idx),:),1),'-*' )
    %                 hold on
    %             end
    %             grid on
    %             hold off
    
    %             for idx = 0:M
    %                 plot(0:M,mean(param_debug_testing.subopt{SNR_idx,1}(find(param_debug_testing.sources_true{SNR_idx,1}==idx),:),1),'--*' )
    %                 hold on,
    %             end
    
    %             set(gca,'TickLabelInterpreter','Latex')
    %             if SNR_idx ==1
    %                 title([  num2str(SNR_testing(SNR_idx))'' ' dB      --- Optimal, -- -- Suboptimal     Testing'],'interpreter','Latex'  )
    %             elseif SNR_idx == length(SNR_testing)
    %                 xlabel('k','interpreter','none')
    %                 ylabel('-2L','interpreter','none')
    %                 set(gca,'TickLabelInterpreter','Latex')
    %                 h_legend = legend('q = 0','q = 1','q = 2','q = 3','q = 4','q = 5','q = 6','Location','NorthEast');
    %                 set(h_legend,'Interpreter','Latex')
    %                 title([  num2str(SNR_testing(SNR_idx))'' ' dB '],'interpreter','Latex'  )
    %
    %             else
    %                 title([  num2str(SNR_testing(SNR_idx))'' ' dB'],'interpreter','Latex'  )
    %             end
    %
    %             grid on
    %             % TESTING PLOTS WITHOUT NORMALIZATION
    %             figure(34)
    %             subplot(length(SNR_testing),2,2*SNR_idx)
    %             for idx = 0:M
    
    %                 if param.norm_allign_zero ==1
    %                     plot(0:M,mean(param_debug_testing.opt{SNR_idx,1}(find(param_debug_testing.sources_true{SNR_idx,1}==idx),:)-repmat(norm_term_optimal,length(find(param_debug_testing.sources_true{SNR_idx,1}==idx)),1),1),'-*' )
    %                 else
    %                     plot(0:M,mean(param_debug_testing.opt{SNR_idx,1}(find(param_debug_testing.sources_true{SNR_idx,1}==idx),:)-repmat(opt_norm_term,length(find(param_debug_testing.sources_true{SNR_idx,1}==idx)),1),1),'-*' )
    %                 end
    %
    %                 hold on,
    %             end
    %             hold on
    %             for idx = 0:M
    %
    %                 if param.norm_allign_zero ==1
    %                     plot(0:M,mean(param_debug_testing.subopt{SNR_idx,1}(find(param_debug_testing.sources_true{SNR_idx,1}==idx),:)-repmat(norm_term_suboptimal,length(find(param_debug_testing.sources_true{SNR_idx,1}==idx)),1),1),'--*' )
    %                 else
    %                     plot(0:M,mean(param_debug_testing.subopt{SNR_idx,1}(find(param_debug_testing.sources_true{SNR_idx,1}==idx),:),1),'--*' )
    %                 end
    
    %                 hold on,
    %             end
    %
    %             if SNR_idx ==1
    %                 title([  num2str(SNR_testing(SNR_idx))'' ' dB    Without Normalization'],'interpreter','Latex'  )
    %             else
    %                 title([  num2str(SNR_testing(SNR_idx))'' ' dB'],'interpreter','Latex'  )
    %             end
    %             grid on
    %         end
    
    % AVERAGE CASE
    AVG_subopt = AVG_subopt + percentage_method_subopt_all{SNR_idx};
    AVG_opt = AVG_opt + percentage_method_opt_all{SNR_idx};
  end
  
  % ----------------------------------------------------------
  %           AVERAGE CASE
  % ----------------------------------------------------------
  AVG_subopt = AVG_subopt./length(SNR_testing);
  AVG_opt = AVG_opt./length(SNR_testing);
  
  if stat_test ==1
    % Plot percentage correct for sub suboptimal methods (average case)
    figure(4); subplot(2,2,4)
    
    plot(0:M , AVG_subopt(:,2),'m'),hold on
    plot(0:M , AVG_subopt(:,3),'-.r')
    plot(0:M , AVG_subopt(:,4),'-.g')
    plot(0:M , AVG_subopt(:,5),'--b')
    if   AICc_both ==1
      plot(0:M , AVG_subopt(:,8),'-ob')
    end
    plot(0:M ,AVG_subopt(:,6),'--c.')
    plot(0:M , AVG_subopt(:,7),'-*k')
    
    plot(0:M , AVG_subopt(:,1),'-*g')
    
    xlim([0 M])
    ylim([0 100])
    grid on
    
    if SNR_idx == length(SNR_testing)
      xlabel('Number of sources, q','interpreter','none')
      
      %             set(gca,'TickLabelInterpreter','Latex')
      %         h_legend = legend('AIC', 'HQ', 'MDL', 'AICc', 'KICvc', 'WIC','NT','Location','NorthEast');
      %         set(h_legend,'Interpreter','Latex')
      %
    end
    title('ALL','interpreter','Latex'  )
    
    % Plot percentage correct for sub optimal methods (average case)
    figure(5),subplot(2,2,4)
    
    plot(0:M , AVG_opt(:,2),'m'), hold on
    plot(0:M , AVG_opt(:,3),'-.r')
    plot(0:M , AVG_opt(:,4),'-.g')
    plot(0:M , AVG_opt(:,5),'--b')
    if   AICc_both ==1
      plot(0:M , AVG_opt(:,8),'-ob')
    end
    plot(0:M ,AVG_opt(:,6),'--c.')
    plot(0:M , AVG_opt(:,7),'-*k')
    plot(0:M , AVG_opt(:,1),'-*g')
    
    xlim([0 M])
    ylim([0 100])
    grid on
    if SNR_idx == length(SNR_testing)
      xlabel('Number of sources, q','interpreter','none')
      
      %             set(gca,'TickLabelInterpreter','Latex')
      %         h_legend = legend('AIC', 'HQ', 'MDL', 'AICc', 'KICvc', 'WIC','NT','Location','NorthEast');
      %         set(h_legend,'Interpreter','Latex')
      
      
      title('ALL','interpreter','Latex'  )
    end
    
    figure(4);
    subplot(2,2,1)
    ylabel('% Correct','interpreter','none')
    h = suptitle('Suboptimal MOE');
    set(h,'FontSize',11,'FontWeight','normal');
    
    figure(5);
    subplot(2,2,1)
    ylabel('% Correct','interpreter','none')
    h =  suptitle('Optimal MOE');
    set(h,'FontSize',11,'FontWeight','normal')
  end
  toc
  %% To plot Slice model
  
  % slice = 11;
  % surface_z = results.z_grid-results.param.monte.target_param{1}.z.mean;
  % surface_z_groung_truth =results.surf_model.z-results.param.monte.target_param{1}.z.mean;
  % figure(5); clf;  plot(results.surf_model.y, surface_z_groung_truth(:,slice),'b');
  % hold on
  % plot(results.surf_model.y, surface_z(:,slice),'r');
  % xlim([-1500 1500])
  % ylim([-200 250])
  % title('Slice - surface model');
  % xlabel('Cross-track (m)');
  % ylabel('WGS84-Elevation (m)');
  % hold off
  % legend('Ground-truth surface','Actual surface');
  % grid on
  %
  % % Slice - Range Bin v/s DOA
  % figure(6),clf
  % scatter(results.tomo.doa(:,1,slice)*(180/pi),results.array_param.bins, 20 , 10*log10(results.tomo.power(:,1,slice)),'fill');
  % colorbar
  % hold on
  % scatter(results.tomo.doa(:,2,slice)*(180/pi),results.array_param.bins, 20 , 10*log10(results.tomo.power(:,2,slice)),'fill');
  % colorbar
  % set(gca,'Ydir','reverse')
  % xlim([-60 60])
  % ylim([1 100])
  % title('Slice');
  % xlabel('DOA (deg)');
  % ylabel('Range bin');
  % cc = caxis;
  % h_cb = colorbar;
  % set(get(h_cb,'YLabel'),'String','Relative power (dB)');
  % grid on
  %
  
  if param.debug_level >= 3
    return
  end
  
end
