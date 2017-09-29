% script Sim.run_doa
%
% Description: Run script for Sim.doa (Calls Sim.doa_cluster).
clear('param_override');
param_override = [];
param = [];
physical_constants;
%Multiple tasks (i.e. jobs)
param_override.sched.type = 'custom_torque';

%Single task (i.e. job)
% param_override.sched.type = 'no scheduler';

param_override.sched.cluster_size = inf;
param_override.sched.rerun_only   = false;
param_override.sched.stop_on_fail = false;

param_override.sched.num_submissions = 170;
param_override.sched.group_size      = 1;
param_override.sched.group_walltime  = 2*86400;
param_override.debug_level     = 3;
param_override.sched.submit_arguments = '-l nodes=1:ppn=1,pmem=14000mb,walltime=48:00:00'; % 1 process
% param_override.sched.submit_arguments = '-l nodes=1:ppn=1,pmem=8000mb,walltime=120:00';

%% Automated loading section
% =========================================================================
global gRadar;

if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%% Radar and source parameters
% =========================================================================
if 1
  fprintf('=======================================================================\n');
  fprintf('Running sim.run_doa example #2\n');
  fprintf('  Wideband radar depth sounder simulation\n');
  fprintf('  Linear array in y-dimension\n');
  fprintf('=======================================================================\n');
  %% Setup simulation parameters
  param = [];
  
  % Debug level of 3 causes this function to stop early and output
  % simulated data, simulation parameters, and array processing results
  % in a "results" structure.
  param.debug_level = 3;
  
  %% Source parameters
  fc = 335e6;
  BW = 370e6;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.t0                      = 2*(500+300-20)/c;
  param.src.t1                      = 2*(800+300+20)/c;
  param.src.ft_func                 = @(t) tukeywin_cont(t * BW);
  param.src.ft_wind                 = @(N) hanning(N);
  
  % TO DO: Use lever_arm instead of lever_arm_example
  %=====================================================
  % lever_arm.m function needs to be called for each channel.
  % lever_arm_exampl.m can return the whole array at once
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  
  % Nsep: number of lambda/4 steps between sensors
  Nsep = 1;
  %TO DO: Use all subarrays, that is Nc=24
  %--------------------------------------
  % Nc: number of sensors
  Nc = 8;
  % Arguments for a linear array in y-dimension
  % param.src.lever_arm.fh_args: Arguments for a linear array in y-dimension
  param.src.lever_arm.fh_args       = {[], 1, 1:Nc, [0; Nsep*c/fc/2; 0]};
  param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
  
  % DOA method parameters
  param.method.list                   = [7];
  
  %% Simulation Runs Setup
  
  % Cross track monte carlo setup
  param.monte.target_func = @sim.surface_gen;
  param.monte.runs = 1;
  param.monte.random_seed_offset = 0;
  
  % Target surface parameters
  surf_param = [];
  surf_param.z.mean = -(500+300);
  surf_param.z.rms_height = 10;
  surf_param.z.corr_length_x = 600;
  surf_param.z.corr_length_y = 600;
  surf_param.rcs_in.mean = 0;
  surf_param.rcs_in.var = 100;
  surf_param.dy = 10;
  surf_param.y_range = [-2000 2000];
  surf_param.dx = 10;
  surf_param.x_range = [-2000 2000];
  surf_param.x = [-500:40:500];
  surf_param.y = [-400:10:400].';
  param.monte.target_param{1} = surf_param;
  
  %% Array Processing parameters
  array_param = [];
  fs = param.src.f1-param.src.f0;
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
  %   array_param.Nsv = {'theta', asin(ky/k)};
  array_param.Nsv = {'theta', fftshift(theta)}; %Never use theta from - to +.
  
  array_param.sv_fh = @array_proc_sv;
  
  array_param.dbin = 1;
  array_param.dline = 1;
  
  %   array_param.bin_rng = -1:1;
  array_param.rline_rng = -10:10;
  
  array_param.Nsig = 2;
  
  array_param.init = 'ap';
  array_param.theta_guard = (max(theta)-min(theta))/(4*Nc);
  
  array_param.W  = 1;
  array_param.NB = 1; %number of bands
  
  %Note that length(array_param.bin_rng) MUST be >= array_param.NB,
  %otherwise, DCM will be all zeros, and will get an error.
  if param.method.list == 9 && array_param.NB>1
    if mod(array_param.NB,2)==0
      array_param.NB = array_param.NB+1;
      warning('Number of sub-bands (array_param.NB) you entered is even, and it should be odd number, so automatically increased by 1')
    end
    array_param.bin_rng = -floor((array_param.NB-1)/2):floor(array_param.NB/2);
  else
    array_param.bin_rng = -1:1;
  end
  
  if length(array_param.bin_rng) < array_param.NB
    error('length(array_param.bin_rng) MUST be >= doa_param.nb_filter_banks (i.e. array_param.NB)')
  end
  
  
  dt = 1/fs;
  array_param.imp_resp.time_vec = -3*array_param.W*dt : dt/8 : 3*array_param.W*dt;
  BW = abs(param.src.f1 - param.src.f0);
  array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / BW);
  
  for idx = 1:array_param.Nsig
    array_param.doa_constraints(idx).method = 'fixed';
    array_param.doa_constraints(idx).init_src_limits = [min(theta) max(theta)]*180/pi;
    array_param.doa_constraints(idx).src_limits = [min(theta) max(theta)]*180/pi;
  end
  
  param.array_param = array_param;
  clear array_param;
  
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
  
  
  %% Run the simulation
  doa_cluster(param,param_override);

% task_param = param;
% fh  = @crosstrack;
% arg = {task_param};
% results=fh(arg{:});

if param.debug_level >= 3
  return
end

  % TO DO: Find a way to store the output of crosstrack.m (results) to use
  % it in sim.crosstrack_rmse(param,results)
  %========================================================================
  if 0
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
    
  end
end
