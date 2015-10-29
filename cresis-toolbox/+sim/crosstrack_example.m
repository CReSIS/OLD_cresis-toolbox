% script sim.crosstrack_example
%
% Example setup scripts for the sim.crosstrack function. To run, enable one
% of the if {0|1} blocks below and then run the script.
%
% Author: Sean Holloway, John Paden, Theresa Stumpf
%
% See also: sim.crosstrack.m, sim.crosstrack_data.m, sim.crosstrack_example.m

physical_constants;

if 1
  fprintf('=======================================================================\n');
  fprintf('Running sim.crosstrack_example example #1\n');
  fprintf('  Narrowband altimeter (e.g. snow radar) simulation \n');
  fprintf('=======================================================================\n');
  %% Setup simulation parameters
  param = [];
  
  % Debug level of 3 causes this function to stop early and output
  % simulated data, simulation parameters, and array processing results
  % in a "results" structure.
  param.debug_level = 3;
  
  %% Source parameters
  fc = 360e9;
  BW = 6000e6;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.t0                      = 2*(500-5)/c;
  param.src.t1                      = 2*(500+5)/c;
  param.src.ft_func                 = @(t) tukeywin_cont(t * BW);
  param.src.ft_wind                 = @(N) hanning(N);
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  Nc = 16;
  param.src.lever_arm.args          = {[],ones(1,Nc),1:Nc,Nc,c/fc/4};
  param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
  
  % DOA method parameters
  param.method.list                   = [7];
  param.method.Nsv                    = 32;
  param.method.OneD_Nsv               = 256;
  param.method.theta_guard            = 0.75/180*pi;
  
  %% Simulation Runs Setup
  
  % Cross track monte carlo setup
  param.monte.target_func = @sim.surface_gen;
  param.monte.runs = 1;
  param.monte.random_seed_offset = 0;
  
  % Target surface parameters
  surf_param = [];
  surf_param.z.mean = -500;
  surf_param.z.rms_height = 0.15;
  surf_param.z.corr_length_x = 40;
  surf_param.z.corr_length_y = 40;
  surf_param.rcs_in.mean = 0;
  surf_param.rcs_in.var = 100;
  surf_param.dy = 1;
  surf_param.y_range = [-1500 1500];
  surf_param.dx = 1;
  surf_param.x_range = [-1500 1500];
  surf_param.x = [-50:1:50];
  surf_param.y = [-50:1:50].';
  param.monte.target_param{1} = surf_param;
  
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
  
  %% Array Processing parameters
  array_param = [];
  
  fc = (param.src.f0+param.src.f1)/2;
  fs = param.src.f1-param.src.f0;
  
  Nsep = 16;
  NN = Nc*Nsep;
  lambda = fc/3e8;
  k = 2*pi/lambda;
  dy = lambda/4;
  My = 16;
  dky = 2*pi / (NN*dy) / My;
  
  ky = dky * ifftshift(-floor(My*Nc/2) : floor((My*Nc-1)/2));
  theta = fftshift(asin(ky/k));
  array_param.Nsv = {'theta', asin(ky/k)};
  
  array_param.sv_fh = @array_proc_sv;
  
  array_param.dbin = 1;
  array_param.dline = 1;
  
  array_param.bin_rng = -1:1;
  array_param.rline_rng = -10:10;
  
  array_param.Nsig = 2;
  
  array_param.init = 'ap';
  array_param.theta_guard = (max(theta)-min(theta))/64;
  
  array_param.W = 1;
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
  results = sim.crosstrack(param);
  
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
  
  return;
end

if 0
  fprintf('=======================================================================\n');
  fprintf('Running sim.crosstrack_example example #2\n');
  fprintf('  Narrowband radar depth sounder simulation\n');
  fprintf('=======================================================================\n');
  %% Setup simulation parameters
  param = [];
  
  % Debug level of 3 causes this function to stop early and output
  % simulated data, simulation parameters, and array processing results
  % in a "results" structure.
  param.debug_level = 3;
  
  %% Source parameters
  fc = 195e6;
  BW = 10e6;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.t0                      = 2*(1500-50)/c;
  param.src.t1                      = 2*(1500+500)/c;
  param.src.ft_func                 = @(t) tukeywin_cont(t * BW);
  param.src.ft_wind                 = @(N) hanning(N);
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  Nc = 8;
  param.src.lever_arm.args          = {[],ones(1,Nc),1:Nc,Nc,c/fc/4};
  param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
  
  % DOA method parameters
  param.method.list                   = [7];
  param.method.Nsv                    = 9;
  param.method.OneD_Nsv               = 256;
  param.method.theta_guard            = 0.75/180*pi;
  
  %% Simulation Runs Setup
  
  % Cross track monte carlo setup
  param.monte.target_func = @sim.surface_gen;
  param.monte.runs = 1;
  param.monte.random_seed_offset = 0;
  
  % Target surface parameters
  surf_param = [];
  surf_param.z.mean = -500;
  surf_param.z.rms_height = 0.15;
  surf_param.z.corr_length_x = 40;
  surf_param.z.corr_length_y = 40;
  surf_param.rcs_in.mean = 0;
  surf_param.rcs_in.var = 100;
  surf_param.dy = 1;
  surf_param.y_range = [-1500 1500];
  surf_param.dx = 1;
  surf_param.x_range = [-1500 1500];
  surf_param.x = [-50:1:50];
  surf_param.y = [-50:1:50].';
  param.monte.target_param{1} = surf_param;
  
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
  
  %% Array Processing parameters
  array_param = [];
  
  fc = (param.src.f0+param.src.f1)/2;
  fs = param.src.f1-param.src.f0;
  
  Nsep = 16;
  NN = Nc*Nsep;
  lambda = fc/3e8;
  k = 2*pi/lambda;
  dy = lambda/4;
  My = 16;
  dky = 2*pi / (NN*dy) / My;
  
  ky = dky * ifftshift(-floor(My*Nc/2) : floor((My*Nc-1)/2));
  theta = fftshift(asin(ky/k));
  array_param.Nsv = {'theta', asin(ky/k)};
  
  array_param.sv_fh = @array_proc_sv;
  
  array_param.dbin = 1;
  array_param.dline = 1;
  
  array_param.bin_rng = -1:1;
  array_param.rline_rng = -10:10;
  
  array_param.Nsig = 2;
  
  array_param.init = 'ap';
  array_param.theta_guard = (max(theta)-min(theta))/64;
  
  array_param.W = 1;
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
  results = sim.crosstrack(param);
  
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
  
  return;
end
