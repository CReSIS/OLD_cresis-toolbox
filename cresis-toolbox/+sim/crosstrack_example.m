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
  %% Mohanad: 2D tutorial
  fprintf('=======================================================================\n');
  fprintf('Running sim.crosstrack_example example #2\n');
  fprintf('  Narrowband radar depth sounder simulation\n');
  fprintf('  Linear array in y-dimension\n');
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
  param.src.t0                      = 2*(1500-500)/c;
  param.src.t1                      = 2*(1500+1000)/c;
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
  param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
  
  % DOA method parameters
  array_proc_methods;
  param.method.list                   = [MLE_METHOD];
%   param.method.method_mode = 'estimator';
  
  %% Simulation Runs Setup
  
  % Cross track monte carlo setup
  param.monte.target_func = @sim.surface_gen;
  param.monte.runs = 1;
  param.monte.random_seed_offset = 0;
  
  % Target surface parameters
  surf_param = [];
  surf_param.z.mean = -1500;
  surf_param.z.rms_height = 0;
  surf_param.z.corr_length_x = 400;
  surf_param.z.corr_length_y = 400;
  surf_param.rcs_in.mean = 0;
  surf_param.rcs_in.var = 100;
  surf_param.dy = 10;
  surf_param.y_range = [-2500 2500];
  surf_param.dx = 10;
  surf_param.x_range = [-2500 2500];
%   surf_param.x = [-1100:-1000];
  surf_param.x = [-1020:-1000];
%   surf_param.y = [-1500:20:1500].';
  surf_param.y = abs(surf_param.z.mean)*tan([-10  10]'*pi/180);
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
  array_param.theta = asin(ky/k);
  
  array_param.sv_fh = @array_proc_sv;
  
  array_param.dbin = 1;
  array_param.dline = 1;
  
  array_param.bin_rng = 0;
  array_param.line_rng = -10:10;
  
  array_param.Nsrc = 2;
  
  array_param.init = 'ap';
  array_param.doa_theta_guard = 1.5 * pi/180;
%   array_param.doa_theta_guard = (max(theta)-min(theta))/(4*Nc);
  
  array_param.Nsubband = 1;
  dt = 1/fs;
  array_param.imp_resp.time_vec = -3*array_param.Nsubband*dt : dt/8 : 3*array_param.Nsubband*dt;
  BW = abs(param.src.f1 - param.src.f0);
  array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / BW);
  
  for idx = 1:array_param.Nsrc
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
  
  %% To plot Slice model
  
  slice = 1;
  surface_z = results.z_grid-results.param.monte.target_param{1}.z.mean;
  surface_z_groung_truth =results.surf_model.z-results.param.monte.target_param{1}.z.mean;
  figure(9); clf;  plot(results.surf_model.y, surface_z_groung_truth(:,slice),'b');
%   hold on
%   plot(results.surf_model.y, surface_z(:,slice),'r');
  xlim([-1500 1500])
  ylim([-200 250])
  title('Slice - surface model');
  xlabel('Cross-track (m)');
  ylabel('WGS84-Elevation (m)');
  hold off
  legend('Ground-truth surface','Actual surface');
  grid on
  
  % Slice - Range Bin v/s DOA
  figure(11),clf
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
  fprintf('Running sim.crosstrack_example example #1\n');
  fprintf('  Narrowband altimeter (e.g. snow radar) simulation \n');
  fprintf('  Linear array in y-dimension\n');
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
  % Nsep: number of lambda/4 steps between sensors
  Nsep = 8;
  % Nc: number of sensors
  Nc = 32;
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
  array_param.theta = asin(ky/k)*180/pi;
  
  array_param.sv_fh = @array_proc_sv;
  
  array_param.dbin = 1;
  array_param.dline = 1;
  
  array_param.bin_rng = -1:1;
  array_param.line_rng = -10:10;
  
  array_param.Nsrc = 2;
  
  array_param.init = 'ap';
  array_param.doa_theta_guard = (max(theta)-min(theta))/length(ky);
  
  array_param.Nsubband = 1;
  dt = 1/fs;
  array_param.imp_resp.time_vec = -3*array_param.Nsubband*dt : dt/8 : 3*array_param.Nsubband*dt;
  BW = abs(param.src.f1 - param.src.f0);
  array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / BW);
  
  for idx = 1:array_param.Nsrc
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
  fprintf('  Linear array in y-dimension\n');
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
  param.src.t0                      = 2*(1500-500)/c;
  param.src.t1                      = 2*(1500+1000)/c;
  param.src.ft_func                 = @(t) tukeywin_cont(t * BW);
  param.src.ft_wind                 = @(N) hanning(N);
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  % Nsep: number of lambda/4 steps between sensors
  Nsep = 1;
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
  surf_param.z.mean = -1500;
  surf_param.z.rms_height = 100;
  surf_param.z.corr_length_x = 400;
  surf_param.z.corr_length_y = 400;
  surf_param.rcs_in.mean = 0;
  surf_param.rcs_in.var = 100;
  surf_param.dy = 10;
  surf_param.y_range = [-2500 2500];
  surf_param.dx = 10;
  surf_param.x_range = [-2500 2500];
  surf_param.x = [-500:5:500];
  surf_param.y = [-1500:20:1500].';
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
  array_param.line_rng = -10:10;
  
  array_param.Nsrc = 2;
  
  array_param.init = 'ap';
  array_param.doa_theta_guard = (max(theta)-min(theta))/(4*Nc);
  
  array_param.Nsubband = 1;
  dt = 1/fs;
  array_param.imp_resp.time_vec = -3*array_param.Nsubband*dt : dt/8 : 3*array_param.Nsubband*dt;
  BW = abs(param.src.f1 - param.src.f0);
  array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / BW);
  
  for idx = 1:array_param.Nsrc
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
  %% Multipath and Mutual coupling effects: New implementation
  % ----------------------------------------------------------------------
  % Here is how we simulate the effect of the MPCs created by a target 
  % at some range-bin on a target at another range-bin:
  % 1) Set the location of the main target that is going to create the MPCs.
  % 2) Create the reflectors locations by defining DOA_MP and z_MP then
  % calculating y_MP.
  % The MPC may not affect the same range-bin that created this MPC, unless
  % the MPC arrives at the same time (within the range-resolution) as the
  % main path.

    %% Setup simulation parameters
  param = [];
  
  % Debug level of 3 causes this function to stop early and output
  % simulated data, simulation parameters, and array processing results
  % in a "results" structure.
  param.debug_level = 3;
  
  %% Source parameters
  fc = 195e6;
  BW = 30e6;
  % Nc: number of sensors
  Nc = 7;
  flight_height  = 0;
  % Set the DoAs and weights of the main paths and MPCs
    % Orthogonal DoAs (to test if it is actually working)
    lambda = c/fc;
    d_y    = lambda/2;
    L = Nc*d_y;
    k_spacing = [-(Nc-1)/2:1:(Nc-1)/2]*(lambda/L);
    DOA_orthogonal = asin(k_spacing);
%     DOA = DOA_orthogonal(5).';
%     DOA_MP = DOA_orthogonal([1 7]).';
    DOA = [30]*pi/180;
    DOA_MP =[-5 8].'*pi/180;
    
    % Setting the MPC weight to -inf ignores the MPC (controls number of
    % MPCs).
    MP_weights = [-inf -inf]; % in dB
    
    % Main target location
    z_pos = -1500;
    Range = abs(z_pos)/cos(DOA); % Range to the main target
    y_pos = Range*sin(DOA);      % y-axis of the main target
    
    % MPCs location
    z_pos_MP = [z_pos+400, z_pos+480];
    y_pos_MP = z_pos_MP.*tan(DOA_MP.');
    Range_MP = abs(z_pos_MP)./cos(DOA_MP.');

%     Range_MP = Range*ones(1,length(DOA_MP));
%     z_pos_MP = -Range_MP.*cos(DOA_MP.');
%     y_pos_MP = Range_MP.*sin(DOA_MP.');
%     R_main_target = Range*ones(1,length(DOA_MP));
    
    target_to_reflector_dist   = sqrt((z_pos-z_pos_MP).^2 + (y_pos_MP-y_pos).^2);
    reflector_to_receiver_dist = sqrt(z_pos_MP.^2 + y_pos_MP.^2);
    total_distance_MP = Range + target_to_reflector_dist + reflector_to_receiver_dist;
    tau_MP = (total_distance_MP - 2*Range)./c; % Relative MPCs delay
    
    % Place a target at the same range as the MPCs (this target doesn't
    % creat a MPC). Make sure that this target and the one created the MPCs
    % are not within one range-resolution distance from each other.
%     Range_2 = total_distance_MP(1)/2;
%     y_pos_2 = sqrt(Range_2^2 - z_pos^2);
%     DOA_2 = sign(y_pos_2)*atan(abs(y_pos_2/z_pos));
    
%     DOA = [DOA DOA_2].';
    DOA = DOA;
  
  if 0
    % Debug: plot problem geometry
    figure(100);clf
    hold on
    h1 = plot(y_pos,z_pos,'*b','LineWidth',10);             % Main target           
    h2 = plot(y_pos_MP(1),z_pos_MP(1),'*r','LineWidth',10); % Reflector 1
    h3 = plot(y_pos_MP(2),z_pos_MP(2),'*m','LineWidth',10); % Reflector 2
    h4 = plot(0,0,'sk','LineWidth',10);                     % Radar
    
    plot([y_pos,y_pos_MP(1)],[z_pos,z_pos_MP(1)],'r','LineWidth',1.5) % Target to reflector 1
    plot([y_pos,y_pos_MP(2)],[z_pos,z_pos_MP(2)],'m','LineWidth',1.5) % Target to reflector 2
    plot([y_pos 0],[z_pos 0],'b','LineWidth',1.5)                     % Radar to target
    plot([0 y_pos_MP(1)],[0 z_pos_MP(1)],'r','LineWidth',1.5)         % Radar to reflector 1
    plot([0 y_pos_MP(2)],[0 z_pos_MP(2)],'m','LineWidth',1.5)         % Radar to reflector 2
    plot([0 0],[0 z_pos],'-.k')                                       % nadir
    
    grid on
    xlabel('y-axis')
    ylabel('z-axis')
    title('Multipath problem geometry')
    legend([h1 h2 h3 h4],'Main target','Reflector 1','Reflector 2','Radar','Location','southwest')
  end
    
  % MP parameters associated with each target
    MP_params{1}.y_pos_MP = y_pos_MP;
    MP_params{1}.z_pos_MP = z_pos_MP;
    MP_params{1}.w_MP     = MP_weights;
    
    MP_params{2}.y_pos_MP = y_pos_MP;
    MP_params{2}.z_pos_MP = z_pos_MP;
    MP_params{2}.w_MP     = [-inf -inf]; % Second target doesn't create MPCs
    
  param.MP_params = MP_params;
  
  % Set the range gate based on the flight height and targets DoA
  max_range      = max(Range + target_to_reflector_dist + reflector_to_receiver_dist);
  min_range      = min(Range);
  
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.t0                      = 2*(min_range-500)/c;
  param.src.t1                      = 2*(max_range+1000)/c;
  if param.src.t0 < 0
    param.src.t0 = 0;
  end
  param.src.ft_func                 = @(t) tukeywin_cont(t * BW);
  param.src.ft_wind                 = @(N) hanning(N);
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  
  % Nsep: number of lambda/4 steps between sensors
  Nsep = 1;
  % Arguments for a linear array in y-dimension
  % param.src.lever_arm.fh_args: Arguments for a linear array in y-dimension
  param.src.lever_arm.fh_args       = {[], 1, 1:Nc, [0; Nsep*c/fc/2; 0]};
  phase_center = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
  
  % Change reference sensor if you want to
  ref_chan = 1;
  phase_center(2,:) = phase_center(2,:) - phase_center(2,ref_chan);
  phase_center(3,:) = phase_center(3,:) - phase_center(3,ref_chan);
 
  % Account for roll angle. Default is 0
  theta_roll = 0 * pi/180;
  phase_center(2,:) = phase_center(2,:) * cos(theta_roll);
  for chan_idx = 1:Nc
  d_chan_ref      = sqrt(abs(diff(phase_center(2,[ref_chan  chan_idx])))^2 + abs(diff(phase_center(3,[ref_chan  chan_idx])))^2);
  phase_center(3,chan_idx) = phase_center(3,chan_idx) + d_chan_ref * sin(theta_roll);
  end
  param.src.phase_center = phase_center;
  ypc = -phase_center(2,:).';
  zpc = -phase_center(3,:).';
  
  dt = 1/BW;
  % Nt: number of range bins (fast time samples)
  Nt = floor((param.src.t1-param.src.t0)/dt);
  
  % Time = fast time axis (sec)
  Time = param.src.t0 + dt*(0:Nt-1).';
  %   param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
  param.src.noise_power             = zeros(1,Nc);
  
  % DOA method parameters
  param.method.list                   = [7];
  
  %% Simulation Runs Setup
  
  % Cross track monte carlo setup
  param.monte.target_func = @sim.surface_gen;
  param.monte.runs = 1;
  param.monte.random_seed_offset = 0;
  
  % Target surface parameters
  surf_param = [];
  surf_param.z.mean = z_pos;
  surf_param.z.rms_height = 0;
  surf_param.z.corr_length_x = 400;
  surf_param.z.corr_length_y = 400;
  surf_param.rcs_in.mean = 0;
  surf_param.rcs_in.var = 1e2;
  surf_param.dy = 10;
  %   surf_param.y_range = [-2500 2500];
  surf_param.dx = 10;
  surf_param.x_range = [-2500 2500];
  surf_param.x = [-1500:-1500+100]; % 199...100 range-lines
  %   surf_param.y = [-1500:20:1500].';
%   surf_param.y = [y_pos y_pos_2].';
  surf_param.y = y_pos;
  
  % y_range should >= maximum y (i.e. targets shoulb be inside the imaged
  % swath)
  y_range = max(target_to_reflector_dist + reflector_to_receiver_dist+Range);
  if y_range == 0
    surf_param.y_range = [-1000 1000];
  else
    surf_param.y_range = [-1.5*y_range 1.5*y_range];
  end
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
  theta = [theta,pi/2];
  array_param.Nsv = {'theta', asin(ky/k)};
  
  array_param.sv_fh = @array_proc_sv;
  
  array_param.dbin = 1;
  array_param.dline = 1;
  
  array_param.bin_rng = 0;
  array_param.line_rng = -50:50;
  
  array_param.Nsrc = 2;
  Nsrc = array_param.Nsrc;
  
  array_param.init = 'ap';
  array_param.doa_theta_guard = 2*pi/180; %(max(theta)-min(theta))/(4*Nc);
  
  array_param.Nsubband = 1;
  dt = 1/fs;
  array_param.imp_resp.time_vec = -3*array_param.Nsubband*dt : dt/8 : 3*array_param.Nsubband*dt;
  BW = abs(param.src.f1 - param.src.f0);
  array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / BW);
  
  for idx = 1:array_param.Nsrc
    array_param.doa_constraints(idx).method = 'fixed';
    array_param.doa_constraints(idx).init_src_limits = [min(theta) max(theta)]*180/pi;
    array_param.doa_constraints(idx).src_limits = [min(theta) max(theta)]*180/pi;
  end
  
  % -----------------------------------------------------------------------
  % Mutual coupling effect
  % -----------------------------------------------------------------------
  % Define S-parameters matrix (symmetric and Toepletz) assuming matched
  % dipoles
  if 0
    S = [0.000   0.200  -0.40  -0.2j   0.10   0.05j   0.010;...
      0.200   0.000   0.20  -0.40  -0.2j   0.100   0.05j;...
      -0.400   0.200   0.00   0.20  -0.40  -0.20j   0.100;...
      -0.20j  -0.400   0.20   0.00   0.20   0.400  -0.20j;...
      0.100  -0.20j  -0.40   0.20   0.00   0.200  -0.400;...
      0.05j   0.100   0.2j  -0.40   0.20   0.000   0.200;...
      0.010   0.05j   0.10  -0.2j  -0.40   0.200   0.000];
    % Define Mutual coupling matrix
    C = eye(Nc) - S;
    %   C = eye(Nc);
    param.src.mutual_coup_mtx = C;
  else
    param.src.mutual_coup_mtx = [];
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
  results = crosstrack(param);
  
  sim_data = squeeze(results.sim_data{1}); % Nt-by-Nx-by-Nc  
%   est_doa         = results.tomo.doa;
  actual_doa_cell = results.actual_doa;
  array_param     = results.array_param;
  
  % Convert actual_doa from cell array into a matrix (maximum Nsrc targets
  % per range-bin).
  actual_doa = NaN(length(array_param.bins),Nsrc,length(array_param.lines));
  for lineIdx = 1:length(array_param.lines)
    for binIdx = 1:length(array_param.bins)
      if ~isempty(actual_doa_cell{binIdx,lineIdx})
        doa_tmp = actual_doa_cell{binIdx,lineIdx};
        if length(doa_tmp)<Nsrc
          doa = NaN(Nsrc,1);
          doa(1:length(doa_tmp)) = doa_tmp;
        else
          doa = doa_tmp(1:Nsrc);
        end
        actual_doa(binIdx,1:Nsrc,lineIdx) = actual_doa_cell{binIdx,lineIdx}; % Nt*Nsrc*Nx
      end
    end
  end
  
  if 1
    % Calculate the range-bins to which the MPCs blong to. The range is
    % calculated as half of the total distance from the send to receive.
    range = Time*c/2;
    dr = c/(2*BW);
    MPC_rbins = [];
    for lineIdx = 1:length(array_param.lines)
      for MP_idx = 1:length(DOA_MP)
        %       rbins = find(abs(Range_MP(MP_idx)-range) < dr/2*2);
        [~,min_idx] = min(abs(total_distance_MP(MP_idx)/2 - range));
        MPC_rbins(MP_idx,lineIdx) = min_idx;
      end
    end
  end
  
  % Calculate the eigenvalues of each test range-bin
  eigenvalue_all  = [];
  sample_data_all = [];
  sample_data     = [];
  for lineIdx_idx = 1:length(array_param.lines)
    line_idx = array_param.lines(lineIdx_idx);
    good_rbins = sort(find(~isnan(actual_doa(:,2,lineIdx_idx))),'ascend');
%     first_rbin = MPC_rbins(1);
    rqd_rbin = good_rbins(1);%min(good_rbins);
    range_bin_idxs = rqd_rbin + [0];
    for binIdx_idx = 1:length(range_bin_idxs)
      bin_idx = range_bin_idxs(binIdx_idx);
      sample_data = squeeze(sim_data(bin_idx,line_idx+array_param.line_rng,:)).';
      
      if 0
        % FB averaging and spatial smoothing
        DCM = 1/size(sample_data,2) * sample_data*sample_data';
        Rxx1 = 1/(size(sample_data(1:end-1,:),2))*sample_data(1:end-1,:) * sample_data(1:end-1,:)';
        Rxx2 = 1/(size(sample_data(2:end,:),2))*sample_data(2:end,:) * sample_data(2:end,:)';
        
        % Apply FB averaging for each subarray
        reflect_mtx = flipud(eye(size(Rxx1)));
        Rxx1_fb = (1/2)*(Rxx1 + reflect_mtx * conj(Rxx1) * reflect_mtx);
        Rxx2_fb = (1/2)*(Rxx2 + reflect_mtx * conj(Rxx2) * reflect_mtx);
        
        % Average the two DCMs
        Rxx_ss = 1/2*(Rxx1_fb+Rxx2_fb); % (Nc-1)-by-(Nc-1)
        
        % Handle the lost sensor (or DOF) such that the final DCM is Nc-by-Nc
        Rxx_tmp = zeros(Nc);
        Rxx_tmp(1:end-1,1:end-1) = Rxx_ss;
        Rxx_tmp(end,:) = DCM(end,:);
        Rxx_tmp(:,end) = DCM(:,end);
        
        DCM = Rxx_tmp; % Nc-by-Nc matrix
      elseif 0
        % FB averaging only
        DCM = 1/size(sample_data,2) * sample_data*sample_data';
        reflect_mtx = flipud(eye(size(Rxx_calib)));
        DCM = (1/2)*(DCM + reflect_mtx * conj(DCM) * reflect_mtx);
      else
        DCM = 1/size(sample_data,2) * sample_data*sample_data';
      end
      if bin_idx == rqd_rbin
        DCM_tmp = DCM;
      end
      
      eigenvalue = real(sort(eig(DCM),'descend'));
      eigenvalue_all(:,end+1) = eigenvalue;
      
      sample_data_all  = [sample_data_all;sample_data(:)];
    end
  end
  mean_eigval = mean(real(eigenvalue_all),2);
  mean_eigval = mean_eigval./max(mean_eigval);
  
  if 1
    % Debug
    % Estimate the SNR
    SNR = 10*log10(mean(abs(sample_data_all).^2));

    TBP = tau_MP/(1/BW);
    fract_BW = BW/fc;
    
%     good_rbins_1 = find(~isnan(results.tomo.doa(:,1,lineIdx_idx)));
%     good_rbins_2 = find(~isnan(results.tomo.doa(:,2,lineIdx_idx)));
%     est_doa_1 = results.tomo.doa(good_rbins_1,1,lineIdx_idx)*180/pi;
%     est_doa_2 = results.tomo.doa(good_rbins_2,2,lineIdx_idx)*180/pi;
%     est_doa = sort(unique([est_doa_1' est_doa_2']),'ascend');
   
    if length(DOA_MP) == 2
    sprintf('\nTBP = [%2.4f  %2.4f]\n \nFractional BW = %2.4e\n \nSNR = %2.1f\n \n DoA = %2.2f deg.\n \n MPCs DoA = [%2.2f %2.2f] deg.\n',...
      TBP(1), TBP(2), fract_BW, SNR, DOA(1)*180/pi,DOA_MP(1)*180/pi,DOA_MP(2)*180/pi)

    elseif length(DOA_MP) == 1
      sprintf('\nTBP = %2.4f\n \nFractional BW = %2.4e\n \nSNR = %2.1f\n \n DoA = %2.2f\n \n MPCs DoA = %2.2f\n',...
      TBP(1), fract_BW, SNR, DOA*180/pi,DOA_MP*180/pi)
    end
    
    range_bins = good_rbins
    MPC_rbins = MPC_rbins(:,1)
  end
  
  if 0
    % Phase difference between sensors that we expect to see in the phase
    % of the DCM. It is calculated from geometry.
%     ypc = -phase_center(2,:).';
%     zpc = -phase_center(3,:).';
    clear max_expected_phase
    for chan_idx = 1:Nc
      % Distance from the a given sensor to the reference sensor
      L_chan      = sqrt(abs(diff(ypc([ref_chan  chan_idx])))^2 + abs(diff(zpc([ref_chan  chan_idx])))^2);
      % Array length in the direction of DOA
      L_doa   = L_chan*cos(DOA);
      % Extra delay with respect to the reference sensor
      L_extra = L_chan*sin(DOA);
      max_expected_phase(chan_idx,1) = 4*pi/lambda * L_extra * 180/pi;
    end
    max_expected_phase 
    
    DCM_phase_UTM = zeros(Nc);
    for chan_idx = 1:Nc     
      DCM_phase_UTM(chan_idx,chan_idx:end) = max_expected_phase(1:Nc-chan_idx+1);
    end
    DCM_phase_LTM = -triu(DCM_phase_UTM).';
    expected_DCM_phase = DCM_phase_UTM + DCM_phase_LTM;
  end
  %% Plots
  % Plot the magnitude and phase of the DCM
  % ---------------------------------------
  figure(3000);clf
  suptitle(sprintf('2D sim: DOA = %2.2f deg.\n',DOA(1)*180/pi))
  subplot(211)
  imagesc(10*log10(abs(DCM_tmp)./max(abs(DCM_tmp(:)))))
%   xlabel('Sensor index')
  ylabel('Sensor index')
  h = colorbar;
  ylabel(h,'Normalized |R| (dB)')
  colormap parula
  
  subplot(212)
  DCM_phase = angle(DCM_tmp);
%   DCM_phase = unwrap(DCM_phase,[],2);
%   DCM_phase = DCM_phase - diag(diag(DCM_phase));
  imagesc(DCM_phase*180/pi);
%   imagesc(DCM_phase*180/pi)
%   xlabel('Sensor index')
  ylabel('Sensor index')
  h = colorbar;
  ylabel(h,'\angle R (\circ)')
  colormap parula
   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  if ~isempty(Ticks)
    h.Ticks = Ticks;
  end
  
  % Plot the phase of the DCM derived from array geometry
  % -----------------------------------------------------
  if 0
  figure(3001);clf
  suptitle(sprintf('From geometry: DOA = %2.2f deg.\n',DOA(1)*180/pi))
  subplot(211)
  imagesc(expected_DCM_phase)
%   xlabel('Sensor index')
  ylabel('Sensor index')
  h = colorbar;
  ylabel(h,'\angle R (\circ) - Before unwrapping')
  colormap parula
   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  if ~isempty(Ticks)
    h.Ticks = Ticks;
  end
  
  subplot(212)
  imagesc(wrapToPi(expected_DCM_phase*pi/180)*180/pi)
  xlabel('Sensor index')
  ylabel('Sensor index')
  h = colorbar;
  ylabel(h,'\angle R (\circ) - After unwrapping')
  colormap parula
   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  if ~isempty(Ticks)
    h.Ticks = Ticks;
  end  
  end
  
  % Plot the eigenvalues
  % --------------------
  figure(3002);clf
  stem(10*log10(mean_eigval),'filled','LineWidth',2)
  xlabel('Eigenvalue index')
  ylabel('Normalized eigenvalue (dB)')
  title(sprintf('2D sim: DOA = %2.2f deg.',DOA(1)*180/pi))
%   title('Mean eigenvalues of R')
  grid on
  grid minor
  
  % Plot eigenvaectors spectrum (or MUSIC pseudo-spectrum)
  % ----------------------------------------------------
  theta_sv = linspace(-90,90,2048)*pi/180;
  
  % Create steering vectors  
  sv_params = [];
%   phase_center = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
  sv_params.src.y_pc = -phase_center(2,:).';
  sv_params.src.z_pc = -phase_center(3,:).';
  sv_params.src.fc   = fc;
  
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs = A(theta_sv,sv_params)./sqrt(Nc);
  
  % MUSIC pesudo-spectrum
  [V, D] = eig(DCM);
  [D, D_idxs] = sort(real(diag(D)),'descend');
  V = V(:,D_idxs);
  
  N = 1;
  f_music = 1./sum(abs(V(:,N+1:end)' * SVs).^2,1);
  f_music = 10*log10(f_music./max(abs(f_music)));
  
  figure(3003);clf
  plot(theta_sv*180/pi,f_music,'b')
  xlim(DOA(1)*180/pi+[-60  60])
  xlabel('\theta ^\circ')
  ylabel('Pseudo-power (dB)')
  title(sprintf('2D sim: DOA = %2.2f deg.',DOA(1)*180/pi))
  grid on
  grid minor 
  
  % Plot array geometry
  % -------------------
  figure(3004);clf
%   ypc = -phase_center(2,:).';
%   zpc = -phase_center(3,:).';
  plot(ypc,zpc ,'b*','LineWidth',2)
  if max(ypc) > min(ypc)
    xlim([min(ypc)  max(ypc)])
  end
  if max(zpc) > min(zpc)
    ylim([min(zpc)  max(zpc)])
  end
  
  grid on
  xlabel('y-axis')
  ylabel('z-axis')
  title(sprintf('Phase centers of the sensors for %2.2f deg. roll angle',theta_roll*180/pi))  
  
  if 0
    % Plot the eigenvalues similar to 1D simulator (i.e. different
    % range-bins are the snapshots, even though there is correlation
    % between range snapshots here, which is ignored in the 1D case)
    sample_data = squeeze(sim_data(:,1,:)).';
    DCM = 1/size(sample_data,2) * sample_data*sample_data';
    eigenvalue = real(sort(eig(DCM),'descend'));
   
    figure(3);clf
    stem(10*log10(eigenvalue),'filled','LineWidth',2)
    xlabel('Eigenvalue index')
    ylabel('Normalized eigenvalue (dB)')
    title('1D-equivalent eigenvalues of R')
    grid on
    grid minor
  end
  
  if param.debug_level >= 3
    return
  end
  return
  
end

if 0
  %% Studying the effect roll and the actual phase of the DCM: for array calibration prurpos
  % ----------------------------------------------------------------------------------------
    %% Setup simulation parameters
  param = [];
  
  % Debug level of 3 causes this function to stop early and output
  % simulated data, simulation parameters, and array processing results
  % in a "results" structure.
  param.debug_level = 3;
  
  %% Source parameters
  fc = 195e6;
  BW = 30e6;
  % Nc: number of sensors
  Nc = 7;
  flight_height  = 0;
 
    lambda = c/fc;

    DOA = [20]*pi/180;
    
    % Main target location
    z_pos = -1500;
    Range = abs(z_pos)/cos(DOA); % Range to the main target
    y_pos = Range*sin(DOA);      % y-axis of the main target
  
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.t0                      = 2*(abs(z_pos)-500)/c;
  param.src.t1                      = 2*(abs(z_pos)+1000)/c;
  if param.src.t0 < 0
    param.src.t0 = 0;
  end
  param.src.ft_func                 = @(t) tukeywin_cont(t * BW);
  param.src.ft_wind                 = @(N) hanning(N);
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  
  % Nsep: number of lambda/4 steps between sensors
  Nsep = 1;
  param.src.lever_arm.fh_args       = {[], 1, 1:Nc, [0; Nsep*c/fc/2; 0]};
  phase_center = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
  
  % Change reference sensor if you want to
  ref_chan = 1;
  phase_center(2,:) = phase_center(2,:) - phase_center(2,ref_chan);
  phase_center(3,:) = phase_center(3,:) - phase_center(3,ref_chan);
 
  % Account for roll angle. Default is 0
%   theta_roll = DOA; 
  theta_roll = 0 * pi/180;
  for chan_idx = 1:Nc
    d_chan_ref = sqrt(abs(diff(phase_center(2,[ref_chan  chan_idx])))^2 + abs(diff(phase_center(3,[ref_chan  chan_idx])))^2);
    phase_center(2,chan_idx) = d_chan_ref .* cos(theta_roll);
    phase_center(3,chan_idx) = d_chan_ref .* sin(theta_roll);
  end
  
  param.src.phase_center = phase_center;
  ypc = -phase_center(2,:).';
  zpc = -phase_center(3,:).';
  
  dt = 1/BW;
  % Nt: number of range bins (fast time samples)
  Nt = floor((param.src.t1-param.src.t0)/dt);
  
  % Time = fast time axis (sec)
  Time = param.src.t0 + dt*(0:Nt-1).';
  %   param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
  param.src.noise_power             = zeros(1,Nc);
  
  % DOA method parameters
  param.method.list                   = [7];
  
  %% Simulation Runs Setup
  
  % Cross track monte carlo setup
  param.monte.target_func = @sim.surface_gen;
  param.monte.runs = 1;
  param.monte.random_seed_offset = 0;
  
  % Target surface parameters
  surf_param = [];
  surf_param.z.mean = z_pos;
  surf_param.z.rms_height = 0;
  surf_param.z.corr_length_x = 400;
  surf_param.z.corr_length_y = 400;
  surf_param.rcs_in.mean = 0;
  surf_param.rcs_in.var = 1e2;
  surf_param.dy = 10;
  %   surf_param.y_range = [-2500 2500];
  surf_param.dx = 10;
  surf_param.x_range = [-2500 2500];
  surf_param.x = [-1500:-1500+100]; % 199...100 range-lines
  %   surf_param.y = [-1500:20:1500].';
  surf_param.y = y_pos;
  
  % y_range should >= maximum y (i.e. targets shoulb be inside the imaged
  % swath)
  y_range = 2*Range;
  if y_range == 0
    surf_param.y_range = [-1000 1000];
  else
    surf_param.y_range = [-1.5*y_range 1.5*y_range];
  end
  param.monte.target_param{1} = surf_param;
  
  %% Array Processing parameters
  array_param = [];
  
  fc = (param.src.f0+param.src.f1)/2;
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
  theta = [theta,pi/2];
  array_param.Nsv = {'theta', asin(ky/k)};
  
  array_param.sv_fh = @array_proc_sv;
  
  array_param.dbin = 1;
  array_param.dline = 1;
  
  array_param.bin_rng = 0;
  array_param.line_rng = -50:50;
  
  array_param.Nsrc = 2;
  Nsrc = array_param.Nsrc;
  
  array_param.init = 'ap';
  array_param.doa_theta_guard = 2*pi/180; %(max(theta)-min(theta))/(4*Nc);
  
  array_param.Nsubband = 1;
  dt = 1/fs;
  array_param.imp_resp.time_vec = -3*array_param.Nsubband*dt : dt/8 : 3*array_param.Nsubband*dt;
  BW = abs(param.src.f1 - param.src.f0);
  array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / BW);
  
  for idx = 1:array_param.Nsrc
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
  results = crosstrack(param);
  
  sim_data = squeeze(results.sim_data{1}); % Nt-by-Nx-by-Nc  
%   est_doa         = results.tomo.doa;
  actual_doa_cell = results.actual_doa;
  array_param     = results.array_param;
  
  % Convert actual_doa from cell array into a matrix (maximum Nsrc targets
  % per range-bin).
  actual_doa = NaN(length(array_param.bins),Nsrc,length(array_param.lines));
  for lineIdx = 1:length(array_param.lines)
    for binIdx = 1:length(array_param.bins)
      if ~isempty(actual_doa_cell{binIdx,lineIdx})
        doa_tmp = actual_doa_cell{binIdx,lineIdx};
        if length(doa_tmp)<Nsrc
          doa = NaN(Nsrc,1);
          doa(1:length(doa_tmp)) = doa_tmp;
        else
          doa = doa_tmp(1:Nsrc);
        end
        actual_doa(binIdx,1:Nsrc,lineIdx) = actual_doa_cell{binIdx,lineIdx}; % Nt*Nsrc*Nx
      end
    end
  end
  
  % Calculate the eigenvalues of each test range-bin
  eigenvalue_all  = [];
  sample_data_all = [];
  sample_data     = [];
  for lineIdx_idx = 1:length(array_param.lines)
    line_idx = array_param.lines(lineIdx_idx);
    good_rbins = sort(find(~isnan(actual_doa(:,2,lineIdx_idx))),'ascend');
%     first_rbin = MPC_rbins(1);
    rqd_rbin = good_rbins(1);%min(good_rbins);
    range_bin_idxs = rqd_rbin + [0];
    for binIdx_idx = 1:length(range_bin_idxs)
      bin_idx = range_bin_idxs(binIdx_idx);
      sample_data = squeeze(sim_data(bin_idx,line_idx+array_param.line_rng,:)).';
      
      if 0
        % FB averaging and spatial smoothing
        DCM = 1/size(sample_data,2) * sample_data*sample_data';
        Rxx1 = 1/(size(sample_data(1:end-1,:),2))*sample_data(1:end-1,:) * sample_data(1:end-1,:)';
        Rxx2 = 1/(size(sample_data(2:end,:),2))*sample_data(2:end,:) * sample_data(2:end,:)';
        
        % Apply FB averaging for each subarray
        reflect_mtx = flipud(eye(size(Rxx1)));
        Rxx1_fb = (1/2)*(Rxx1 + reflect_mtx * conj(Rxx1) * reflect_mtx);
        Rxx2_fb = (1/2)*(Rxx2 + reflect_mtx * conj(Rxx2) * reflect_mtx);
        
        % Average the two DCMs
        Rxx_ss = 1/2*(Rxx1_fb+Rxx2_fb); % (Nc-1)-by-(Nc-1)
        
        % Handle the lost sensor (or DOF) such that the final DCM is Nc-by-Nc
        Rxx_tmp = zeros(Nc);
        Rxx_tmp(1:end-1,1:end-1) = Rxx_ss;
        Rxx_tmp(end,:) = DCM(end,:);
        Rxx_tmp(:,end) = DCM(:,end);
        
        DCM = Rxx_tmp; % Nc-by-Nc matrix
      elseif 0
        % FB averaging only
        DCM = 1/size(sample_data,2) * sample_data*sample_data';
        reflect_mtx = flipud(eye(size(Rxx_calib)));
        DCM = (1/2)*(DCM + reflect_mtx * conj(DCM) * reflect_mtx);
      else
        DCM = 1/size(sample_data,2) * sample_data*sample_data';
      end
      if bin_idx == rqd_rbin
        DCM_tmp = DCM;
      end
      
      eigenvalue = real(sort(eig(DCM),'descend'));
      eigenvalue_all(:,end+1) = eigenvalue;
      
      sample_data_all  = [sample_data_all;sample_data(:)];
    end
  end
  mean_eigval = mean(real(eigenvalue_all),2);
  mean_eigval = mean_eigval./max(mean_eigval);
  
  if 1
    % Phase difference between sensors that we expect to see in the phase
    % of the DCM. It is calculated from geometry.
%     ypc = -phase_center(2,:).';
%     zpc = -phase_center(3,:).';
    
% DOA_roll is the DOA wrt the range vector connecting the array and target
% (i.e. the new nadir).
    if theta_roll >= DOA
       DOA_roll = theta_roll - DOA;  
    elseif theta_roll < DOA
       DOA_roll = DOA - theta_roll; 
    end
    
    clear max_expected_phase
    for chan_idx = 1:Nc
      % Distance from the a given sensor to the reference sensor
      L_chan      = sqrt(abs(diff(ypc([ref_chan  chan_idx])))^2 + abs(diff(zpc([ref_chan  chan_idx])))^2);
      % Array length in the direction of DOA
      L_doa   = L_chan*cos(DOA_roll);
      % Extra delay with respect to the reference sensor
      L_extra = L_chan*sin(DOA_roll);
      max_expected_phase(chan_idx,1) = 4*pi/lambda * L_extra * 180/pi;
    end
    max_expected_phase 
    
    DCM_phase_UTM = zeros(Nc);
    for chan_idx = 1:Nc     
      DCM_phase_UTM(chan_idx,chan_idx:end) = max_expected_phase(1:Nc-chan_idx+1);
    end
    DCM_phase_LTM = -triu(DCM_phase_UTM).';
    expected_DCM_phase = DCM_phase_UTM + DCM_phase_LTM;
  end
  
  % Calculate TBP across the array in the direction of DOA
  TBP = 2*L_doa/c * BW;
  
  sprintf('TBP across the array in the direction of %2.1f deg. is %2.2f',DOA*180/pi,TBP)
  
  %% Plots
  % Plot the magnitude and phase of the DCM
  % ---------------------------------------
  figure(3000);clf
  suptitle(sprintf('2D sim: DOA = %2.2f deg.\n',DOA(1)*180/pi))
  subplot(211)
  imagesc(10*log10(abs(DCM_tmp)./max(abs(DCM_tmp(:)))))
%   xlabel('Sensor index')
  ylabel('Sensor index')
  h = colorbar;
  ylabel(h,'Normalized |R| (dB)')
  colormap parula
  
  subplot(212)
  DCM_phase = angle(DCM_tmp);
  % Unwrap along both dimsnsions
  DCM_phase = unwrap(DCM_phase,[],2);
  DCM_phase = DCM_phase - diag(diag(DCM_phase));
  
  DCM_phase = unwrap(DCM_phase,[],1);
  DCM_phase = DCM_phase - diag(diag(DCM_phase));
  
  imagesc(DCM_phase*180/pi);
%   imagesc(DCM_phase*180/pi)
%   xlabel('Sensor index')
  ylabel('Sensor index')
  h = colorbar;
  ylabel(h,'\angle R (\circ)')
  colormap parula
   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  if ~isempty(Ticks)
    h.Ticks = Ticks;
  end
  
  % Plot the phase of the DCM derived from array geometry
  % -----------------------------------------------------
  figure(3001);clf
  suptitle(sprintf('From geometry: DOA = %2.2f deg.\n',DOA(1)*180/pi))
  subplot(211)
  imagesc(expected_DCM_phase)
%   xlabel('Sensor index')
  ylabel('Sensor index')
  h = colorbar;
  ylabel(h,'\angle R (\circ) - Before unwrapping')
  colormap parula
   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  if ~isempty(Ticks)
    h.Ticks = Ticks;
  end
  
  subplot(212)
  imagesc(wrapToPi(expected_DCM_phase*pi/180)*180/pi)
  xlabel('Sensor index')
  ylabel('Sensor index')
  h = colorbar;
  ylabel(h,'\angle R (\circ) - After unwrapping')
  colormap parula
   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  if ~isempty(Ticks)
    h.Ticks = Ticks;
  end  
  
  % Plot the eigenvalues
  % --------------------
  figure(3002);clf
  stem(10*log10(mean_eigval),'filled','LineWidth',2)
  xlabel('Eigenvalue index')
  ylabel('Normalized eigenvalue (dB)')
  title(sprintf('2D sim: DOA = %2.2f deg.',DOA(1)*180/pi))
%   title('Mean eigenvalues of R')
  grid on
  grid minor
  
  % Plot eigenvaectors spectrum (or MUSIC pseudo-spectrum)
  % ----------------------------------------------------
  theta_sv = linspace(-90,90,2048)*pi/180;
  
  % Create steering vectors  
  sv_params = [];
%   phase_center = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
  sv_params.src.y_pc = ypc;
  sv_params.src.z_pc = zpc;
  sv_params.src.fc   = fc;
  
  A = @(theta,param)steering_mtx(theta,sv_params);
  SVs = A(theta_sv,sv_params)./sqrt(Nc);
  
  % MUSIC pesudo-spectrum
  [V, D] = eig(DCM);
  [D, D_idxs] = sort(real(diag(D)),'descend');
  V = V(:,D_idxs);
  
  N = 1;
  f_music = 1./sum(abs(V(:,N+1:end)' * SVs).^2,1);
  f_music = 10*log10(f_music./max(abs(f_music)));
  
  figure(3003);clf
  plot(theta_sv*180/pi,f_music,'b')
  xlim(DOA(1)*180/pi+[-60  60])
  xlabel('\theta ^\circ')
  ylabel('Pseudo-power (dB)')
  title(sprintf('2D sim: DOA = %2.2f deg.',DOA(1)*180/pi))
  grid on
  grid minor 
  
  % Plot array geometry
  % -------------------
  figure(3004);clf
%   ypc = -phase_center(2,:).';
%   zpc = -phase_center(3,:).';
  plot(ypc,zpc ,'b*','LineWidth',2)
  if max(ypc) > min(ypc)
    xlim([min(ypc)  max(ypc)])
  end
  if max(zpc) > min(zpc)
    ylim([min(zpc)  max(zpc)])
  end
  
  grid on
  xlabel('y-axis')
  ylabel('z-axis')
  title(sprintf('Phase centers of the sensors for %2.2f deg. roll angle',theta_roll*180/pi))  
  
  if 0
    % Plot the eigenvalues similar to 1D simulator (i.e. different
    % range-bins are the snapshots, even though there is correlation
    % between range snapshots here, which is ignored in the 1D case)
    sample_data = squeeze(sim_data(:,1,:)).';
    DCM = 1/size(sample_data,2) * sample_data*sample_data';
    eigenvalue = real(sort(eig(DCM),'descend'));
   
    figure(3);clf
    stem(10*log10(eigenvalue),'filled','LineWidth',2)
    xlabel('Eigenvalue index')
    ylabel('Normalized eigenvalue (dB)')
    title('1D-equivalent eigenvalues of R')
    grid on
    grid minor
  end
  
  if param.debug_level >= 3
    return
  end
  return
  
end

%%    Array Calibration: Effect of array errors on estimated DoAs
% =======================================================================
if 0
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
  
  return
end

%%    Array Calibration: The new implementation
% =======================================================================
if 0
  param = [];
  
  % Debug level of 3 causes this function to stop early and output
  % simulated data, simulation parameters, and array processing results
  % in a "results" structure.
  param.debug_level = 3;
  
  % -----------------------------------------------------------------------
  %                          Source parameters
  % -----------------------------------------------------------------------
  fc = 195e6;
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
targets_doa(abs(targets_doa)>40*pi/180) = [];
 
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
  phase_center(3,:) = z_pc;
  param.src.phase_center = phase_center;
  
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
  surf_param.rcs_in.var = 1e2;
  surf_param.dy = 10;
  %   surf_param.y_range = [-2500 2500];
  surf_param.dx = 10;
  surf_param.x_range = [-2500 2500];
  %   surf_param.x = [-500:1:500]; % Defined later
  %   surf_param.y = [-1500:20:1500].';
  surf_param_y = range_vec .* sin(targets_doa);
  
  % y_range should >= maximum y (i.e. targets shoulb be inside the imaged
  % swath)
  surf_param.y_range = [-1.5*max(surf_param_y) 1.5*max(surf_param_y)];
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
  array_param.line_rng = -10:10;
  
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
               ;-inf  inf; 0  0; 0  0];

%    
%   error_ypc      = [-0.1943    0.0452   -0.1972         0   -0.1958   -0.1943    0.0875]';%* lambda; % In meters units
%   error_zpc      = [-0.1924    0.0452    0.1944         0    0.1965   -0.0128   -0.1859]';%* lambda; % In meters units
%   error_phase    = [-1.1228    0.6236    1.6716         0    1.6941    0.1720   -1.0669]';
%   error_g_s      = [-3.7225    4.5863    4.0187         0    9.8245    4.0963    3.7311]';
%   error_g_p      = [-1.4270   -1.6135   -1.7050         0   -0.0157   -1.9984   -1.5337]';
%   error_g_offset = [-9.8911    6.6751    8.8213         0   -9.7071    9.6572    9.7167]';
    
    % y_pc and z_pc errors are generated as percentage fo the actual y_pc
    % and z_pc values. This guarantees that the errors are sensible.
    y_pc_err_percentage = [1 5 0 4 9 2 8].' ./100;
    z_pc_err_percentage = [1  7  1  1  3  5  2].' ./100;
    error_ypc      = y_pc .* y_pc_err_percentage;%  [0.02 0.01 -0.03 0 0.009 0.002 -0.001]';% * lambda; % In meters units
    error_zpc      = z_pc .* z_pc_err_percentage;% [0 0.001 0.001 -0.02 0.005 0.01 0.001]';% *lambda ; % In meters units
    error_phase    = zeros(Nc,1);
    error_g_s      = [0.2 0.8 1 0.9 1 0.8 1]';
    error_g_p      = zeros(Nc,1);
    error_g_offset = zeros(Nc,1);
    
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
    
  % -----------------------------------------------------------------------
  %                       Run the simulation
  % -----------------------------------------------------------------------  
  est_doa_before_calib = [];
  Rxx_all = [];
  actual_doa = [];
  for target_i = 1:length(targets_doa)
    doa_tmp = targets_doa(target_i);
    
    surf_param.y = surf_param_y(target_i);
    param.monte.target_param{1} = surf_param;
    
    results = crosstrack(param);
    
     
    good_doa_i = find(~isnan(results.tomo.doa(:)));
    if ~isempty(good_doa_i)
      est_doa_before_calib(end+1) = results.tomo.doa(good_doa_i);
      sim_data = squeeze(results.sim_data{1}(good_doa_i,:,:)).'; % Nt*Nx*Nc data matrix
      Data{target_i} = sim_data;
      Rxx_all{end+1} = 1/size(sim_data,2) * sim_data * sim_data';
      actual_doa(end+1) = targets_doa(target_i);
    end
  end
   actual_doa_deg = actual_doa*180/pi;

   
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
  ac_cost_params.y_pc       = repmat(y_pc,[1 length(actual_doa)]);
  ac_cost_params.z_pc       = repmat(z_pc,[1 length(actual_doa)]);
  ac_cost_params.fc         = fc;
  ac_cost_params.BW = BW;
  ac_cost_params.Rxx_all = Rxx_all;
%   ac_cost_params.sim_data   = sim_data;
  ac_cost_params.array_param = array_param;
%   ac_cost_params.array_gain_model_fh = param.array_gain_model_fh;
  
  tic
  % -----------------------------------------------------------------------
  %                          Run the optimizer
  % -----------------------------------------------------------------------
  LB = repmat(error_bounds(:,1).',[Nc 1]);
  UB = repmat(error_bounds(:,2).',[Nc 1]);
  
  ref_chan = 1;
  LB(ref_chan,:) = 0;
  UB(ref_chan,:) = 0;
%   LB(ceil(Nc/2),:) = 0;
%   UB(ceil(Nc/2),:) = 0;
  
  LB = LB(:);
  UB = UB(:);
  initial_ac = 0*ones(size(LB)); % Nc*6 matrix
  
  [nadir_doa,nadir_doa_idx] = min(abs(targets_doa));
  
  param.fc = fc;
  param.nadir_y_pc = y_pc;
  param.nadir_z_pc = z_pc;
  param.nadir_doa  = nadir_doa;
  param.doa = targets_doa.';
  array_calib_nonlcon_fh = []; %@(est_errors) array_calib_nonlcon(est_errors,param);
  
  if 1
    % Local solver (fmincon solver)
    options =  optimoptions(@fmincon,'TolX',1e-6,'TolFun',1e-6,'MaxIter',10^7,'MaxFunEvals',10^7, ...
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
    calib_zpc = uncalib_zpc + calib_params.calib_zpc;
    
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
      phase_deg = (calib_params.calib_ypc*k*sin(actual_doa(doa_idx)) - calib_params.calib_zpc*k*cos(actual_doa(doa_idx)) + calib_params.calib_phase)*180/pi;
      
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
      suptitle(sprintf('Eigenvectors pattern: Nsrc=%2d',N))

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
      suptitle(sprintf('Eigenvectors pattern: Nsrc=%2d',N))
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
        repmat(calib_params.calib_phase,[1 length(theta_RP)]))*180/pi;
      
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
end

%%                Model order estimation
% =========================================================================
if 1
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
    LL_results = crosstrack(param);
%     LL_results = sim.crosstrack(param);
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
      LL_results = crosstrack(param);
%       LL_results = sim.crosstrack(param);
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
    
    results = crosstrack(param);
%     results = sim.crosstrack(param);
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
  return;
end


%% Particle Filter
% =========================================================================
if 0
  fprintf('=======================================================================\n');
  fprintf('Running Particle filter simulation\n');
  fprintf('=======================================================================\n');
  
line_rng =[-10:10];% round([-1*linspace(1,125,51).' , linspace(1,125,51).']);   % Ntests by 2 (start and end of range-lins range)
Ntests = size(line_rng,1);

% Bf = [1e-4:0.2:10].';
% Ntests = length(Bf);
% Ntests = 1;
test_param = [];
rmse_tests = [];
  for test_idx = 1:Ntests
  %% Setup simulation parameters
  param = [];
  
  % Debug level of 3 causes this function to stop early and output
  % simulated data, simulation parameters, and array processing results
  % in a "results" structure.
  param.debug_level = 3;
  
  %% Source parameters
  fc = 195e9;
  BW = 32e6;
%   fc = BW/Bf(test_idx);
  flight_height                     = 1500;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.t0                      = 2*(flight_height-500)/c;
  param.src.t1                      = 2*(flight_height+1000)/c;
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
  param.src.noise_power             = zeros(1,Nc);
%   param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
  
  % DOA method 10 is the particle filter
  param.method.list = [10];
  method.name       = 'PF';
  
 %% Place targets on the surface at specific angles
% ang = [-1.2154 -1.0654 -0.9484 -0.8481 -0.7580 -0.6751 -0.5974 -0.5236 -0.4528 -0.3844 -0.3178 -0.2527 -0.1886 -0.1253 -0.0625...         0
%      0 0.0625 0.1253 0.1886 0.2527 0.3178 0.3844 0.4528 0.5236 0.5974 0.6751 0.7580 0.8481 0.9484 1.0654 1.2154]';
% ang = [0 9 16 19 22 26 30]'*pi/180;
% ang = [-[17 15 11 8 5], [0 5 8 11 15 17]]'*pi/180;
ang = [-[17 15 13 11 9 7 5 3], [0 3 5 7 9 11 13 15 17]]'*pi/180;
% ang = [-10 10]'*pi/180;

  %% Simulation Runs Setup
  
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
  surf_param.dy = 10;
  surf_param.y_range = [-2500 2500];
  surf_param.dx = 10;
  surf_param.x_range = [-2500 2500];
%   surf_param.x = [-500:5:500];
%   surf_param.y = [-1500:20:1500].';
   surf_param.y = flight_height*tan(ang);
%   param.monte.target_param{1} = surf_param;
  
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
%   array_param.line_rng = -10:10;
array_param.line_rng = round([line_rng(test_idx,1) : line_rng(test_idx,2)]);
  
  N_skipped_rlines = 1;  
  N_reqd_rlines    = 1;  
  surf_param.x = [-2000:N_skipped_rlines:-2000+N_skipped_rlines*(N_reqd_rlines+2*max(array_param.line_rng))-1];
  param.monte.target_param{1} = surf_param;
  % In the PF contest, Nsrc is the maximum number of targets.  
  array_param.Nsrc = 2;
  
%   array_param.init = 'ap';
  array_param.doa_theta_guard = (max(theta)-min(theta))/(4*Nc);
  
  array_param.Nsubband = 1;
  dt = 1/fs;
  array_param.imp_resp.time_vec = -3*array_param.Nsubband*dt : dt/8 : 3*array_param.Nsubband*dt;
  BW = abs(param.src.f1 - param.src.f0);
  array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / BW);
  
  for idx = 1:array_param.Nsrc
    array_param.doa_constraints(idx).method = 'fixed';
    array_param.doa_constraints(idx).init_src_limits = [min(theta) max(theta)]*180/pi;
    array_param.doa_constraints(idx).src_limits = [min(theta) max(theta)]*180/pi;
  end
  
  param.array_param = array_param;
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
  
  %% Particle filter parameters
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
  
  return;
end

