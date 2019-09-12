% script sim.crosstrack_example_rds
%
% Example narrowband radar depth sounder simulations. Calls crosstrack.m and
% crosstrack_rmse.m.
%
% Author: Mohanad Al-Ibadi, Sravya Athinarapu, Sean Holloway, John Paden,
% Theresa Stumpf
%
% See also: sim.crosstrack.m, sim.crosstrack_data.m,
% sim.crosstrack_example_*.m, sim.crosstrack_rmse.m

sim_number = 2; % Sim 1 or Sim 2
num_layers = 2; % 1 or 2
scene_size = 'small'; % 'large' or 'small'

%% Sim 1: Radar Depth Sounder (Multichannel)
% =========================================================================
% Narrowband
% Uniform linear array
if sim_number == 1
  
  physical_constants;
  
  %% Sim 1: Setup simulation parameters
  param = [];
  
  % Debug level of 3 causes this function to stop early and output
  % simulated data, simulation parameters, and array processing results
  % in a "results" structure.
  param.debug_level = 3;
  
  %% Sim 1: Source parameters
  fc = 195e6;
  BW = 30e6;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.t0                      = 2*(1500-500)/c;
  param.src.t1                      = 2*(1500+2500)/c;
  param.src.ft_func                 = @(t) tukeywin_cont(t * BW);
  param.src.ft_wind                 = @(N) hanning(N);
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  % Nsep: number of lambda/4 steps between sensors
  Nsep = 1;
  % Nc: number of sensors
  Nc = 7;
  % param.src.lever_arm.fh_args: Arguments for a linear array in y-dimension
  param.src.lever_arm.fh_args       = {[], 1, 1:Nc, [0; Nsep*c/fc/2; 0]};
  param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
  
  % DOA method parameters
  array_proc_methods;
  param.method.list                   = [MUSIC_METHOD];
  
  %% Sim 1: Simulation Runs Setup
  
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
  surf_param.dy = 5;
  surf_param.y_range = [-5000 5000];
  surf_param.dx = 5;
  surf_param.x_range = [-5000 5000];
  if strcmpi(scene_size,'small')
    surf_param.x = [-100:5:100];
    surf_param.y = [-500:5:500].';
  else
    surf_param.x = [-2500:5:2500];
    surf_param.y = [-2500:5:2500].';
  end
  param.monte.target_param{1} = surf_param;
  
  if num_layers == 2
    surf_param = [];
    surf_param.z.mean = -2250;
    surf_param.z.rms_height = 100;
    surf_param.z.corr_length_x = 400;
    surf_param.z.corr_length_y = 400;
    surf_param.rcs_in.mean = 0;
    surf_param.rcs_in.var = 100;
    surf_param.dy = 5;
    surf_param.y_range = [-5000 5000];
    surf_param.dx = 5;
    surf_param.x_range = [-5000 5000];
    if strcmpi(scene_size,'small')
      surf_param.x = [-100:5:100];
      surf_param.y = [-500:5:500].';
    else
      surf_param.x = [-2500:5:2500];
      surf_param.y = [-2500:5:2500].';
    end
    param.monte.target_param{2} = surf_param;
  end
  
  %% Sim 1: Array Processing parameters
  param.array = [];
  
  fc = (param.src.f0+param.src.f1)/2;
  fs = param.src.f1-param.src.f0;
  
  % NN: total length of the sensor array (in lambda/4 units)
  NN = Nc*Nsep;
  % lambda: wavelength
  lambda = c/fc;
  % k: wavenumber
  k = 2*pi/(lambda/2);
  % My: over sampling factor
  My = 20;
  % dy: phase center spacing
  dy = Nsep*lambda/4;
  % dky and ky: y-component of wavenumber (spacing and axis)
  dky = 2*pi / (Nc*dy) / My;
  ky = dky * ifftshift(-floor(My*Nc/2) : floor((My*Nc-1)/2));
  % theta: theta values associated with ky axis
  theta = fftshift(asin(ky/k));
  param.array.theta = asin(ky/k)*180/pi;
  
  param.array.sv_fh = @array_proc_sv;
  
  param.array.dbin = 1;
  param.array.dline = 5;
  
  param.array.bin_rng = 0;
  param.array.line_rng = -10:10;
  
  param.array.Nsrc = 2;
  
  param.array.last_fprintf_time_delay = 5;
  param.array.init = 'ap';
  param.array.doa_theta_guard = (max(theta)-min(theta))/length(ky);
  
  param.array.Nsubband = 1;
  dt = 1/fs;
  param.array.imp_resp.time_vec = -3*param.array.Nsubband*dt : dt/8 : 3*param.array.Nsubband*dt;
  BW = abs(param.src.f1 - param.src.f0);
  param.array.imp_resp.vals = tukeywin_cont(param.array.imp_resp.time_vec / BW);
  
  for idx = 1:param.array.Nsrc
    param.array.doa_constraints(idx).method = 'fixed';
    param.array.doa_constraints(idx).init_src_limits = [min(theta) max(theta)]*180/pi;
    param.array.doa_constraints(idx).src_limits = [min(theta) max(theta)]*180/pi;
  end
  
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
  
  %% Sim 1: Run the simulation
  results = sim.crosstrack(param);
  
  if param.debug_level >= 3
    return
  end
  
  %% Sim 1: Process and save the outputs
  
  % Save Outputs
%   out_fn = ct_filename_tmp([],'crosstrack_example_rds_sim1_.fig');
%   fprintf('Saving %s\n', out_fn);
%   ct_saveas(1,out_fn);
%   if ~exist(out_fn_dir,'dir')
%     mkdir(out_fn_dir);
%   end
%   out_fn = fullfile(out_fn_dir,out_fn_name);
%   savefig(1,[out_fn '_src1.fig']);
%   savefig(2,[out_fn '_src2.fig']);
%   save([out_fn '.mat'],'param','results')
end

%% Sim 2: Radar Depth Sounder (Multichannel)
% =========================================================================
% Narrowband
% Uniform linear array
if sim_number == 2
  %% Sim 2: Setup simulation parameters
  param = [];
  
  % Debug level of 3 causes this function to stop early and output
  % simulated data, simulation parameters, and array processing results
  % in a "results" structure.
  param.debug_level = 3;
  
  %% Sim 2: Source parameters
  fc = 195e6;
  BW = 30e6;
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

  %% Sim 2: Simulation Runs Setup
  
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
  surf_param.dy = 5;
  surf_param.y_range = [-5000 5000];
  surf_param.dx = 5;
  surf_param.x_range = [-5000 5000];
  if strcmpi(scene_size,'small')
    surf_param.x = [-100:5:100];
    surf_param.y = [-500:5:500].';
  else
    surf_param.x = [-2500:5:2500];
    surf_param.y = [-2500:5:2500].';
  end
  param.monte.target_param{1} = surf_param;
  
  if num_layers == 2
    surf_param = [];
    surf_param.z.mean = -2250;
    surf_param.z.rms_height = 100;
    surf_param.z.corr_length_x = 400;
    surf_param.z.corr_length_y = 400;
    surf_param.rcs_in.mean = 0;
    surf_param.rcs_in.var = 100;
    surf_param.dy = 5;
    surf_param.y_range = [-5000 5000];
    surf_param.dx = 5;
    surf_param.x_range = [-5000 5000];
    if strcmpi(scene_size,'small')
      surf_param.x = [-100:5:100];
      surf_param.y = [-500:5:500].';
    else
      surf_param.x = [-2500:5:2500];
      surf_param.y = [-2500:5:2500].';
    end
    param.monte.target_param{2} = surf_param;
  end
  
  %% Sim 2: Array Processing parameters
  param.array = [];
  
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
  param.array.Nsv = {'theta', asin(ky/k)};
  
  param.array.sv_fh = @array_proc_sv;
  
  param.array.dbin = 1;
  param.array.dline = 5;
  
  param.array.bin_rng = 0;
  param.array.line_rng = -10:10;
  
  param.array.Nsrc = 2;
  
  param.array.last_fprintf_time_delay = 5;
  param.array.init = 'ap';
  param.array.doa_theta_guard = (max(theta)-min(theta))/(4*Nc);
  
  param.array.Nsubband = 1;
  dt = 1/fs;
  param.array.imp_resp.time_vec = -3*param.array.Nsubband*dt : dt/8 : 3*param.array.Nsubband*dt;
  BW = abs(param.src.f1 - param.src.f0);
  param.array.imp_resp.vals = tukeywin_cont(param.array.imp_resp.time_vec / BW);
  
  for idx = 1:param.array.Nsrc
    param.array.doa_constraints(idx).method = 'fixed';
    param.array.doa_constraints(idx).init_src_limits = [min(theta) max(theta)]*180/pi;
    param.array.doa_constraints(idx).src_limits = [min(theta) max(theta)]*180/pi;
  end
  
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
  
  %% Sim 2: Run the simulation
  results = sim.crosstrack(param);
  
  
  %% Sim 2: To plot Slice model
  
  layer = 1;
  grid_slice = 21;
  surface_z = results.z_grid-results.param_crosstrack.monte.target_param{1}.z.mean;
  surface_z_groung_truth = results.surf_model.z(results.surf_model.layer==layer,grid_slice)-results.param_crosstrack.monte.target_param{1}.z.mean;
  figure(5); clf;  plot(results.surf_model.y(results.surf_model.layer==layer), surface_z_groung_truth,'b');
  hold on
  plot(results.surf_model.y, surface_z(:,grid_slice),'r');
  xlim([-750 750])
  elev = results.surf_model.z(results.surf_model.layer==layer,grid_slice);
  ylim([min(elev) max(elev)]-results.param_crosstrack.monte.target_param{1}.z.mean)
  title('Slice - surface model');
  xlabel('Cross-track (m)');
  ylabel('WGS84-Elevation (m)');
  hold off
  legend('Labeled surface','Estimated surface');
  grid on
  
  % Slice - Range Bin v/s DOA
  slice = 3;
  figure(6),clf
  scatter(results.tomo.tomo.theta(:,1,slice)*(180/pi),results.param_crosstrack.array_proc.bins, 20 , 10*log10(results.tomo.tomo.img(:,1,slice)),'fill');
  colorbar
  hold on
  scatter(results.tomo.tomo.theta(:,2,slice)*(180/pi),results.param_crosstrack.array_proc.bins, 20 , 10*log10(results.tomo.tomo.img(:,2,slice)),'fill');
  colorbar
  set(gca,'Ydir','reverse')
  xlim([-60 60])
  twtt = 2/c*sqrt(results.surf_model.z(results.surf_model.layer==layer,grid_slice).^2+results.surf_model.y(results.surf_model.layer==layer).^2);
  min_twtt = min(twtt(:));
  max_twtt = max(twtt(:));
  time = results.param_crosstrack.array_proc.wfs.time(results.param_crosstrack.array_proc.bins);
  first_bin = find(time<min_twtt,1,'last');
  if isempty(first_bin)
    first_bin = 1;
  end
  last_bin = find(time>max_twtt,1);
  if isempty(last_bin)
    last_bin = length(time);
  end
  ylim([first_bin last_bin]);
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
  
  %% Sim 2: Process and save the outputs
%   out_fn_dir = 'D:\tmp\RADAR_CONF\';
%   out_fn_name = 'example1';
%   
%   RMSE = sim.crosstrack_rmse(param,results);
%   
%   figure(1); clf
%   plot(param.monte.SNR(:,1)+10*log10(3),RMSE(:,:,1).','.','LineWidth',2);
%   grid on
%   xlabel('Source SNR (dB)')
%   ylabel('RMS error ( \circ )')
%   legend({'MUSIC','MLE','WB','WBMLE'})
%   title('Figure 2 from Wax and Ziskind 1988')
%   
%   figure(2); clf
%   plot(param.monte.SNR(:,2)+10*log10(3),RMSE(:,:,2).','.','LineWidth',2);
%   grid on
%   xlabel('Source SNR (dB)')
%   ylabel('RMS error ( \circ )')
%   legend({'MUSIC','MLE','WB','WBMLE'})
%   title('Complement to figure 2 from Wax and Ziskind (source 2)')
%   
%   % Save Outputs
%   if ~exist(out_fn_dir,'dir')
%     mkdir(out_fn_dir);
%   end
%   out_fn = fullfile(out_fn_dir,out_fn_name);
%   savefig(1,[out_fn '_src1.fig']);
%   savefig(2,[out_fn '_src2.fig']);
%   save([out_fn '.mat'],'param','results')
end
