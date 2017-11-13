% script sim.crosstrack_example
%
% Example setup scripts for the sim.crosstrack function. To run, enable one
% of the if {0|1} blocks below and then run the script.
%
% Author: Sean Holloway, John Paden, Theresa Stumpf
%
% See also: sim.crosstrack.m, sim.crosstrack_data.m, sim.crosstrack_example.m


clear
close
load params_debug param_debug
param_debug.l_norm_all = [];
param_debug.l_all = [];
param_debug.subopt = [];
save params_debug param_debug

physical_constants;

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
  array_param.Nsv = {'theta', asin(ky/k)};  
  
  array_param.sv_fh = @array_proc_sv;
  
  array_param.dbin = 1;
  array_param.dline = 1;
  
  array_param.bin_rng = -1:1;
  array_param.rline_rng = -10:10;
  
  
  %%%%%%%%%%%%
  array_param.Nsig = 1:Nc-1;
  
  array_param.init = 'ap';
  array_param.theta_guard = (max(theta)-min(theta))/length(ky);
  
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

if 1
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
  %BW = 10e6;
  
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
  
  
  %CHECK NOISE POWER FORMULATION
   %param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
  param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
  
  
%   param.src.noise_power             = 10*log10(ones(1,Nc));
  
% from crosstrack_data  
%   for chan = 1:Nc
%   noise_data = 10.^(param.src.noise_power(chan)/20)*(randn(100,1,1) ...
%     +1i*randn(100,1,1))/sqrt(2);
% 
% end
% 
%  var(noise_data(:,1,1)) 
 
 
 
 % I THINK /20 IS WRONG INSTEAD WE SHOULD DO /10 above
 
  
  % DOA method parameters
  param.method.list                   = [7];
  
  %% Simulation Runs Setup
  
  % Cross track monte carlo setup
  param.monte.target_func = @sim.surface_gen;
  
  
  %RUNS
  param.monte.runs = 1;
  
  
  
  param.monte.random_seed_offset = 0;
  
  % Target surface parameters
  surf_param = [];
  surf_param.z.mean = -1500;
  surf_param.z.rms_height = 2;
  surf_param.z.corr_length_x = 400;
  surf_param.z.corr_length_y = 400;
  surf_param.rcs_in.mean = 0;
  
  
    noise_power = 10.^(param.src.noise_power(1)/10)
    SNR_db = 20
 
  surf_param.rcs_in.var = noise_power * 10.^(SNR_db./10)   
  

  
  surf_param.dy = 10;
  surf_param.y_range = [-2500 2500];
  surf_param.dx = 10;
  surf_param.x_range = [-2500 2500];
  surf_param.x = [-500:5:500];
  surf_param.y = [-1500:60:1500].';
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
  array_param.rline_rng = -10:10;
  
  
  %SRAVYA
   array_param.Nsig  = 1:Nc-1;
  % array_param.Nsig  = 2;
  array_param.init = 'ap';
  array_param.theta_guard = (max(theta)-min(theta))/(4*Nc);
  
  array_param.W = 1;
  dt = 1/fs;
  array_param.imp_resp.time_vec = -3*array_param.W*dt : dt/8 : 3*array_param.W*dt;
  BW = abs(param.src.f1 - param.src.f0);
  array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / BW);
  
  
  %SRAVYA
  %for idx = 1:2
  for idx = 1:max(array_param.Nsig)
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
  [results, sources_true_all] = sim.crosstrack(param);
  
  
  
  
  
%%   %HISTOGRAM FOR A METHOD W.R.T MODEL ORDER

% %   Nest_suboptimal= results.tomo.doa.model_order_results_suboptimal; 
% %   Nest_optimal = results.tomo.doa.model_order_results ;
% %   
% %   RLINES = 3;
% %   
% %   for method_idx = 1:13
% %   
% %       
% %       
% %       switch method_idx
% %           
% %           case 1
% %               
% %               Nest = Nest_suboptimal.AIC.Nest;
% %           case 2
% %               
% %               Nest = Nest_suboptimal.HQ.Nest;
% %           case 3
% %               
% %               Nest = Nest_suboptimal.MDL.Nest;
% %           case 4
% %               
% %               Nest = Nest_suboptimal.AICc.Nest;
% %           case 5
% %               
% %               Nest = Nest_suboptimal.KICvc.Nest;
% %           case 6
% %               
% %               Nest = Nest_suboptimal.WIC.Nest;             
% %               
% %           case 7
% %               Nest = Nest_suboptimal.AIC.Nest;
% %           case 8
% %               
% %               Nest = Nest_suboptimal.HQ.Nest;
% %           case 9
% %               
% %               Nest = Nest_suboptimal.MDL.Nest;
% %           case 10
% %               
% %               Nest = Nest_suboptimal.AICc.Nest;
% %           case 11
% %               
% %               Nest = Nest_suboptimal.KICvc.Nest;
% %           case 12
% %               
% %               Nest = Nest_suboptimal.WIC.Nest;
% %               
% %           case 13
% %               
% %               Nest = Nest_suboptimal.MDL_sravya.Nest;
% %               
% %           otherwise
% %               error('Not supported')
% %       end
% %       
% %       
% %       
% %   for k_idx = 0:Nc-1
% %   percentage(k_idx+1) = numel(find(Nest(find (sources_true_all(:,1:RLINES) == k_idx)) == k_idx)) / numel(find (sources_true_all(:,1:RLINES) == k_idx)); 
% %  
% %   end
% %   
% %   %performance can't be tested for the below case. 
% %   %since no number of sorces equal the number we are testing in true sources
% %   %percentage(find(~isnan(percentage) == 1)) = ? ;
% %    
% %   percentage_method(:,method_idx) = percentage;
% %   end
% %   

  
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
