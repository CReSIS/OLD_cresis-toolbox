% script sim.doa_example
%
% Example setup scripts for the sim.doa function. To run, enable one
% of the if {0|1} blocks below and then run the script.
%
% Author: John Paden, Theresa Stumpf
%
% See also: doa.m

physical_constants;

%%   Mohanad: For 1D tutorial
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  fc = 195e9;
  BW = 32e6;
  Nc = 7;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.ft_wind                 = @(N) boxcar(N);
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:Nc],[0; c/fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [7];
  param.method.Nsv                    = 9;
  param.method.OneD_Nsv               = 128;
%   param.method.src_limits             = {[-70]/180*pi,[70]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'grid';
  param.method.wb_td.init             = 'grid';
  param.method.wb_td.widening_factor  = 1;
  param.method.wb_fd.init             = 'grid';
  param.method.wb_fd.filter_banks     = 1;
  
  % DOA monte carlo setup
  % Two equal power sources, sweep SNR
  % One source fixed at 0 deg, another fixed at 20 deg
  % Snapshots fixed at 10
  DOA = [-10 10]; % Row vector
  param.monte.SNR   = repmat(linspace(10,25,16).' - 10*log10(Nc), [1 length(DOA)]);
  num_tests = size(param.monte.SNR,1);
  param.monte.DOA   = repmat(DOA,[num_tests 1]);
  param.monte.Nsnap = repmat(101,[num_tests 1]);
  param.monte.runs  = 10;
  param.monte.random_seed_offset = 0;
  
  param.method.src_limits = [];
  for idx = 1:size(param.monte.DOA,2)
    param.method.src_limits{idx}     = [-70  70]*pi/180; 
  end

  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
%   out_fn_dir = 'E:\tmp\TSP_DOA';
%   out_fn_name = 'wax_ziskind_fig2';
  
  RMSE = sim.doa_rmse(param,results); % Nmethods-by-Ntests-by-Nsrc
  
  figure(1); clf
  plot(param.monte.SNR(:,2)+10*log10(Nc),RMSE(:,:,2).','xb','LineWidth',2);
  grid on
  xlabel('Source SNR (dB)')
  ylabel('RMS error ( \circ )')
  title('MLE: SNR vs RMSE of source 2')
  
  if param.method.list == 2
    legend('MUSIC')
  elseif param.method.list == 7
    legend('MLE')
  elseif param.method.list(1) == 2 && param.method.list(2) == 7
  legend({'MUSIC','MLE'})
  end

  
if 0
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
end
  return;
end

%% Wax and Ziskind 1988 Fig 2
% =======================================================================
if 1
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  fc = 195e6;
  BW = 30e6;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.ft_wind                 = @hanning;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:3],[0; c/fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [2 7];
  param.method.Nsv                    = 9;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[-20 40]/180*pi,[-20 40]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'grid';
  param.method.wb_td.init             = 'grid';
  param.method.wb_td.widening_factor  = 1;
  param.method.wb_fd.init             = 'grid';
  param.method.wb_fd.filter_banks     = 1;
  
  % DOA monte carlo setup
  % Two equal power sources, sweep SNR
  % One source fixed at 0 deg, another fixed at 20 deg
  % Snapshots fixed at 10
  param.monte.SNR   = repmat(linspace(10,25,16).' - 10*log10(3), [1 2]);
  num_tests = size(param.monte.SNR,1);
  param.monte.DOA   = repmat([-10 10],[num_tests 1]);
  param.monte.Nsnap = repmat(11,[num_tests 1]);
  param.monte.runs  = 1000;
  param.monte.random_seed_offset = 0;

  %% Run the simulation
  results = doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'wax_ziskind_fig2';
  
  RMSE = sim.doa_rmse(param,results);
  
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
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
%% Wax and Ziskind 1988 Fig 3
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  fc = 312.5e6;
  BW = 1e6;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:3],[0; c/fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [2 7 8 9];
  param.method.Nsv                    = 9;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[-20 40]/180*pi,[-20 40]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'grid';
  param.method.wb_td.init             = 'grid';
  param.method.wb_td.widening_factor  = 1;
  param.method.wb_fd.init             = 'grid';
  param.method.wb_fd.filter_banks     = 1;
  
  % DOA monte carlo setup
  % Two equal power 20 dB sources
  % One source fixed at 0 deg, another fixed at 20 deg
  % Sweep snapshots 10 to 1000
  param.monte.Nsnap   = round(logspace(log10(10),log10(1000),21).');
  num_tests = size(param.monte.Nsnap,1);
  param.monte.SNR   = repmat([20 20] - 10*log10(3),[num_tests 1]);
  param.monte.DOA   = repmat([0 20],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;

  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'wax_ziskind_fig3';
  
  RMSE = sim.doa_rmse(param,results);
  
  figure(1); clf
  semilogx(param.monte.Nsnap,RMSE(:,:,1).','.','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('Figure 3 from Wax and Ziskind 1988')
  
  figure(2); clf
  semilogx(param.monte.Nsnap,RMSE(:,:,2).','.','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('Complement to figure 3 from Wax and Ziskind (source 2)')
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
%% Wax and Ziskind 1988 Fig 4
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  fc = 312.5e6;
  BW = 1e6;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:3],[0; c/fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [2 7 8 9];
  param.method.Nsv                    = 9;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[-20 40]/180*pi,[-20 40]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'grid';
  param.method.wb_td.init             = 'grid';
  param.method.wb_td.widening_factor  = 1;
  param.method.wb_fd.init             = 'grid';
  param.method.wb_fd.filter_banks     = 1;
  
  % DOA monte carlo setup
  % Two equal power 20 dB sources
  % One source fixed at 0 deg, another fixed at 5 deg
  % Sweep snapshots 100 to 1000
  param.monte.Nsnap   = round(logspace(log10(100),log10(1000),11).');
  num_tests = size(param.monte.Nsnap,1);
  param.monte.SNR   = repmat([20 20] - 10*log10(3),[num_tests 1]);
  param.monte.DOA   = repmat([0 5],[num_tests 1]);
  param.monte.runs  = 500;
  param.monte.random_seed_offset = 0;

  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'wax_ziskind_fig4';
  
  RMSE = sim.doa_rmse(param,results);
  
  figure(1); clf
  plot(param.monte.Nsnap,RMSE(:,:,1).','.','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('Figure 4 from Wax and Ziskind 1988')
  
  figure(2); clf
  plot(param.monte.Nsnap,RMSE(:,:,2).','.','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('Complement to figure 4 from Wax and Ziskind (source 2)')
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
%% Wax and Ziskind 1988 Fig 6
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  fc = 312.5e6;
  BW = 1e6;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:7],[0; c/fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [2 7 8 9];
  param.method.Nsv                    = 64;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[-20 20]/180*pi,[-20 20]/180*pi,[-20 20]/180*pi,[-20 20]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 1;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 1;
  
  % DOA monte carlo setup
  % Four equal power 20 dB sources
  % Fixed at -8, -1, 5, and 15 deg
  % Sweep snapshots 100 to 1000
  param.monte.Nsnap   = round(logspace(log10(50),log10(500),11).');
  num_tests = size(param.monte.Nsnap,1);
  param.monte.SNR   = repmat([20 20 20 20],[num_tests 1]);
  param.monte.DOA   = repmat([-8 -1 5 15],[num_tests 1]);
  param.monte.runs  = 500;
  param.monte.random_seed_offset = 0;

  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'wax_ziskind_fig6';
  
  RMSE = sim.doa_rmse(param,results);
  
  figure(1); clf
  plot(param.monte.Nsnap,RMSE(:,:,3).','.','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('Figure 6 from Wax and Ziskind 1988')
  
  figure(2); clf
  plot(param.monte.Nsnap,RMSE(:,:,1).','.','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('Complement to figure 6 from Wax and Ziskind (source 1)')
  
  figure(3); clf
  plot(param.monte.Nsnap,RMSE(:,:,2).','.','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('Complement to figure 6 from Wax and Ziskind (source 2)')
  
  figure(4); clf
  plot(param.monte.Nsnap,RMSE(:,:,4).','.','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('Complement to figure 6 from Wax and Ziskind (source 4)')
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  saveas(3,[out_fn '_src3.fig']);
  saveas(4,[out_fn '_src4.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
%% TSP_DOA Figure 1
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; 0.24; 0]};
  % DOA method parameters
  param.method.list                   = [2 7 8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[-20 20]/180*pi,[-20 20]/180*pi,[-20 20]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 11;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat(linspace(0,20,11).', [1 3]);
  num_tests = size(param.monte.SNR,1);
  param.monte.DOA   = repmat([-5 0 5],[num_tests 1]);
  param.monte.Nsnap = repmat(5*11,[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;

  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);

  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','.','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('RMSE Source 1')
  
  figure(2); clf
  plot(param.monte.SNR(:,2),RMSE(:,:,2).','.','LineWidth',2);
  grid on
  xlabel('Source SNR (dB)')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('RMSE Source 2')
  
  figure(3); clf
  plot(param.monte.SNR(:,3),RMSE(:,:,3).','.','LineWidth',2);
  grid on
  xlabel('Source SNR (dB)')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('RMSE Source 3')
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  saveas(3,[out_fn '_src3.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
%% TSP_DOA Figure 2
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; 0.24; 0]};
  % DOA method parameters
  param.method.list                   = [2 7 8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[10 40]/180*pi,[45 75]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 5;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 20;
  
  % DOA monte carlo setup
  % Two sources, SNR swept from 0 to 20 dB, second source is 10 dB
  % Fixed at 35 and 55 deg
  % Snapshots fixed at 9 fast-time and 11 slow-time samples (i.e. 99)
  param.monte.SNR   = [linspace(0,20,11).', 10*ones(11,1)];
  num_tests = size(param.monte.SNR,1);
  param.monte.DOA   = repmat([25 60],[num_tests 1]);
  param.monte.Nsnap = repmat(100,[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;

  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig2';
  
  RMSE = sim.doa_rmse(param,results);

  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  grid on
  xlabel('Source SNR (dB)')
  ylabel('RMS error ( \circ )')
  %legend({'MUSIC','MLE','WDOA','WBMLE'})
  legend({'WDOA','WBMLE'})
  title('RMSE Source 1')
  
  figure(2); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,2).','+','LineWidth',2);
  grid on
  xlabel('Source SNR (dB)')
  ylabel('RMS error ( \circ )')
  %legend({'MUSIC','MLE','WDOA','WBMLE'})
  legend({'WDOA','WBMLE'})
  title('RMSE Source 2')
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
    
%% TSP_DOA Figure 3
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; 0.24; 0]};
  % DOA method parameters
  param.method.list                   = [2 7 8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[30 40]/180*pi,[50 60]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 5;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 32;
  
  % DOA monte carlo setup
  % Two sources, SNR swept from 0 to 20 dB, second source is 10 dB
  % Fixed at 35 and 55 deg
  % Snapshots fixed at 16 fast-time and 11 slow-time samples (i.e. 176)
  param.monte.SNR   = [linspace(0,20,11).', 10*ones(11,1)];
  num_tests = size(param.monte.SNR,1);
  param.monte.DOA   = repmat([35 55],[num_tests 1]);
  param.monte.Nsnap = repmat(9*11,[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;

  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig3';
  
  RMSE = sim.doa_rmse(param,results);

  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','.','LineWidth',2);
  grid on
  xlabel('Source SNR (dB)')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('RMSE Source 1')
  
  figure(2); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,2).','.','LineWidth',2);
  grid on
  xlabel('Source SNR (dB)')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('RMSE Source 2')

  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
%% TSP_DOA Figure 4
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; 0.24; 0]};
  % DOA method parameters
  param.method.list                   = [2 7 8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[-20 20]/180*pi,[-20 20]/180*pi,[-20 20]/180*pi,[-20 20]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 3;
  
  % DOA monte carlo setup
  % Four equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -9, -4, 0, 5 deg
  % Snapshots fixed at 3 fast-time and 21 slow-time samples (i.e. 63)
  param.monte.SNR   = repmat(linspace(5,25,11).', [1 4]);
  num_tests = size(param.monte.SNR,1);
  param.monte.DOA   = repmat([-9 -4 0 5],[num_tests 1]);
  param.monte.Nsnap = repmat(3*21,[num_tests 1]);
  param.monte.runs  = 300;
  param.monte.random_seed_offset = 0;

  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig4';
  
  RMSE = sim.doa_rmse(param,results);

  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','.','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('RMSE Source 1')
  
  figure(2); clf
  plot(param.monte.SNR(:,2),RMSE(:,:,2).','.','LineWidth',2);
  grid on
  xlabel('Source SNR (dB)')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('RMSE Source 2')
  
  figure(3); clf
  plot(param.monte.SNR(:,3),RMSE(:,:,3).','.','LineWidth',2);
  grid on
  xlabel('Source SNR (dB)')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('RMSE Source 3')
  
  figure(4); clf
  plot(param.monte.SNR(:,4),RMSE(:,:,4).','.','LineWidth',2);
  grid on
  xlabel('Source SNR (dB)')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('RMSE Source 4')
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  saveas(3,[out_fn '_src3.fig']);
  saveas(4,[out_fn '_src4.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

  %% 1D vs 2D comparison 
% =======================================================================
% THE RMSE CALCULATIONS HERE IS WRONG I THINK ... CHECK
if 0
  %% Setup simulation parameters
  param = [];
  physical_constants
  % Source parameters
  param.fc     = 195e9;%320e6;
  param.BW     = 10e6;%32e6;
  param.fs     = param.BW;
  param.src.f0 = param.fc-param.BW/2;
  param.src.f1 = param.fc+param.BW/2;
  Nc = 7;
  
  lambda = c/param.fc;
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
  tx_weights    = hanning(Nc);   
  if ~exist('tx_weights','var')
    tx_weights = ones(Nc,1);
  end
  
  steering_delay = (y_pc*sin(steer_ang) - z_pc*cos(steer_ang))/c;
  w = tx_weights.*exp(-1i*4*pi*param.fc*steering_delay);
   
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

  rmse_tmp = NaN(length(theta_left)+1,2);
  
  SNR = 30;
  Nruns = 20;
  % DoA = 0
  param.monte.SNR   = SNR*ones(1,1);
  param.monte.Nsnap = 3*21*ones(1,1);
  param.monte.runs  = Nruns;
  param.monte.random_seed_offset = 0;
  
  param.monte.DOA   = [0]*180/pi;
  
  results = sim.doa(param);
  
  RMSE = sim.doa_rmse(param,results);
  rmse_tmp(1,1) = squeeze(RMSE);
  rmse_tmp(1,2) = NaN;

  est_doa0 = mean(results.theta_est{param.method.list}(:,1))*180/pi;
  
  % All other DoAs
  param.monte.SNR   = SNR*ones(1,2);
  param.monte.Nsnap = 3*21*ones(1,2);
  param.monte.runs  = Nruns;
  param.monte.random_seed_offset = 0;
  
  est_doa1 = [];
  est_doa2 = [];
  for doa_idx = 2:length(theta_left)+1
    param.monte.DOA   = [theta_left(doa_idx-1) theta_right(doa_idx-1)]*180/pi;    
    results = sim.doa(param);
    
    param.monte.DOA(param.monte.DOA < plot_doa_lims(1)) = NaN;
    param.monte.DOA(param.monte.DOA > plot_doa_lims(2)) = NaN;
    RMSE = sim.doa_rmse(param,results);
    
    rmse_tmp(doa_idx,:) = squeeze(RMSE);
    est_doa1(doa_idx-1) = mean(results.theta_est{param.method.list}(:,1))*180/pi;
    est_doa2(doa_idx-1) = mean(results.theta_est{param.method.list}(:,2))*180/pi;
  end
  est_doa = [fliplr(est_doa1), est_doa0, est_doa2].';
  est_doa(est_doa<plot_doa_lims(1)) = NaN;
  est_doa(est_doa>plot_doa_lims(2)) = NaN;
  
  rmse_all = [flipud(rmse_tmp(:,1));rmse_tmp(2:end,2)];
  rmse_all(isnan(est_doa)) = NaN;    
    
  true_doa = theta*180/pi;
  
  threshold = 2 *std(rmse_all(~isnan(rmse_all)));
  rmse_all(rmse_all>threshold) = NaN;
  
  figure(10);clf;
  scatter(theta*180/pi,rmse_all,20,'fill');
  xlim([plot_doa_lims(1)-5 plot_doa_lims(2)+5])
  xlabel('True DoA (deg.)')
  ylabel('RMSE (deg.)')
  title('1D simulator')
  grid on
  
  figure(101);clf;
  scatter(true_doa,est_doa,20,'fill')
  xlim([plot_doa_lims(1)-5 plot_doa_lims(2)+5])
  ylim([plot_doa_lims(1)-5 plot_doa_lims(2)+5])
  xlabel('True DoA (deg.)')
  ylabel('Estimated DoA (deg.)')
  title('1D simulator')
   grid on
   
   return
end

%% Array calibration  
% =========================================================================
if 0
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
      [~,DCM{doa_idx+1}] = doa_wideband_data(param);
      
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
     [ac_est_error,~,exitflag] = fmincon(@(est_errors) array_calibration_cost_1D(est_errors,ac_cost_params),...
        initial_ac, [],[],[],[],LB,UB,[],options);
      exitflag = exitflag
    elseif 0
      % Local solver (patternsearch solver) .. Too slow for TolMesh>1e-2
      options =  psoptimset('TolX',1e-6,'TolFun',1e-6,'TolMesh',1e-4,'InitialMeshSize',0.001,'MaxIter',10^8,'MaxFunEvals',10^8,...
        'PlotFcn',@psplotbestf);
      % options =  psoptimset('TolX',1e-15,'InitialMeshSize',0.001,'MaxIter',10^100,'MaxFunEvals',10^100,'TolMesh',1e-15);
      ac_est_error = patternsearch(@(est_errors) array_calibration_cost_1D(est_errors,ac_cost_params),...
        initial_ac, [],[],[],[],LB,UB,[],options);

    elseif 0
      % Global solver
      % Sometimes it fails and recover .. Wast of time
      gs = GlobalSearch;
      objFun = @(est_errors) array_calibration_cost_1D(est_errors,ac_cost_params);
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
 A = @(theta,param)steering_mtx(theta,sv_params);
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
A = @(theta,param)steering_mtx(theta,sv_params);
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
  return

end

%%     Effect of array errors (uncalibrated array) on estimated DoAs
% =======================================================================

if 0
% ------------------------------------------------------------------------
%                           Setup simulation parameters
% ------------------------------------------------------------------------
  param = [];
  physical_constants
  % Source parameters
  param.fc     = 195e9;%320e6;
  param.BW     = 10e6;%32e6;
  param.fs     = param.BW;
  param.src.f0 = param.fc-param.BW/2;
  param.src.f1 = param.fc+param.BW/2;
  Nc = 7;
  
  lambda = c/param.fc;
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
  w = tx_weights.*exp(-1i*4*pi*param.fc*steering_delay);
   
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
   
   return
end

%% Multipath and Mutual coupling effects 
% =========================================================================
  if 0
  param = [];
  physical_constants
  
  % -----------------------------------------------------------------------
  % Parameters
  % -----------------------------------------------------------------------
  % Source parameters
  Nc                = 7;
  param.src.fc      = 195e6;%195e9; 
  param.src.BW      = 10e6;%10e6; 
  param.src.fs      = param.src.BW;
  param.src.f0      = param.src.fc-param.src.BW/2;
  param.src.f1      = param.src.fc+param.src.BW/2;
  param.src.ft_wind = @(N) hanning(N);
  param.src.ft_wind  = @boxcar;  

  % DOA method parameters
  param.method.list                   = [7];
  param.method.Nsv                    = 32;
  param.method.OneD_Nsv               = 128;
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 1;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 1;
   
  % -----------------------------------------------------------------------
  % Beam steering and 3dB beamwidth
  % -----------------------------------------------------------------------
  % theta contains the DoAs from a single range-bin. The goal here is to
  % test the effect of multipath on data structure. So, 1 range-bin is
  % enough.
  lambda = c/param.src.fc;
  d_y    = lambda/2;
    
  phase_center      =  zeros(3,Nc);
  phase_center(2,:) = d_y/2*(0:Nc-1);
  
  y_pc = phase_center(2,:).';
  z_pc = phase_center(3,:).';
  
  param.src.y_pc       = y_pc;
  param.src.z_pc       = z_pc;
  k = 2*pi/(lambda/2);
  
  % Transmit beamforming
  if 0
    % Steering martix function handle
    %   A = @(theta) sim.steering_mtx(theta,param);
    A = @(theta) steering_mtx(theta,param);
    
    % Test if the SVs are really orthogonal if orthogonal DoAs were used
    if 0
      SVs = A(DOA_orthogonal);
      for idx = 1:size(SVs,2)
        orth_check(idx) = SVs(:,1)'*SVs(:,idx);
      end
      orth_check = orth_check
    end
    
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
    
    idxs_3dB  = find(RP_dB>=-3.5);
    RP_3dB    = RP_dB(idxs_3dB);
    theta_3dB = theta_RP(idxs_3dB);
    
    %   beam_doa_lims = [min(theta_3dB) max(theta_3dB)];
  end
      
  % -----------------------------------------------------------------------
  % Multipath effect
  % -----------------------------------------------------------------------
  % Set the multipath componens, MPCs, associated with each actual DoA. Number of
  % multipath components can be different for each actual DoA. Assume same
  % number for now.
  
  if 1
    % Orthogonal DoAs (to test that MP is working fine)
    % Also use it to set angles and delays rather than locations of targets
    L = Nc*d_y;
    k_spacing = [-(Nc-1)/2:1:(Nc-1)/2]*(lambda/L);
    %     k_spacing = [-Nc/2:1:Nc/2-1]*(lambda/L); % John's
    DOA_orthogonal = asin(k_spacing);
%     theta = DOA_orthogonal(4);
%     theta_MP = DOA_orthogonal([3 5]);
    theta = [10]*pi/180;
    theta_MP = [-50 40]*pi/180;
    z_pos = -500;    % z-axis of the main target
    Range = abs(z_pos)/cos(theta); % Range to the main target
    y_pos = Range*sin(theta); % y-axis of the main target
    
    z_pos_MP = [z_pos+400, z_pos+480];
    y_pos_MP = z_pos_MP.*tan(theta_MP);
    
    % MPCs delay relative to the BW (or relative to the main path delay)
%     tau_MP = 1/param.src.BW*ones(1,2);
    target_to_reflector_dist   = sqrt((z_pos-z_pos_MP).^2 + (y_pos_MP-y_pos).^2);
    reflector_to_receiver_dist = sqrt(z_pos_MP.^2 + y_pos_MP.^2);
    total_distance_MP = target_to_reflector_dist + reflector_to_receiver_dist;
    tau_MP = (total_distance_MP - Range)./c; % Relative MPCs delay
  else
    theta = [+5]*pi/180;
    
    % Main path location and DoA
    z_pos = -500;    % z-axis of the main target
    Range = abs(z_pos)/cos(theta); % Range to the main target
    y_pos = Range*sin(theta); % y-axis of the main target
    
    % MPCs location, delay, and DoA
    z_pos_MP = [z_pos+400, z_pos+480]; % z-axis of the reflectors
    y_pos_MP = [y_pos+200, y_pos+200]; % y-axis of the reflectors
    Range_MP = sqrt(z_pos_MP.^2 + y_pos_MP.^2); % Range to the reflectors
    
    target_to_reflector_dist   = sqrt(z_pos_MP.^2 + abs((y_pos_MP-y_pos)).^2);
    reflector_to_receiver_dist = sqrt((z_pos-z_pos_MP).^2 + y_pos_MP.^2);
    total_distance_MP = target_to_reflector_dist + reflector_to_receiver_dist;
    tau_MP = (total_distance_MP - Range)./c; % Relative MPCs delay
    
    theta_MP = sign(y_pos_MP).*atan(y_pos_MP./(abs(z_pos_MP))); % MPCs DoA relative to the nadir
  end
  
  TBP = tau_MP/(1/param.src.BW);
  fract_BW = param.src.BW/param.src.fc;

  if 1
    % Debug: plot problem geometry
    figure(100);clf
    hold on
    h1 = plot(y_pos,z_pos,'*b','LineWidth',10);  % Main target           
    h2 = plot(y_pos_MP(1),z_pos_MP(1),'*r','LineWidth',10); % Reflector 1
    h3 = plot(y_pos_MP(2),z_pos_MP(2),'*m','LineWidth',10); % Reflector 2
    h4 = plot(0,0,'sk','LineWidth',10); % Radar
    
    plot([y_pos,y_pos_MP(1)],[z_pos,z_pos_MP(1)],'r') % Target to reflector 1
    plot([y_pos,y_pos_MP(2)],[z_pos,z_pos_MP(2)],'m') % Target to reflector 2
    plot([y_pos 0],[z_pos 0],'b')                     % Radar to target
    plot([0 y_pos_MP(1)],[0 z_pos_MP(1)],'r')         % Radar to reflector 1
    plot([0 y_pos_MP(2)],[0 z_pos_MP(2)],'m')         % Radar to reflector 2
    plot([0 0],[0 z_pos],'-.k')                       % nadir
    
    grid on
    xlabel('y-axis')
    ylabel('z-axis')
    title('Multipath problem geometry')
    legend([h1 h2 h3 h4],'Target','Reflector 1','Reflector 2','Radar','Location','southwest')
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
  
  % -----------------------------------------------------------------------
  % Run the simulator
  % -----------------------------------------------------------------------
  
  % DOA monte carlo setup 
  SNR   = 70;
  Nsnap = 1*101;
  param.monte.runs  = 1;
  param.monte.random_seed_offset = 0;
  
  param.monte.DOA = theta*180/pi;
  param.monte.SNR = SNR*ones(length(theta),1);
  param.monte.Nsnap = Nsnap*ones(length(theta),1);
  comp_tx_weight = ones(length(theta)+length(theta_MP),1);
  
  param.monte.DOA_MP = theta_MP*180/pi;
  delay_MP = [0 tau_MP];
  
  % Setting the MPC weight to -inf ignores the MPC (controls number of
  % MPCs). In dB.
  w_MP     = [0, -10 -10];

    % Transmit beamforming weights
    param.src.tx_weights = comp_tx_weight;
    % MPCs relative delays
    param.src.tau_MP = delay_MP;
    % MPCs weights (in dB)
    param.monte.w_MP = w_MP;
    % Set MCPs DoA limits
    param.method.src_limits = [];
    for MPC_idx = 1:length(param.monte.DOA)+length(param.monte.DOA_MP)
      param.method.src_limits{MPC_idx} = [-90 +90]*pi/180;
%       param.method.src_limits{idx} = [beam_doa_lims(1) beam_doa_lims(2)];
    end
   
    if 1
      clear doa
      results = doa(param);
      est_doa = squeeze(results.theta_est{param.method.list}(1,1,:));
    end
    
%     param.src.Nsnap = param.monte.Nsnap(1);
    param.src.SNR   = param.monte.SNR;
    param.src.Nsnap = param.monte.Nsnap;
    param.src.DOAs  = param.monte.DOA;
    param.src.DOAs_MP  = param.monte.DOA_MP;
    
    [Data,DCM] = doa_wideband_data(param);
    
    eigval = real(sort(eig(DCM),'descend')).';
    
%     sigma_n = mean(eigval(4:end));
%     sigma_s = mean(abs(Data(:)).^2);
%     SNR = 10*log10(sigma_s/sigma_n - 1);
     
    fprintf('\n\n TBP = [%2.4f %2.4f]\n \nFractional BW = %2.4f\n \nDoA = %2.2f\n \nMPCs DOA = [%2.2f %2.2f]\n\n',...
      TBP(1),TBP(2),fract_BW,theta*180/pi,theta_MP(1)*180/pi,theta_MP(2)*180/pi)

    if 0
    % Plot true vs estimated DoA
    figure(101);clf;
    doa_true = param.src.DOAs.';
    doa_hat = est_doa*180/pi;
    scatter(doa_true,doa_hat,20,'fill');%,'MarkerFaceColor','b')
    
    xlim([-60 +60])
    ylim([-60 +60])
    xlabel('True DoA (deg.)')
    ylabel('Estimated DoA (deg.)')
    title('1D simulator: multipath')
    grid on
    end
    
    % Plot DCM eigen values in dB
    figure(1);clf
    stem(10*log10(abs(eigval./max(eigval))),'filled','b','LineWidth',1.5)
    xlabel('Eigenvalue index')
    ylabel('Eigenvalue (dB)')
    title('Normalized eigenvalues of R')
    grid on
    
    % Plot the magnitude and phase of the DCM
    figure(2);clf
    subplot(211)
    imagesc(10*log10(abs(DCM)./max(abs(DCM(:)))))
    xlabel('Sensor index')
    ylabel('Sensor index')
    h = colorbar;
    ylabel(h,'Normalized |R| (dB)')
    colormap parula
    
    subplot(212)
    imagesc(angle(DCM)*180/pi)
    xlabel('Sensor index')
    ylabel('Sensor index')
    h = colorbar;
    ylabel(h,'Angle of R (deg.)')
    colormap parula
    
  % -----------------------------------------------------------------------
  % Plot true vs estimated DoA
  % -----------------------------------------------------------------------
%   figure(101);clf;
%   hold on
%   for doa_idx = 1:length(actual_doa)
%     doa_true = theta_multipath{doa_idx}.'*180/pi;
%     doa_hat = est_doa{doa_idx}(:)*180/pi;
%     scatter(doa_true,doa_hat,20,'fill');%,'MarkerFaceColor','b')
%   end
%   
%   xlim([-60 +60])
%   ylim([-60 +60])
%   xlabel('True DoA (deg.)')
%   ylabel('Estimated DoA (deg.)')
%   title('1D simulator: multipath')
%    grid on
  return
  end

%%  Model Order Estimation 
% =========================================================================
if 1
  % Author: John Paden, Theresa Stumpf , Sravya Athinarapu, and Mohanad Al-Ibadi
  physical_constants;
  tic
  %%
  % Set norm_saved/penalty_saved to 1 if they are already generated and saved
  norm_saved    = 0;%1;
  penalty_saved = 0;% 1;
  save_eigenvalues = 1;
  
  % opt_norm=1 if  normalization is used (normalize optimal methods). So, if
  % the normalization coefficients are not already saved, generate them.
  opt_norm = ~ norm_saved;
  
  % norm_allign_zero=1 is the best case.It means the reference for normalizing
  % loglikelihoods is q=0 case. It is done separately for suboptimal and
  % optimal.
  param_MOE.norm_allign_zero = 1;
  
  % If suboptimal/optimal_test=1, then NT will be activated (i.e.
  % optimizer_NT will be called).
  suboptimal_test = 1;
  optimal_test    = 1;
  
  % Run the optimizer to generate penalty coeffcients if you don't already
  % have them.
  D1_optimizer = ~penalty_saved;
  AICc_both    = 0;
  
  %% plots control parameters
  likelihood_plots = 1;
  stat_test        = 1;% set to 0 when do not want to do comparision of model order estimators but just test or use one of them but
  IMAGESC_plots    = 1;
  
  %%
  if AICc_both ==0
    methods = 0:6;
  elseif AICc_both == 1
    methods = 0:7;
  end
  
  % suboptimal_methods = 0:6;
  % optimal_methods = 0:6;
  
  suboptimal_methods = methods;
  optimal_methods    = methods;
  %%
  %used to save log-likelihood functions generated in optimizer_NT using
  %doa_example_NT (suboptimal) and doa_example_NT_opt (optimal)
  results =[];
  
  param_debug.eigval_all = [];
  param_debug.opt        = [];
  param_debug.subopt     = [];
  
  param_debug_NT.subopt = [];
  param_debug_NT.opt    = [];
  
  param_debug_all     = [];
  opt_log_func_all    = [];
  subopt_log_func_all = [];
  %% Setup simulation parameters
  param = [];
  Nc = 7;
  % Source parameters
  fc = 195e6;
  BW = 30e6;

  %M = Nc-1; % max number of sources we are trying to estimate (it can go max upto Nc-1)
  M = 6;%6;
  param.Nc = Nc;
  param.M = M;
  param.src.f0             = fc-BW/2;
  param.src.f1             = fc+BW/2;
  param.src.ft_wind        = @boxcar;
  param.src.lever_arm.fh   = @sim.lever_arm_example;
  param.src.lever_arm.args = {[],1,[1:Nc],[0; c/fc/2; 0]};
  % DOA method parameters
  
  if optimal_test == 0
    param.method.list               = []; %only suboptimal
  elseif optimal_test == 1
    param.method.list               = [7];
  end
  
  param.method.Nsv                    = 3*Nc; % sampling grid for mle cost
  param.method.OneD_Nsv               = 128;
  %param.method.src_limits             = {[-20 40]/180*pi,[-20 40]/180*pi};
  %param.method.src_limits             =  repmat({[-20 40]/180*pi} , [1,Nc-1])
  param.method.src_limits             =  repmat({[-80 80]/180*pi} , [1,Nc-1]);
  
  % since it Avoid evaluating cost function around previous sources and
  % hence search range decreases and can't be able to evaluate doa of more
  % sources case.
  
  %TAKE CARE OF THETA GUARD AND SPACING OF SOURCES
  %keyboard;
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 1;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 1;
  
  Nb = param.method.wb_fd.filter_banks;
  if Nb > 1 && param.method.list ~= 9
    warning('filter_banks (or number of subbans) parameter is > 1, but method.list is not 9 (or WBMLW). Either change number of subbands to 1 or change the method to WBMLE')
    keyboard
  end
  
  %param.monte.SNR   = repmat(linspace(10,25,16).' - 10*log10(3), [1 2]);
  
  num_tests = 1;
  %   param.monte.Nsnap = repmat(101,[num_tests 1]);   %%% sample size 100 used in the CHEN paper
  param.monte.Nsnap = repmat(11*Nb,[num_tests 1]);
  param.monte.runs  = 10000;      %% RUNS 10000 used in the CHEN paper
  
  % doa accessed using doa_axample , doa_axample_suboptimal using
  % this varable
  % set to 1 if the model order estimation done using doa_axample
  % it is set to 0 if the model order estimation done using
  % doa_axample_suboptimal
  
  param.doa_example = 1;
  
  % -----------------------------------------------------------------------------------
  %% Training (Q=0): generated normalizing term of the loglikelihood function
  % -----------------------------------------------------------------------------------
  % The main aprameters here are:
  % 1) param_debug.opt/subopt: contains the loglikelihoods from all runs,
  % but for a single training SNR for Q=0 case.
  % 2)  param_debug_Q_0.opt/subopt: contains the loglikelihoods from all runs
  % and all training SNRs for Q=0 case.
  % 3) param_MOE.opt_norm_term/norm_term_optimal/norm_term_suboptimal:
  % normalization coefficients wrt Q=0 case suboptimal, but applied to optimal
  % /wrt zeros applied to optimal/wrt zeros applied to suboptimal.
  
  SNR_training_Q_0 = -inf;
  
  % Generate normalization coefficients if they are not already saved
  if norm_saved == 0
    if opt_norm == 1
      param.monte.random_seed_offset = 0;
      
      param_MOE.opt_norm_term        = zeros(1,M+1);
      param_MOE.norm_term_suboptimal = zeros(1,M+1);
      param_MOE.norm_term_optimal    = zeros(1,M+1);
      
      for SNR_testing_idx =1: length(SNR_training_Q_0)
        q_idx= 0;
        for q = 0
          q_idx = q_idx +1; %USED TO SAVE PARAMS DEBUG
          DOA_true = [];
          
          param.monte.SNR   = -inf ;
          param.monte.DOA   = [];
          
          % M : Maximum number of signals
          % possible range of k (number of signals)
          k=1:M ;
          
          for Nsig_tmp = k
            % k is 1:M WHEN WE WANTED TO RUN OPTIMAL METHODS
            % for each k we should run doa.m
            param.Nsig_tmp = Nsig_tmp;
            
            %% Run the simulation
            % This section used in doa_example1
            
            %   if Nsig_tmp == 1
            %   [results{Nsig_tmp}, DCM_runs] = sim.doa(param,[]);
            %   elseif Nsig_tmp > 1
            %   doa_prev = squeeze(results{1,Nsig_tmp-1}.theta_est{1,param.method.list});
            %   [results{Nsig_tmp}, DCM_runs] = sim.doa(param, doa_prev(:,1:Nsig_tmp-1));
            %   end
            %
            clear doa
            [results{Nsig_tmp}, DCM_runs] = doa(param);
            
            %% Process and save the outputs {method}(run_idx,test_idx,1:Nsig_tmp)
            doa_tmp(1:param.monte.runs,:) = squeeze(results{1, Nsig_tmp}.theta_est{1, param.method.list});
            
            doa_mle_all{Nsig_tmp} = doa_tmp(:,1:Nsig_tmp); % each cell contain DOA for all runs. row indicates the run number.
          end
          Nruns = param.monte.runs;  % for now
          for Run_Idx = 1:Nruns;
            clear eigvals eigvecs
            if 1   % 1 for optimal
              %arranging all doa for possible k for a single run.
              for Nsig_tmp = 1:M  %%%%%CHANGE
                
                doa_mle{Nsig_tmp} = doa_mle_all{Nsig_tmp}(Run_Idx,:) ;
              end
            end
            
            for nb_idx = 1:Nb
              clear Rxx
              % DCM of this subband
              Rxx = DCM_runs{Run_Idx};
              Rxx = Rxx(Nc*(nb_idx-1)+1:nb_idx*Nc,:);
              [V,D] = eig(Rxx);
              
              [eigval,index] = sort(abs(real(diag(D))),'descend');
              
              % DCM is symmetric and always we get real eigen values. due to some
              % rounding errors in the data generated by matlab we got complex eigen
              % values
              
%               eigval_all(Run_Idx,:) = eigval;
              eigvec = V(:,index);
              
              eigvals(:,nb_idx) = eigval;
              eigvecs(Nc*(nb_idx-1)+1:nb_idx*Nc,:) = eigvec;
            end           
            %% Model Order Estimators
            % Suboptimal Methods
            model_order_suboptimal_param.Nc         = Nc;
            model_order_suboptimal_param.Nsnap      = param.monte.Nsnap;
            model_order_suboptimal_param.eigval     = eigvals;
            model_order_suboptimal_param.penalty_NT = zeros(1,M+1);
            model_order_suboptimal_param.param_MOE  = param_MOE;
            model_order_suboptimal_param.filter_banks = param.method.wb_fd.filter_banks;
            
            for model_order_method = 0%:6 <== For training, only method=0 is required
              clear sources_number doa
              model_order_suboptimal_param.method       = model_order_method;
              [~, log_penalty_cost_subopt] = model_order_suboptimal(model_order_suboptimal_param);
            end
            
            %Optimal Methods
            
            phase_center = param.src.lever_arm.fh(param.src.lever_arm.args{:});
            model_order_optimal_param.y_pc       = phase_center(2,:).';
            model_order_optimal_param.z_pc       = phase_center(3,:).';
            model_order_optimal_param.fc         = fc;
            model_order_optimal_param.Nc         = Nc;
            model_order_optimal_param.Nsnap      = param.monte.Nsnap;
            model_order_optimal_param.eigval     = eigvals;
            model_order_optimal_param.eigvec     = eigvecs;
            model_order_optimal_param.penalty_NT_opt = zeros(1,M+1);
            model_order_optimal_param.param_MOE  = param_MOE;
            model_order_optimal_param.doa_mle    = doa_mle;
            model_order_optimal_param.filter_banks = param.method.wb_fd.filter_banks;
            model_order_optimal_param.fs = BW;
            
            for model_order_method = 0;%:6 <== For training, only method=0 is required
              clear sources_number doa
              model_order_optimal_param.method = model_order_method;
              [~,~,log_penalty_cost_opt] = model_order_optimal(model_order_optimal_param);
            end
            
            % For each run, collect the eigenvalues and
            % penalty, log-likelihood, and cost (sum of
            % liglikelihood and penalty) from optimal and
            % suboptimal methods.
            param_debug.eigval_all(Run_Idx,:,:) = eigvals; % Nruns*Nc*Nb
            param_debug.opt{Run_Idx}            = log_penalty_cost_opt; %[log_fun (1 by M+1), penalty (1 by M+1), cost (1 by M+1)]
            param_debug.subopt{Run_Idx}         = log_penalty_cost_subopt;
            
          end
          
          % For each training SNR, for the case of Q=0, collect the loglikelihood
          % values from all runs.Each cell array is Nruns by
          % number of targets (M+1).
          for run_idx = 1: param.monte.runs
            opt_log_func_all{SNR_testing_idx}(run_idx,:)    =  param_debug.opt{run_idx}([0:M]+1);
            subopt_log_func_all{SNR_testing_idx}(run_idx,:) =  param_debug.subopt{run_idx}([0:M]+1);
          end
        end
        
        % ======== Now normalize loglikelihood values ========
        % Case 1 normalize the optimal loglikelihood values wrt
        % suboptimal loglikelihood values when Q=0 (SNR=-inf).
        opt_norm_term        = mean(subopt_log_func_all{1},1) - mean(opt_log_func_all{1},1);
        
        % Case 2 (the best case) normalize the suboptimal loglikelihood
        % values wrt suboptimal loglikelihood values when Q=0 (SNR=-inf).
        norm_term_suboptimal = - mean(subopt_log_func_all{1},1);
        
        % Case 3 (the best case) normalize the optimal loglikelihood
        % values wrt optimal loglikelihood values when Q=0 (SNR=-inf).
        norm_term_optimal    = - mean(opt_log_func_all{1},1);
        
        % Collect logliklihood values and DoAs for Q=0 case from all runs.
        param_debug_Q_0.eigval_all{SNR_testing_idx} = param_debug.eigval_all;
        param_debug_Q_0.opt{SNR_testing_idx}        = opt_log_func_all{SNR_testing_idx};
        param_debug_Q_0.subopt{SNR_testing_idx}     = subopt_log_func_all{SNR_testing_idx};
        param_debug_Q_0.doa_mle{SNR_testing_idx}    = doa_mle_all;
      end
      clear opt_log_func_all subopt_log_func_all
      param_MOE.opt_norm_term        = opt_norm_term;
      param_MOE.norm_term_suboptimal = norm_term_suboptimal;
      param_MOE.norm_term_optimal    = norm_term_optimal;
    else
      param_MOE.opt_norm_term        = zeros(1,M+1);
      param_MOE.norm_term_suboptimal = zeros(1,M+1);
      param_MOE.norm_term_optimal    = zeros(1,M+1);
    end
  else
    % write already saved values
    load('N_101_1D.mat','param_MOE')
    %load('N_11_all_SNR_NB_1D.mat','param_MOE')
  end
  opt_norm_term        = param_MOE.opt_norm_term ;
  norm_term_suboptimal = param_MOE.norm_term_suboptimal ;
  norm_term_optimal    = param_MOE.norm_term_optimal ;
  
  param_debug_Q_0.opt_norm_term        = param_MOE.opt_norm_term;
  param_debug_Q_0.norm_term_suboptimal = param_MOE.norm_term_suboptimal;
  param_debug_Q_0.norm_term_optimal    = param_MOE.norm_term_optimal;
  
  % ---------------------------------------------------------------------
  %% TRAINING THE OPTIMIZER (Numerical Tuning of penalty term)
  % ---------------------------------------------------------------------
  % The main aprameters here are:
  % 1) param_debug_NT.opt/subopt: contains the loglikelihoods from all runs,
  % but for a single training SNR for Q=0 case.
  % 2) param_MOE: defined in the previous section.
  param.param_MOE = param_MOE;
  
  % SNR CASES FOR OPTIMIZER  (TRAINING DATA)
  param.SNR_db = [10 20 30] ;
  %param.SNR_db = [30] ;
  if penalty_saved==0
    if D1_optimizer
      % Generate and save optimizer results
      data_idx = 1;
      data     = [];
      
      % param.opt for accessing doa_example_NT doa_example_NT_opt called
      % using same optimizer_NT
      if suboptimal_test && suboptimal_methods(1)==0  %Numerical tuning is required
        param.opt = 0;
        [penalty_sravya, ~, param_debug_NT_subopt] = optimizer_NT(param) ;
      else
        penalty_sravya = zeros(1,Nc);
      end
      
      if optimal_test &&optimal_methods(1)==0
        param.opt = 1;
        [penalty_sravya_opt, ~, param_debug_NT_opt] = optimizer_NT(param);
      else
        penalty_sravya_opt = zeros(1,Nc);
      end
      % penalty_sravya_opt = [[3.17810000000000,334.428099999999,668.856200000000,1003.28430000000,1337.71240000000,1672.14050000000,3171.11508971558]]
    else
      penalty_sravya     = zeros(1,Nc);
      penalty_sravya_opt = zeros(1,Nc);
    end
    
    % param_debug_NT_subopt/opt is a cell of M elements. It has two
    % fields log_func_all (for each SNR) and penalty_min_jump_SNR (Don't know why
    % Sravya needed this second field).
    param_debug_NT.subopt = param_debug_NT_subopt;
    param_debug_NT.opt    = param_debug_NT_opt;
  else
    % Load an already saved values. For example
    load('N_101_1D.mat','penalty_sravya','penalty_sravya_opt')
    %  load('N_11_all_SNR_NB_1D.mat','penalty_sravya','penalty_sravya_opt')
  end
  
  % ---------------------------------------------------------------------
  %% TESTING
  % ---------------------------------------------------------------------
  % param_debug.opt/subopt used in this section for storing the
  % log-likelihoods for each run in the testing mode.
  param_debug = [];
  param.monte.random_seed_offset = 2;
  SNR_testing = param.SNR_db;  % SNR TEST CASES
  % SNR_testing = [3 15 40];   % SNR EDGE CASES TESTING
  
  AIC_RESULT_mean   = zeros(M+1,1);
  MDL_RESULT_mean   = zeros(M+1,1);
  HQ_RESULT_mean    = zeros(M+1,1);
  AICc_RESULT_mean  = zeros(M+1,1);
  KICvc_RESULT_mean = zeros(M+1,1);
  WIC_RESULT_mean   = zeros(M+1,1);
  NT_RESULT_mean    = zeros(M+1,1);
  
  AIC_RESULT_mean_opt   = zeros(M+1,1);
  MDL_RESULT_mean_opt   = zeros(M+1,1);
  HQ_RESULT_mean_opt    = zeros(M+1,1);
  AICc_RESULT_mean_opt  = zeros(M+1,1);
  KICvc_RESULT_mean_opt = zeros(M+1,1);
  WIC_RESULT_mean_opt   = zeros(M+1,1);
  NT_RESULT_mean_opt    = zeros(M+1,1);
  
  % lambda: wavelength
  lambda = c/fc;
  % k: wavenumber
  k = 2*pi/(lambda/2);
  p = Nc;
  d= lambda/2;
  L =p*d;
  
  k_spacing = [-(p-1)/2:1:(p-1)/2]*(lambda/L); % as implemented in paper for orthogonal sv
  %k_spacing = [-p/2:1:(p/2)-1]*(lambda/L); % John
  DOA_orthogonal = asin(k_spacing)*180/pi;
  
  index_ref = ceil(p/2);
  clear eigenvalues_runs_snr clear eigenvalues_mean
  for SNR_testing_idx =1: length(SNR_testing)
    q_idx= 0;
    for q = 0:M
      target_idx = q+1;   % TO SAVE ALL MODEL ORDERS
      q_idx = q_idx +1; %USED TO SAVE PARAMS DEBUG
      if q == 0
        % angs_elecs: 1 x N matrix, angle of arrival for each source (deg)
        DOA_true = [];
      else
        if mod(q,2)==1 %odd
          index =((index_ref)-(q-1)/2):1:(index_ref)+(q-1)/2;
        else
          index =((index_ref)-(q)/2):1:(index_ref)+(q)/2;
          index(find(index==index_ref))= [];
        end
        DOA_true = DOA_orthogonal(index);
        clear index
      end
      %             if isempty(DOA_true)==0   %% SRAVYA ==> this is wrong
      if isempty(DOA_true)   %% Mohanad
        param.monte.SNR   = -inf ;
        param.monte.DOA   = [];
      else
        param.monte.SNR   = repmat(SNR_testing(SNR_testing_idx ),1, length(DOA_true));  %% SNR 10dB used in the CHEN paper
        param.monte.DOA   = repmat(DOA_true,[num_tests 1]);
      end
      % M : Maximum number of signals
      % possible range of k (number of signals)
      if suboptimal_test ==1
        k=1 ;   % FOR Suboptimal COMPARISION WE DO NOT NEED DOA SO ONE SIMULATION RESULTS ENOUGH FROM doa WHICH IS DCM
      end
      if optimal_test ==1
        k=1:M ;
      end
      
      for Nsig_tmp = k
        % k is 1:M WHEN WE WANTED TO RUN OPTIMAL METHODS
        % for each k we should run doa.m whereas only suboptimal methods
        % are required just k= 1 is enough to generate data and as we do not
        % need doa estimation(no need to run doa.m multiple times)
        param.Nsig_tmp = Nsig_tmp;
        %% Run the simulation
        % This section used in doa_example1
        
        %   if Nsig_tmp == 1
        %   [results{Nsig_tmp}, DCM_runs] = sim.doa(param,[]);
        %   elseif Nsig_tmp > 1
        %   doa_prev = squeeze(results{1,Nsig_tmp-1}.theta_est{1,param.method.list});
        %   [results{Nsig_tmp}, DCM_runs] = sim.doa(param, doa_prev(:,1:Nsig_tmp-1));
        %   end
        %
        clear doa
        [results{Nsig_tmp}, DCM_runs] = doa(param);
        
%         DCM_runs_q{SNR_testing_idx}{q_idx} = DCM_runs;
      
        % Stor the DCM and eigenvalues from each run and SNR (not needed
        % for NT MOE, but for any other purpose that might come out in the
        % future.
        DCM_runs_snr{SNR_testing_idx} = DCM_runs; % {snr_idx}{run_idx}
          
        for DCM_idx = 1:param.monte.runs
          for nb_idx = 1:Nb
            clear Rxx
            % DCM of this subband
            Rxx = DCM_runs{DCM_idx};
            Rxx = Rxx(Nc*(nb_idx-1)+1:nb_idx*Nc,:);
            
            eigenvalues_runs_snr{SNR_testing_idx}(DCM_idx,q_idx,nb_idx,:) = abs(real(sort(eig(Rxx),'descend'))); % Nruns*Ntargets*Nb*Nc
          end
        end
         % Mean eigenvalues over all runs
         for nb_idx = 1:Nb
          eigenvalues_mean(SNR_testing_idx,q_idx,nb_idx,:) = mean(squeeze(eigenvalues_runs_snr{SNR_testing_idx}(:,q_idx,nb_idx,:)),1); %Ntests*Ntargets*Nb*Nc
         end
         
        if optimal_test ==1   % 1 for optimal
          doa_tmp(1:param.monte.runs,:) = squeeze(results{1, Nsig_tmp}.theta_est{1, param.method.list});
          doa_mle_all_k{Nsig_tmp}       = doa_tmp(:,1:Nsig_tmp); % each cell contain DOA for all runs. row indicates the run number.
        end
      end
      %       if optimal_test ==1
      %       param_debug_all.doa_mle{SNR_testing_idx}{q_idx} =  doa_mle_all_k;
      %       end
      %
      Nruns = param.monte.runs;  % for now
      if 1
        for Run_Idx = 1:Nruns;
          if optimal_test ==1   % 1 for optimal
            %arranging all doa for possible k for a single run.
            for Nsig_tmp = 1:M  %%%%%CHANGE
              doa_mle{Nsig_tmp} = doa_mle_all_k{Nsig_tmp}(Run_Idx,:) ;
            end
          end
          
          clear eigvals eigvecs
          for nb_idx = 1:Nb
            clear Rxx
            % DCM of this subband
            Rxx = DCM_runs{Run_Idx};
            Rxx = Rxx(Nc*(nb_idx-1)+1:nb_idx*Nc,:);
            [eigval,index] = sort(eig(Rxx),'descend');
            
            % DCM is symmetric and always we get real eigen values. due to some
            % rounding errors in the data generated by matlab we got complex eigen
            % values
            
            eigval = abs(real(eigval));
%             eigval_all(Run_Idx,:) = eigval;
            [V,D] = eig(Rxx);
            eigvec = V(:,index);
            
            eigvals(:,nb_idx) = eigval;
            eigvecs(Nc*(nb_idx-1)+1:nb_idx*Nc,:) = eigvec;
          end
          %% Suboptimal Methods
          
          %  0 IF SUBOPTIMAL ESTIMATORS ARE NOT NEEDED
          if suboptimal_test
            model_order_suboptimal_param.Nc         = Nc;
            model_order_suboptimal_param.Nsnap      = param.monte.Nsnap;
            model_order_suboptimal_param.eigval     = eigvals;
            model_order_suboptimal_param.penalty_NT = penalty_sravya;%zeros(1,M+1);
            model_order_suboptimal_param.param_MOE  = param_MOE;
            model_order_suboptimal_param.filter_banks = param.method.wb_fd.filter_banks;
            
            for model_order_method = suboptimal_methods
              clear sources_number doa
              
              model_order_suboptimal_param.method = model_order_method;
              [sources_number, log_penalty_cost_subopt] = model_order_suboptimal(model_order_suboptimal_param);
              
              param_debug.eigval_all(Run_Idx,:,:) = eigvals; % Nruns*Nc*Nb
              param_debug.subopt{Run_Idx}(model_order_method+1,:) = log_penalty_cost_subopt;
              
              if model_order_method == 0
                loglikelihood_runs_snr.subopt{SNR_testing_idx}(Run_Idx,q_idx,:) = log_penalty_cost_subopt(1:M+1);
              end
              
              doa = NaN *ones(1,Nc-1);
              
              if optimal_test ==1   % when optimal also implemented since then we have doa
                if sources_number > 0
                  doa(1:sources_number) = doa_mle{sources_number};
                end
                [doa,sort_idxs] = sort(doa);
              end
              
              switch model_order_method
                case 0
                  model_order_results_suboptimal.NT.doa(Run_Idx,:,target_idx)  = doa;
                  model_order_results_suboptimal.NT.Nest(Run_Idx,target_idx)  = sources_number;
                case 1
                  model_order_results_suboptimal.AIC.doa(Run_Idx,:,target_idx)  = doa;
                  model_order_results_suboptimal.AIC.Nest(Run_Idx,target_idx)  = sources_number;
                case 2
                  model_order_results_suboptimal.HQ.doa(Run_Idx,:,target_idx)  = doa;
                  model_order_results_suboptimal.HQ.Nest(Run_Idx,target_idx)  = sources_number;
                case 3
                  model_order_results_suboptimal.MDL.doa(Run_Idx,:,target_idx)  = doa;
                  model_order_results_suboptimal.MDL.Nest(Run_Idx,target_idx)  = sources_number;
                case 4
                  model_order_results_suboptimal.AICc.doa(Run_Idx,:,target_idx)  = doa;
                  model_order_results_suboptimal.AICc.Nest(Run_Idx,target_idx)  = sources_number;
                case 5
                  model_order_results_suboptimal.KICvc.doa(Run_Idx,:,target_idx)  = doa;
                  model_order_results_suboptimal.KICvc.Nest(Run_Idx,target_idx)  = sources_number;
                case 6
                  model_order_results_suboptimal.WIC.doa(Run_Idx,:,target_idx)  = doa;
                  model_order_results_suboptimal.WIC.Nest(Run_Idx,target_idx)  = sources_number;
                otherwise
                  error('Not supported')
              end
            end
          end
          
          %% Model Order Estimators
          if optimal_test ==1
            % Optimal Methods
            phase_center = param.src.lever_arm.fh(param.src.lever_arm.args{:});
            model_order_optimal_param.y_pc       = phase_center(2,:).';
            model_order_optimal_param.z_pc       = phase_center(3,:).';
            model_order_optimal_param.fc         = fc;
            model_order_optimal_param.Nc         = Nc;
            model_order_optimal_param.Nsnap      = param.monte.Nsnap;
            model_order_optimal_param.eigval     = eigvals;
            model_order_optimal_param.eigvec     = eigvecs;
            model_order_optimal_param.penalty_NT_opt = penalty_sravya_opt;%zeros(1,M+1);
            model_order_optimal_param.param_MOE  = param_MOE;
            model_order_optimal_param.doa_mle    = doa_mle;
            model_order_optimal_param.filter_banks = param.method.wb_fd.filter_banks;
            model_order_optimal_param.fs = BW;
            
            for model_order_method = 0:6
              clear sources_number doa
              
              model_order_optimal_param.method = model_order_method;
              
              % log_penalty_cost_opt is 3*M vector for each method.
              % That is,[log_func(1*M) penalty(1*M) cost(1*M)]
              [sources_number,doa,log_penalty_cost_opt] = model_order_optimal(model_order_optimal_param);
              
              param_debug.opt{Run_Idx}(model_order_method+1,:) = log_penalty_cost_opt;
              
              if model_order_method == 0
                loglikelihood_runs_snr.opt{SNR_testing_idx}(Run_Idx,q_idx,:) = log_penalty_cost_opt(1:M+1);
              end
              
              [doa,sort_idxs] = sort(doa);
              
              switch model_order_method
                case 0
                  model_order_results_optimal.NT.doa(Run_Idx,:,target_idx)  = doa;
                  model_order_results_optimal.NT.Nest(Run_Idx,target_idx)  = sources_number;
                case 1
                  model_order_results_optimal.AIC.doa(Run_Idx,:,target_idx)  = doa;
                  model_order_results_optimal.AIC.Nest(Run_Idx,target_idx)  = sources_number;
                case 2
                  model_order_results_optimal.HQ.doa(Run_Idx,:,target_idx)  = doa;
                  model_order_results_optimal.HQ.Nest(Run_Idx,target_idx)  = sources_number;
                case 3
                  model_order_results_optimal.MDL.doa(Run_Idx,:,target_idx)  = doa;
                  model_order_results_optimal.MDL.Nest(Run_Idx,target_idx)  = sources_number;
                case 4
                  model_order_results_optimal.AICc.doa(Run_Idx,:,target_idx) = doa;
                  model_order_results_optimal.AICc.Nest(Run_Idx,target_idx)  = sources_number;
                case 5
                  model_order_results_optimal.KICvc.doa(Run_Idx,:,target_idx) = doa;
                  model_order_results_optimal.KICvc.Nest(Run_Idx,target_idx)  = sources_number;
                case 6
                  model_order_results_optimal.WIC.doa(Run_Idx,:,target_idx) = doa;
                  model_order_results_optimal.WIC.Nest(Run_Idx,target_idx)  = sources_number;
                otherwise
                  error('Not supported')
              end
            end
          end
        end
        dout.model_order_results_optimal    = model_order_results_optimal;
        dout.model_order_results_suboptimal = model_order_results_suboptimal;
        dout_temp{SNR_testing_idx,1} = dout;
        
        clear dout
        
        %% Saving results from all model order estimators for comparision
        if stat_test==1
          %  0 IF SUBOPTIMAL ESTIMATORS ARE NOT NEEDED
          if suboptimal_test == 1     % suboptimal
            %Nest_results(:,1) = length(DOA_true).*ones(bin,1);
            Nest_results(:,1) = model_order_results_suboptimal.AIC.Nest(:,target_idx);
            Nest_results(:,2) = model_order_results_suboptimal.HQ.Nest(:,target_idx);
            Nest_results(:,3) = model_order_results_suboptimal.MDL.Nest(:,target_idx);
            Nest_results(:,4) = model_order_results_suboptimal.AICc.Nest(:,target_idx);
            Nest_results(:,5) = model_order_results_suboptimal.KICvc.Nest(:,target_idx);
            Nest_results(:,6) = model_order_results_suboptimal.WIC.Nest(:,target_idx);
            Nest_results(:,7) = model_order_results_suboptimal.NT.Nest(:,target_idx);
          else
            Nest_results = zeros(Nruns,7);
          end
          
          if optimal_test == 1 %optimal
            Nest_results(:,8) = model_order_results_optimal.AIC.Nest(:,target_idx);
            Nest_results(:,9) = model_order_results_optimal.HQ.Nest(:,target_idx);
            Nest_results(:,10) = model_order_results_optimal.MDL.Nest(:,target_idx);
            Nest_results(:,11) = model_order_results_optimal.AICc.Nest(:,target_idx);
            Nest_results(:,12) = model_order_results_optimal.KICvc.Nest(:,target_idx);
            Nest_results(:,13) = model_order_results_optimal.WIC.Nest(:,target_idx);
            Nest_results(:,14) = model_order_results_optimal.NT.Nest(:,target_idx);
          end
          
        else  % stat_test ==0  (testing to save params_debug_all)
          
          Nest_results(:,1) = model_order_results_suboptimal.NT.Nest(:,target_idx);
          Nest_results(:,2) = model_order_results.NT.Nest(:,target_idx);
        end
        
        methods_num = size(Nest_results,2);
        
        for k =0:M
          for i= 1:methods_num
            Statistics{i}(q+1,k+1) =  (length(find(Nest_results(:,i)==k))/Nruns)*100 ;
          end
        end
        
        % -------------------------------------------------------------------------
        %% SAVE PARAMS DEBUG for each Q
        % -------------------------------------------------------------------------
        param_debug_all.eigval_all{SNR_testing_idx}{q_idx} = param_debug.eigval_all; % Nruns*Nc*Nb
        
        if optimal_test ==1
          param_debug_all.opt{SNR_testing_idx}{q_idx} = param_debug.opt;
          
          for run_idx = 1: param.monte.runs
            opt_log_func_all{SNR_testing_idx}{q_idx}(run_idx,:) =  param_debug.opt{run_idx}([0:M]+1);
          end
          param_debug_all.doa_mle{SNR_testing_idx}{q_idx} =  doa_mle_all_k;
          param_debug.opt = [];
        end
        
        if suboptimal_test ==1
          param_debug_all.subopt{SNR_testing_idx}{q_idx} = param_debug.subopt;
          for run_idx = 1: param.monte.runs
            subopt_log_func_all{SNR_testing_idx}{q_idx}(run_idx,:) =  param_debug.subopt{run_idx}(1,[0:M]+1);
          end
          param_debug.subopt = [];
        end
        param_debug.eigval_all = [];
      end
    end
  end
  
  % dout_tmp is a cell array of length=number of testing SNRs.Each cell
  % array contains the MOE results (Nest and doa) for optimal and suboptimal methods
  % (all 7 methods) ==>dout_temp{1}.model_order_results_optimal.AIC.Nest
  dout = dout_temp;
  SNR_training = param.SNR_db;
  
  % -------------------------------------------------------------------------
  %% Formatting the results
  % -------------------------------------------------------------------------
  for SNR_idx = 1:length(SNR_testing)
    
    subopt_log_func_all_Q = [];
    opt_log_func_all_Q    = [];
    sources_true_Q        = [];
    
    for idx_shape = 1: size(subopt_log_func_all{SNR_idx},2)
      if ~isempty(subopt_log_func_all{SNR_idx}{idx_shape})
        subopt_log_func_all_Q = vertcat(subopt_log_func_all_Q,subopt_log_func_all{SNR_idx}{idx_shape});
        opt_log_func_all_Q    = vertcat(opt_log_func_all_Q,opt_log_func_all{SNR_idx}{idx_shape});
        sources_true_Q        = vertcat(sources_true_Q, (idx_shape-1)*ones(size(subopt_log_func_all{SNR_idx}{idx_shape},1),1));
      end
    end
    param_debug_testing.subopt{SNR_idx}       = subopt_log_func_all_Q;
    param_debug_testing.opt{SNR_idx}          = opt_log_func_all_Q;
    param_debug_testing.sources_true{SNR_idx} = sources_true_Q;
  end
  
  if D1_optimizer==1
    
    for SNR_idx = 1:length(param.SNR_db)
      
      subopt_log_func_all_Q = [];
      opt_log_func_all_Q   = [];
      sources_true_Q       = [];
      
      % Loop over all targets
      for idx_shape = 1: length(param_debug_NT.subopt)
        if ~isempty(param_debug_NT.subopt{idx_shape})
          subopt_log_func_all_Q = vertcat(subopt_log_func_all_Q,param_debug_NT.subopt{idx_shape}.log_func_all{SNR_idx});
          opt_log_func_all_Q    = vertcat(opt_log_func_all_Q,param_debug_NT.opt{idx_shape}.log_func_all{SNR_idx});
          sources_true_Q        = vertcat(sources_true_Q, (idx_shape-1)*ones(size(param_debug_NT.subopt{idx_shape}.log_func_all{SNR_idx},1),1));
        end
      end
      param_debug_nt.subopt{SNR_idx,1} = subopt_log_func_all_Q;
      param_debug_nt.opt{SNR_idx,1} = opt_log_func_all_Q;
      
      param_debug_nt.sources_true{SNR_idx,1} = sources_true_Q;
    end
    
    clear param_debug_NT
    param_debug_NT = param_debug_nt;
    
  end
  
  % Remove NaN from the log-likelihoods
  if 1
    for SNR_idx = 1:length(SNR_testing)
      for idx= [0:M]+1
        test = isnan(param_debug_NT.opt{SNR_idx, 1}(:,idx));
        row_idx{idx} = find(test==1);
        %now we can ignore rows with row_idx as they have NaN cases
        param_debug_NT.opt{SNR_idx, 1}(row_idx{idx},:)          = [];
        param_debug_NT.subopt{SNR_idx, 1}(row_idx{idx},:)       = [];
        param_debug_NT.sources_true{SNR_idx, 1}(row_idx{idx},:) = [];
      end
    end
  end
  
  %% Store DCMs and eigenvalues from all runs and SNRs, as well as simulation parameters 
  if save_eigenvalues
    sim_param.fc           = fc;
    sim_param.BW           = BW;
    sim_param.Nc           = Nc;
    sim_param.Nsnap        = param.monte.Nsnap;
    sim_param.M            = M;
    sim_param.SNR_training = param.SNR_db;
    sim_param.phase_center = phase_center;
    sim_param.Nb = Nb;
    sim_param.Nruns = Nruns;
    sim_param.notes{1} ='Each eigenvalues cell correspond to one test case. The eigenvalues of each test case are arrays with dimension Nruns by Ntargets by Nb by Nc';
    
    penalty.opt    = penalty_sravya_opt;
    penalty.subopt = penalty_sravya;
    
    normalization.norm_allign_zero     = param_MOE.norm_allign_zero;
    normalization.norm_term_optimal    = param_MOE.norm_term_optimal;
    normalization.norm_term_suboptimal = param_MOE.norm_term_suboptimal;
    normalization.opt_norm_term        = param_MOE.opt_norm_term;
    
%    loglikelihood_runs_snr.opt    = param_debug_NT.opt;
%    loglikelihood_runs_snr.subopt = param_debug_NT.subopt;
   sources_true = param_debug_NT.sources_true;
    out_fn_dir = '/users/mohanad/IceSheetProject/MOE work/DCMandEigenvalues/';
    out_fn_name = '1D';
    if ~exist(out_fn_dir,'dir')
      mkdir(out_fn_dir);
    end
    out_fn = fullfile(out_fn_dir,out_fn_name);
    save([out_fn '.mat'],'DCM_runs_snr','eigenvalues_runs_snr','loglikelihood_runs_snr','sim_param','penalty','normalization')
  end
  
  
  %%                      Plot Results
  %---------------------------------------------------------------------
  %   eigenvalues_mean(SNR_testing_idx,q_idx,nb_idx,:)
  if 1
    % Plot eigenvalues
    for nb_idx = 1:Nb
      % eigenvalues_mean is Ntests*Ntargets*Nb*Nc matrix
      eigenvalues_mean_nb = squeeze(eigenvalues_mean(:,:,nb_idx,:));
      figure(1000+nb_idx);clf
      for q_idx = 1:size(eigenvalues_mean_nb,2)
        % Eigenvalues for all the SNRs for a given Nsig
        eigenvalue_q = squeeze(eigenvalues_mean_nb(:,q_idx,:)).';
        eigenvalue_q = 10*log10(eigenvalue_q);%./repmat(max(eigenvalue_q,[],1),[size(eigenvalue_q,1) 1]));
        
        subplot(ceil(size(eigenvalues_mean_nb,2)/2),2,q_idx)
        plot(repmat([0:1:Nc-1].',[1 3]),eigenvalue_q,'-*')
        title(sprintf('q = %d',q_idx-1))
        
        if q_idx == (size(eigenvalues_mean_nb,2))-1 || q_idx == size(eigenvalues_mean_nb,2)
          xlabel('Eigenvalue index')
        end
        
        if q_idx == 1 %mod(q_idx,2) == 1
          ylabel('Eigenvalue (dB)')
        end
        grid on
        if q_idx == 1
          legend('10dB','20dB','30dB','Location','southwest')
        end
        xlim([0 Nc-1])
        ylim([min(eigenvalue_q(:))   max(eigenvalue_q(:))])
      end
      
      if Nb>1
        suptitle(sprintf('1D simulator: Bf = %2.2f percent, subband# %1.0d',(BW/Nb)/fc *100,nb_idx))
      else
        suptitle(sprintf('1D simulator: Bf = %2.2f percent',BW/fc *100))
      end
    end
    
    % Mean eigenvalues over all Nb subbands
    if Nb > 1
      eivenvalues_mean_nb_mean = squeeze(mean(eigenvalues_mean,3));
      for nb_idx = 1:Nb
        figure(1000+Nb+1);clf
        for q_idx = 1:size(eivenvalues_mean_nb_mean,2)
          % Eigenvalues for all the SNRs for a given Nsig
          eigenvalue_q = squeeze(eivenvalues_mean_nb_mean(:,q_idx,:)).';
          eigenvalue_q = 10*log10(eigenvalue_q);%./repmat(max(eigenvalue_q,[],1),[size(eigenvalue_q,1) 1]));
          
          subplot(ceil(size(eivenvalues_mean_nb_mean,2)/2),2,q_idx)
          plot(repmat([0:1:Nc-1].',[1 3]),eigenvalue_q,'-*')
          title(sprintf('q = %d',q_idx-1))
          
          if q_idx == (size(eivenvalues_mean_nb_mean,2))-1 || q_idx == size(eivenvalues_mean_nb_mean,2)
            xlabel('Eigenvalue index')
          end
          
          if q_idx == 1 %mod(q_idx,2) == 1
            ylabel('Eigenvalue (dB)')
          end
          grid on
          if q_idx == 1
            legend('10dB','20dB','30dB','Location','southwest')
          end
          xlim([0 Nc-1])
          ylim([min(eigenvalue_q(:))   max(eigenvalue_q(:))])
        end
           suptitle(sprintf('1D simulator: Bf = %2.2f percent, mean over %1.0d subbands',(BW/Nb)/fc *100,Nb))
      end      
    end
  end
  
%   return
  
  if likelihood_plots == 1
    %% Log likelihood plots
    %         warning('off','MATLAB:legend:IgnoringExtraEntries')
    
    if opt_norm == 1
      for SNR_idx = 1:length(SNR_training_Q_0)
        %TESTING DATA USED TO CALCULATE Q=0 case
        figure(3);clf; subplot(length(SNR_training_Q_0),1,SNR_idx)
        plot(0:M,mean(param_debug_Q_0.opt{SNR_idx,1},1),'b-*' )
        hold on
        plot(0:M,mean(param_debug_Q_0.subopt{SNR_idx,1},1),'r--*' )
        
        xlabel('k','interpreter','none')
        ylabel('-2L','interpreter','none')
        %       set(gca,'TickLabelInterpreter','Latex')
        h_legend = legend('Suboptimal','Optimal','Location','best');
        set(h_legend,'Interpreter','Latex')
        title([  num2str(SNR_training_Q_0(SNR_idx))'' ' dB    Q=0 '],'interpreter','Latex'  )
        grid on
      end
    end
    
    if D1_optimizer==1
      %TRAINING DATA
      figure(1);clf
      figure(2);clf
      
      Color         = {'b','r','k','c','m','g',[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560]};
      Marker_opt    = {'-*','-o','-+','-.','-s','-^','-<','->'};
      Marker_subopt = {'--*','--o','--+','--.','--s','--^','--<','-->'};
      for SNR_idx = 1: length(SNR_training)
        % Include this section later (dealing NaN cases)
%         if 1
%           for idx= [0:M]+1
%             test = isnan(param_debug_NT.opt{SNR_idx, 1}(:,idx));
%             row_idx{idx} = find(test==1);
%             %now we can ignore rows with row_idx as they have NaN cases
%             param_debug_NT.opt{SNR_idx, 1}(row_idx{idx},:)          = [];
%             param_debug_NT.subopt{SNR_idx, 1}(row_idx{idx},:)       = [];
%             param_debug_NT.sources_true{SNR_idx, 1}(row_idx{idx},:) = [];
%           end
%         end
        
        % PLOT LOG-LIKELIHOOD (NOT COST) WITH NORMALIZATION:
        % TRAINING
        %--------------------------------------------------
        % Plot optimal methods
        figure(1);
        subplot(length(SNR_training),2,2*SNR_idx-1)
        for idx = 0:M
          plot(0:M,mean(param_debug_NT.opt{SNR_idx,1}(find(param_debug_NT.sources_true{SNR_idx,1}==idx),:),1),Marker_opt{idx+1},'Color',Color{idx+1} )
          hold on,
        end
        grid on
        
        % Plot suboptimal methods
        figure(2)
        subplot(length(SNR_training),2,2*SNR_idx-1)
        for idx = 0:M
          plot(0:M,mean(param_debug_NT.subopt{SNR_idx,1}(find(param_debug_NT.sources_true{SNR_idx,1}==idx),:),1),Marker_subopt{idx+1},'Color',Color{idx+1} )
          hold on
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
        for idx = 0:M
          if param_MOE.norm_allign_zero ==1
            plot(0:M,mean(param_debug_NT.opt{SNR_idx,1}(find(param_debug_NT.sources_true{SNR_idx,1}==idx),:)-repmat(norm_term_optimal,length(find(param_debug_NT.sources_true{SNR_idx,1}==idx)),1),1),Marker_opt{idx+1},'Color',Color{idx+1})
          else
            plot(0:M,mean(param_debug_NT.opt{SNR_idx,1}(find(param_debug_NT.sources_true{SNR_idx,1}==idx),:)-repmat(opt_norm_term,length(find(param_debug_NT.sources_true{SNR_idx,1}==idx)),1),1),Marker_subopt{idx+1},'Color',Color{idx+1})
          end
          hold on
        end
        grid on
        
        % Plot suboptimal methods
        figure(2);
        subplot(length(SNR_training),2,2*SNR_idx)
        for idx = 0:M
          if param_MOE.norm_allign_zero ==1
            plot(0:M,mean(param_debug_NT.subopt{SNR_idx,1}(find(param_debug_NT.sources_true{SNR_idx,1}==idx),:)-repmat(norm_term_suboptimal,length(find(param_debug_NT.sources_true{SNR_idx,1}==idx)),1),1),Marker_subopt{idx+1},'Color',Color{idx+1} )
          else
            plot(0:M,mean(param_debug_NT.subopt{SNR_idx,1}(find(param_debug_NT.sources_true{SNR_idx,1}==idx),:),1),Marker_subopt{idx+1},'Color',Color{idx+1} )
          end
          hold on,
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
  
  %% HISTOGRAM FOR A METHOD W.R.T MODEL ORDER
  for SNR_idx =1:length(SNR_testing)
    
    Nest_suboptimal= dout{SNR_idx,1}.model_order_results_suboptimal;
    Nest_optimal = dout{SNR_idx,1}.model_order_results_optimal ;
    
    RLINES = size(Nest_suboptimal.AIC.Nest,2);   %ALTER THIS
    
    for method_idx = methods
      switch method_idx
        case 0
          Nest_subopt = Nest_suboptimal.NT.Nest;
          Nest_opt = Nest_optimal.NT.Nest;
        case 1
          
          Nest_subopt = Nest_suboptimal.AIC.Nest;
          Nest_opt = Nest_optimal.AIC.Nest;
        case 2
          
          Nest_subopt = Nest_suboptimal.HQ.Nest;
          Nest_opt = Nest_optimal.HQ.Nest;
        case 3
          
          Nest_subopt = Nest_suboptimal.MDL.Nest;
          Nest_opt = Nest_optimal.MDL.Nest;
        case 4
          
          Nest_subopt = Nest_suboptimal.AICc.Nest;
          Nest_opt = Nest_optimal.AICc.Nest;
        case 5
          
          Nest_subopt = Nest_suboptimal.KICvc.Nest;
          Nest_opt = Nest_optimal.KICvc.Nest;
        case 6
          
          Nest_subopt = Nest_suboptimal.WIC.Nest;
          Nest_opt = Nest_optimal.WIC.Nest;
        case 7
          
          Nest_subopt = Nest_suboptimal.AICc_19.Nest;
          Nest_opt = Nest_optimal.AICc_19.Nest;
          
        otherwise
          error('Not supported')
      end
      
      Nest_subopt = reshape(Nest_subopt,[numel(Nest_subopt),1]); % (Nruns*M)-by-1
      Nest_opt   = reshape(Nest_opt,[numel(Nest_opt),1]);
      
      for k_idx = 0:M
        if 0
          if numel(find (param_debug_testing.sources_true{SNR_idx,1}(:,1:RLINES) == k_idx)) == 0
            percentage_correct_subopt(k_idx+1) = NaN;
            percentage_correct_opt(k_idx+1) = NaN;
          else
            percentage_correct_subopt(k_idx+1) = ...
              (numel(find(Nest_subopt(find (param_debug_testing.sources_true{SNR_idx,1}(:,1:RLINES) == k_idx)) == k_idx)) / numel(find (param_debug_testing.sources_true{SNR_idx,1}(:,1:RLINES) == k_idx)))*100;
            percentage_correct_opt(k_idx+1) = ...
              (numel(find(Nest_opt(find (param_debug_testing.sources_true{SNR_idx,1}(:,1:RLINES) == k_idx)) == k_idx)) / numel(find (param_debug_testing.sources_true{SNR_idx,1}(:,1:RLINES) == k_idx)))*100;
          end
        end
        
        if numel(find (param_debug_testing.sources_true{SNR_idx} == k_idx)) == 0
          percentage_correct_subopt(k_idx+1) = NaN;
          percentage_correct_opt(k_idx+1) = NaN;
        else
          percentage_correct_subopt(k_idx+1) = ...
            (numel(find(Nest_subopt(find (param_debug_testing.sources_true{SNR_idx} == k_idx)) == k_idx)) / numel(find (param_debug_testing.sources_true{SNR_idx} == k_idx)))*100;
          percentage_correct_opt(k_idx+1) = ...
            (numel(find(Nest_opt(find (param_debug_testing.sources_true{SNR_idx} == k_idx)) == k_idx)) / numel(find (param_debug_testing.sources_true{SNR_idx} == k_idx)))*100;
        end
      end
      
      %performance can't be tested for the below case.
      %since no number of sorces equal the number we are testing in true sources
      %percentage(find(~isnan(percentage) == 1)) = ? ;
      
      percentage_method_subopt(:,method_idx+1) = percentage_correct_subopt; % Ntargets*Nmethods
      percentage_method_opt(:,method_idx+1)    = percentage_correct_opt;
      %
      percentage_method_subopt_all{SNR_idx} = percentage_method_subopt;
      percentage_method_opt_all{SNR_idx}    = percentage_method_opt;
      
      if stat_test ==1
        %STEM
        if 0
          figure(50),subplot(2,2,SNR_idx),hold on,stem([0:M]+(0.03)*(method_idx+1), percentage_method_subopt(:,method_idx+1),'filled')
          xlim([0 M+0.5])
          grid on
          
          grid on
          %     set(gca,'TickLabelInterpreter','Latex')
          if SNR_idx == length(SNR_testing)
            xlabel('Number of sources, q','interpreter','none')
            ylabel('% Correct','interpreter','none')
            %       set(gca,'TickLabelInterpreter','Latex')
            h_legend = legend('NT','AIC', 'HQ', 'MDL', 'AICc', 'KICvc', 'WIC','Location','best');
            set(h_legend,'Interpreter','Latex')
          end
          
          if  SNR_idx ==1
            title([  num2str(SNR_testing(SNR_idx))'' ' dB  Suboptimal'],'interpreter','Latex'  )
          else
            title([  num2str(SNR_testing(SNR_idx))'' ' dB'],'interpreter','Latex'  )
          end
          figure(55),subplot(2,2,SNR_idx),hold on,stem([0:M]+(0.03)*(method_idx+1), percentage_method_opt(:,method_idx+1),'filled')
          xlim([0 M+0.5])
          
          grid on
          %     set(gca,'TickLabelInterpreter','Latex')
          if SNR_idx == length(SNR_testing)
            xlabel('Number of sources, q','interpreter','none')
            ylabel('% Correct','interpreter','none')
            %       set(gca,'TickLabelInterpreter','Latex')
            %         h_legend = legend('NT','AIC', 'HQ', 'MDL', 'AICc', 'KICvc', 'WIC','Location','best');
            %         set(h_legend,'Interpreter','Latex')
          end
          
          if  SNR_idx ==1
            title([  num2str(SNR_testing(SNR_idx))'' ' dB Optimal'],'interpreter','Latex'  )
          else
            title([  num2str(SNR_testing(SNR_idx))'' ' dB'],'interpreter','Latex'  )
          end
        end
        
        Nest_2D_subopt_all(:,method_idx+1) = Nest_subopt;
        Nest_2D_opt_all(:,method_idx+1) = Nest_opt;
      end
    end
    
    if stat_test ==1
      % Plot suboptimal methods MOE results: testing
      figure(4)
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
      grid on
      % set(gca,'TickLabelInterpreter','Latex')
      if SNR_idx == length(SNR_testing)
        xlabel('Number of sources, q','interpreter','none')
        ylabel('% Correct','interpreter','none')
        %   set(gca,'TickLabelInterpreter','Latex')
        
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
      
      % Plot optimal methods MOE results: testing
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
      % set(gca,'TickLabelInterpreter','Latex')
      grid on
      if SNR_idx == length(SNR_testing)
        xlabel('Number of sources, q','interpreter','none')
        ylabel('% Correct','interpreter','none')
        %   set(gca,'TickLabelInterpreter','Latex')
        
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
    
    sources_true_all = param_debug_testing.sources_true;
    
    %reshape may be required
    method_idx = method_idx+1;
    Nest_2D_subopt_all(:,method_idx+1) = sources_true_all{SNR_idx};
    Nest_2D_opt_all(:,method_idx+1) = sources_true_all{SNR_idx}; % (Nruns*Ntargets)-by-Nmethods matrix
    
    %% IMAGESC PLOTS
    if IMAGESC_plots ==1
      if AICc_both==1
        figure,imagesc(Nest_2D_subopt_all(:,[2:5 end-1 6:end-2 1 end]));
        methods_name = {'AIC'  , 'HQ'  , 'MDL'  , 'AICc'  ,'AICc 19' , 'KICvc'  , 'WIC'  , 'NT'  , 'Q'};
      else
        figure(SNR_idx+210);clf
        imagesc([Nest_2D_subopt_all(:,2:end-1),Nest_2D_subopt_all(:,1),Nest_2D_subopt_all(:,end)]);
        methods_name = {'AIC'  , 'HQ'  , 'MDL'  , 'AICc'  , 'KICvc'  , 'WIC'  , 'NT'  , 'Q'};
      end
      
      % set(gca,'TickLabelInterpreter','Latex')
      set(gca,'XtickLabel',methods_name)
      ylabel('Range bins','interpreter','Latex')
      %title('comaparision for a range line')
      cbr = colorbar;
      set(cbr,'YTick',0:1:M)
      title([  num2str(SNR_testing(SNR_idx))'' ' dB   Suboptimal'],'interpreter','Latex'  )
      % title('Suboptimal','interpreter','Latex')
      
      if AICc_both==1
        figure(SNR_idx+110);clf;
        imagesc(Nest_2D_opt_all(:,[2:5 end-1 6:end-2 1 end]));
        methods_name = {'AIC'  , 'HQ'  , 'MDL'  , 'AICc'  ,'AICc 19' , 'KICvc'  , 'WIC'  , 'NT'  , 'Q'};
      else
        figure(SNR_idx+110);clf;
        imagesc([Nest_2D_opt_all(:,2:end-1),Nest_2D_opt_all(:,1),Nest_2D_opt_all(:,end)]);
        methods_name = {'AIC'  , 'HQ'  , 'MDL'  , 'AICc'  , 'KICvc'  , 'WIC'  , 'NT'  , 'Q'};
      end
      
      %   set(gca,'TickLabelInterpreter','Latex')
      set(gca,'XtickLabel',methods_name)
      % Sravya calls the (Nruns*Ntargets) as range-bins to make an
      % analogy with 2D case.
      ylabel('Range bins','interpreter','Latex')
      %title('comaparision for a range line')
      cbr = colorbar;
      set(cbr,'YTick',0:1:M)
      title([  num2str(SNR_testing(SNR_idx))'' ' dB   Optimal'],'interpreter','Latex'  )
      
      %keyboard;
    end
    
    %%
    if likelihood_plots == 1
      
      %TESTING DATA
      
      %% Include this section later (dealing NaN cases)
      if 1
        for idx= [0:M]+1
          test = isnan(param_debug_testing.opt{SNR_idx}(:,idx));
          row_idx{idx} = find(test==1);
          %now we can ignore rows with row_idx as they have NaN cases
          param_debug_testing.opt{SNR_idx}(row_idx{idx},:) = [];
          param_debug_testing.subopt{SNR_idx}(row_idx{idx},:) = [];
          param_debug_testing.sources_true{SNR_idx}(row_idx{idx},:) = [];
        end
      end
      
      Color         = {'b','r','k','c','m','g',[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560]};
      Marker_opt    = {'-*','-o','-+','-.','-s','-^','-<','->'};
      Marker_subopt = {'--*','--o','--+','--.','--s','--^','--<','-->'};
      
      
      % PLOT LOG-LIKELIHOOD (NOT COST) WITH NORMALIZATION:
      % TESTING
      %--------------------------------------------------
      % Plot optimal methods
      figure(6);
      subplot(length(SNR_testing),2,2*SNR_idx-1)
      for idx = 0:M
        plot(0:M,mean(param_debug_testing.opt{SNR_idx}(find(param_debug_testing.sources_true{SNR_idx}==idx),:),1),Marker_opt{idx+1},'Color',Color{idx+1} )
        hold on,
      end
      grid on
      hold off
      
      % Plot suboptimal methods
      figure(7);
      subplot(length(SNR_testing),2,2*SNR_idx-1)
      for idx = 0:M
        plot(0:M,mean(param_debug_testing.subopt{SNR_idx}(find(param_debug_testing.sources_true{SNR_idx}==idx),:),1),Marker_subopt{idx+1},'Color',Color{idx+1} )
        hold on,
      end
      grid on
      hold off
      
      % set(gca,'TickLabelInterpreter','Latex')
      if SNR_idx ==1
        figure(6);
        title([  num2str(SNR_testing(SNR_idx))'' ' dB      --- Optimal  Testing'],'interpreter','Latex'  )
        
        figure(7);
        title([  num2str(SNR_testing(SNR_idx))'' ' dB    -- -- Suboptimal     Testing'],'interpreter','Latex'  )
      end
      
      figure(6)
      ylabel('-2L','interpreter','Latex')
      % set(gca,'TickLabelInterpreter','Latex')
      
      figure(7)
      ylabel('-2L','interpreter','Latex')
      % set(gca,'TickLabelInterpreter','Latex')
      
      if SNR_idx == length(SNR_testing)
        figure(6)
        xlabel('k','interpreter','Latex')
        h_legend = legend('q = 0','q = 1','q = 2','q = 3','q = 4','q = 5','q = 6','Location','best');
        set(h_legend,'Interpreter','Latex')
        
        figure(7)
        xlabel('k','interpreter','Latex')
        h_legend = legend('q = 0','q = 1','q = 2','q = 3','q = 4','q = 5','q = 6','Location','best');
        set(h_legend,'Interpreter','Latex')
      end
      
      if SNR_idx ~= 1
        figure(6)
        title([  num2str(SNR_testing(SNR_idx))'' ' dB'],'interpreter','Latex'  )
        
        figure(7)
        title([  num2str(SNR_testing(SNR_idx))'' ' dB'],'interpreter','Latex'  )
      end
      
      
      grid on
      
      % PLOT LOG-LIKELIHOOD (NOT COST) WITHOUT NORMALIZATION:
      % TESTING
      %--------------------------------------------------
      % Plot optimal methods
      figure(6)
      subplot(length(SNR_testing),2,2*SNR_idx)
      for idx = 0:M
        if param_MOE.norm_allign_zero ==1
          plot(0:M,mean(param_debug_testing.opt{SNR_idx}(find(param_debug_testing.sources_true{SNR_idx}==idx),:)-repmat(norm_term_optimal,length(find(param_debug_testing.sources_true{SNR_idx}==idx)),1),1),Marker_opt{idx+1},'Color',Color{idx+1})
        else
          plot(0:M,mean(param_debug_testing.opt{SNR_idx}(find(param_debug_testing.sources_true{SNR_idx}==idx),:)-repmat(opt_norm_term,length(find(param_debug_testing.sources_true{SNR_idx}==idx)),1),1),Marker_opt{idx+1},'Color',Color{idx+1})
        end
        hold on
      end
      grid on
      
      % Plot suboptimal methods
      figure(7)
      subplot(length(SNR_testing),2,2*SNR_idx)
      for idx = 0:M
        if param_MOE.norm_allign_zero ==1
          plot(0:M,mean(param_debug_testing.subopt{SNR_idx}(find(param_debug_testing.sources_true{SNR_idx}==idx),:)-repmat(norm_term_suboptimal,length(find(param_debug_testing.sources_true{SNR_idx}==idx)),1),1),Marker_subopt{idx+1},'Color',Color{idx+1} )
        else
          plot(0:M,mean(param_debug_testing.subopt{SNR_idx}(find(param_debug_testing.sources_true{SNR_idx}==idx),:),1),Marker_subopt{idx+1},'Color',Color{idx+1})
        end
        hold on
      end
      grid on
      
      if SNR_idx ==1
        figure(6);
        title([  num2str(SNR_testing(SNR_idx))'' ' dB    Without Normalization'],'interpreter','Latex'  )
        
        figure(7);
        title([  num2str(SNR_testing(SNR_idx))'' ' dB    Without Normalization'],'interpreter','Latex'  )
      else
        figure(6);
        title([  num2str(SNR_testing(SNR_idx))'' ' dB'],'interpreter','Latex'  )
        
        figure(7);
        title([  num2str(SNR_testing(SNR_idx))'' ' dB'],'interpreter','Latex'  )
      end
      
      if SNR_idx ~= 1
        figure(6)
        title([  num2str(SNR_testing(SNR_idx))'' ' dB'],'interpreter','Latex'  )
        
        figure(7)
        title([  num2str(SNR_testing(SNR_idx))'' ' dB'],'interpreter','Latex'  )
      end
    end
    
    if SNR_idx==1
      AVG_subopt = zeros(size( percentage_method_subopt_all{SNR_idx}));
      AVG_opt = zeros(size( percentage_method_opt_all{SNR_idx}));
    end
    
    %AVERAGE CASE
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
    figure(4),subplot(2,2,4)
    
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
    grid on
    
    if SNR_idx == length(SNR_testing)
      xlabel('Number of sources, q','interpreter','none')
      
      %         set(gca,'TickLabelInterpreter','Latex')
      %         h_legend = legend('AIC', 'HQ', 'MDL', 'AICc', 'KICvc', 'WIC','NT','Location','best');
      %         set(h_legend,'Interpreter','Latex')
      %
    end
    title('ALL','interpreter','Latex'  )
    
    % Plot percentage correct for optimal methods (average case)
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
    % set(gca,'TickLabelInterpreter','Latex')
    grid on
    if SNR_idx == length(SNR_testing)
      xlabel('Number of sources, q','interpreter','none')
      %          set(gca,'TickLabelInterpreter','Latex')
      %         h_legend = legend('AIC', 'HQ', 'MDL', 'AICc', 'KICvc', 'WIC','NT','Location','best');
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
  
  return;
end
