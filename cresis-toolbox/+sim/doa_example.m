% script sim.doa_example
%
% Example setup scripts for the sim.doa function. To run, enable one
% of the if {0|1} blocks below and then run the script.
%
% Author: John Paden, Theresa Stumpf
%
% See also: doa.m

physical_constants;

if 0
  % =======================================================================
  % Wax and Ziskind 1988 Fig 2
  % =======================================================================
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  fc = 312.5e6;
  BW = 1e6;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],[1 1 1],[1:3],3,c/fc/4};
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
  % Two equal power sources, sweep SNR
  % One source fixed at 0 deg, another fixed at 20 deg
  % Snapshots fixed at 10
  param.monte.SNR   = repmat(linspace(10,25,16).' - 10*log10(3), [1 2]);
  num_tests = size(param.monte.SNR,1);
  param.monte.DOA   = repmat([0 20],[num_tests 1]);
  param.monte.Nsnap = repmat(10,[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;

  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'D:\tmp\TSP_DOA';
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
  savefig(1,[out_fn '_src1.fig']);
  savefig(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
if 0
  % =======================================================================
  % Wax and Ziskind 1988 Fig 3
  % =======================================================================
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  fc = 312.5e6;
  BW = 1e6;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],[1 1 1],[1:3],3,c/fc/4};
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
  out_fn_dir = 'D:\tmp\TSP_DOA';
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
  savefig(1,[out_fn '_src1.fig']);
  savefig(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
if 0
  % =======================================================================
  % Wax and Ziskind 1988 Fig 4
  % =======================================================================
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  fc = 312.5e6;
  BW = 1e6;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],[1 1 1],[1:3],3,c/fc/4};
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
  out_fn_dir = 'D:\tmp\TSP_DOA';
  out_fn_name = 'wax_ziskind_fig4';
  
  RMSE = sim.doa_rmse(param,results);
  
  figure(1); clf
  plot(param.monte.Nsnap,RMSE(:,:,1).','.','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE'})
  title('Figure 3 from Wax and Ziskind 1988')
  
  figure(2); clf
  plot(param.monte.Nsnap,RMSE(:,:,2).','.','LineWidth',2);
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
  savefig(1,[out_fn '_src1.fig']);
  savefig(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
if 1
  % =======================================================================
  % Wax and Ziskind 1988 Fig 6
  % =======================================================================
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  fc = 312.5e6;
  BW = 1e6;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],ones(1,7),[1:7],7,c/fc/4};
  % DOA method parameters
  param.method.list                   = [2 7 8 9];
  param.method.Nsv                    = 21;
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
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;

  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'D:\tmp\TSP_DOA';
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
  savefig(1,[out_fn '_src1.fig']);
  savefig(2,[out_fn '_src2.fig']);
  savefig(3,[out_fn '_src3.fig']);
  savefig(4,[out_fn '_src4.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
if 0
  % =======================================================================
  % TSP_DOA Figure 1
  % =======================================================================
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],ones(1,8),[1:8],8,0.24};
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
  param.method.wb_fd.filter_banks     = 5;
  
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
  out_fn_dir = 'D:\tmp\TSP_DOA';
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
  savefig(1,[out_fn '_src1.fig']);
  savefig(2,[out_fn '_src2.fig']);
  savefig(3,[out_fn '_src3.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
if 0
  % =======================================================================
  % TSP_DOA Figure 2
  % =======================================================================
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],ones(1,8),[1:8],8,0.24};
  % DOA method parameters
  param.method.list                   = [2 7 8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[30 40]/180*pi,[50 60]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 9;
  
  % DOA monte carlo setup
  % Two sources, SNR swept from 0 to 20 dB, second source is 10 dB
  % Fixed at 35 and 55 deg
  % Snapshots fixed at 9 fast-time and 11 slow-time samples (i.e. 99)
  param.monte.SNR   = [linspace(0,20,11).', 10*ones(11,1)];
  num_tests = size(param.monte.SNR,1);
  param.monte.DOA   = repmat([35 55],[num_tests 1]);
  param.monte.Nsnap = repmat(9*11,[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;

  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'D:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig2';
  
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
  savefig(1,[out_fn '_src1.fig']);
  savefig(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
if 0
  % =======================================================================
  % TSP_DOA Figure 3
  % =======================================================================
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],ones(1,8),[1:8],8,0.24};
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
  out_fn_dir = 'D:\tmp\TSP_DOA';
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
  savefig(1,[out_fn '_src1.fig']);
  savefig(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  
if 1
  % =======================================================================
  % TSP_DOA Figure 4
  % =======================================================================
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],ones(1,8),[1:8],8,0.24};
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
  out_fn_dir = 'D:\tmp\TSP_DOA';
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
  savefig(1,[out_fn '_src1.fig']);
  savefig(2,[out_fn '_src2.fig']);
  savefig(3,[out_fn '_src3.fig']);
  savefig(4,[out_fn '_src4.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
  