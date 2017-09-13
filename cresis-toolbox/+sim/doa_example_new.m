% script sim.doa_example
%
% Example setup scripts for the sim.doa function. To run, enable one
% of the if {0|1} blocks below and then run the script.
%
% Author: John Paden, Theresa Stumpf
% Feb 06, 2017
% See also: doa.m

%Note: SNR is the signal-to-noise ratio per antenna element.
%      ASNR = Nsens*SNR is the total signal-to-noise ratio of the antenna array

clearvars -EXCEPT  gRadar fn_dir fn_idx fn_name fns pidx;
close all; clc;

physical_constants;

% =========================================================================
% ===================== Part 1: Theresa's thesis results ==================
% =========================================================================

%% Wax and Ziskind 1988 Fig 2: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =3, Nsnaps =10, SNR = variable, DoAs = [0 20],
%             fc =312.5MHz, BW =1MHz, W =1 , Nbands =1 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 312e6;
  param.src.f1                      = 313e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:3],[0; c/param.src.fc/2; 0]};
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
  param.monte.ASNR = repmat(linspace(10,25,16).',[1 2]);
  param.monte.SNR   = param.monte.ASNR - 10*log10(size(param.src.lever_arm.args{3},2));
  num_tests = size(param.monte.SNR,1);
  param.monte.DOA   = repmat([0 20],[num_tests 1]);
  param.monte.Nsnap = repmat(10,[num_tests 1]);
  param.monte.runs  = 200;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'wax_ziskind_fig2';
  
  RMSE = sim.doa_rmse(param,results);
  
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.ASNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.ASNR(:,1),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Source SNR (dB)')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Figure 2 from Wax and Ziskind 1988')
  
  figure(2); clf
  plot(param.monte.ASNR(:,2),RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.ASNR(:,2),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Source SNR (dB)')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
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

%% Wax and Ziskind 1988 Fig 3: Nsnaps vs RMSE
%
% Parameters: Nsrc =2, Nsens =3, Nsnaps =variable, SNR =[20 20], DoAs =[0 20],
%             fc =312.5MHz, BW =1MHz, W =1 , Nbands =1 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 312e6;
  param.src.f1                      = 313e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:3],[0; c/param.src.fc/2; 0]};
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
  param.monte.ASNR = repmat([20 20],[num_tests 1]);
  param.monte.SNR   = param.monte.ASNR - 10*log10(size(param.src.lever_arm.args{3},2));
  param.monte.DOA   = repmat([0 20],[num_tests 1]);
  param.monte.runs  = 200;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'wax_ziskind_fig3';
  
  RMSE = sim.doa_rmse(param,results);
  
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  semilogx(param.monte.Nsnap,RMSE(:,:,1).','+','LineWidth',2);
  hold on
  semilogx(param.monte.Nsnap,(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Figure 3 from Wax and Ziskind 1988')
  
  figure(2); clf
  semilogx(param.monte.Nsnap,RMSE(:,:,2).','+','LineWidth',2);
  hold on
  semilogx(param.monte.Nsnap,(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
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

%% Wax and Ziskind 1988 Fig 4: Nsnaps vs RMSE
%
% Parameters: Nsrc =2, Nsens =3, Nsnaps =variable, SNR =[20 20], DoAs =[0 5],
%             fc =312.5MHz, BW =1MHz, W =1 , Nbands =1 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 312e6;
  param.src.f1                      = 313e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:3],[0; c/param.src.fc/2; 0]};
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
  param.monte.Nsnap   = round(linspace(50,1000,11).');
  num_tests = size(param.monte.Nsnap,1);
  param.monte.ASNR = repmat([20 20],[num_tests 1]);
  param.monte.SNR   = param.monte.ASNR - 10*log10(size(param.src.lever_arm.args{3},2));
  param.monte.DOA   = repmat([0 5],[num_tests 1]);
  param.monte.runs  = 500;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'wax_ziskind_fig3';
  
  RMSE = sim.doa_rmse(param,results);
  
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.Nsnap,RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.Nsnap,(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Figure 4 from Wax and Ziskind 1988')
  
  figure(2); clf
  plot(param.monte.Nsnap,RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.Nsnap,(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
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

%% Wax and Ziskind 1988 Fig 6: Nsnaps vs RMSE
%
% Parameters: Nsrc =4, Nsens =7, Nsnaps =variable, SNR =[20 20 20 20], DoAs =[-8 -1 5 15],
%             fc =312.5MHz, BW =1MHz, W =1 , Nbands =1 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 312e6;
  param.src.f1                      = 313e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:7],[0; c/param.src.fc/2; 0]};
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
  param.monte.Nsnap   = round(linspace(50,1000,10).');%round(logspace(log10(50),log10(500),10).');
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
  
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.Nsnap,RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.Nsnap,(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Figure 6 from Wax and Ziskind 1988')
  
  figure(2); clf
  plot(param.monte.Nsnap,RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.Nsnap,(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Complement to figure 6 from Wax and Ziskind (source 1)')
  
  figure(3); clf
  plot(param.monte.Nsnap,RMSE(:,:,3).','+','LineWidth',2);
  hold on
  plot(param.monte.Nsnap,(sqrt(CRB_angular(3,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Complement to figure 6 from Wax and Ziskind (source 2)')
  
  figure(4); clf
  plot(param.monte.Nsnap,RMSE(:,:,4).','+','LineWidth',2);
  hold on
  plot(param.monte.Nsnap,(sqrt(CRB_angular(4,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
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

%% Theresa, Figure 4-7: DOA vs RMSE
%
% Parameters: Nsrc =1, Nsens =8, Nsnaps =100, SNR =[20], DoAs =variable,
%             fc =325MHz, BW =250MHz, W =3 , Nbands =11 .
% =======================================================================
if 0
  %% Setup simulation parameters
  DOArange = [0:5:45]; % the left side is the mirror of the right side
  RMSE_Result = [];
  for doa_idx = 1:length(DOArange)
    param = [];
    
    % Source parameters
    param.src.f0                      = 200e6;
    param.src.f1                      = 450e6;
    param.src.fc                      = (param.src.f1+param.src.f0)/2;
    param.src.ft_wind                 = @boxcar;
    param.src.lever_arm.fh            = @sim.lever_arm_example;
    param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};%0.24
    % DOA method parameters
    param.method.list                   = [2 7 8 9];
    param.method.Nsv                    = 24;
    param.method.OneD_Nsv               = 128;
    param.method.src_limits             = {[DOArange(doa_idx)-15 DOArange(doa_idx)+15]/180*pi};%{[min(DOArange) max(DOArange)]/180*pi};
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
    param.monte.SNR   = repmat(linspace(0,20,1).', [1 1]);
    num_tests = size(param.monte.SNR,1);
    param.monte.DOA   = repmat([DOArange(doa_idx)],[num_tests 1]);
    param.monte.Nsnap = repmat(100,[num_tests 1]);
    param.monte.runs  = 100;
    param.monte.random_seed_offset = 0;
    
    %% Run the simulation
    results = sim.doa(param);
    
    %% Process and save the outputs
    
    RMSE = sim.doa_rmse(param,results);  % 4X1X1
    
    [CRB_angular, CRB_spatial] = CRB(param);
    
    RMSE_Result(doa_idx,:) = squeeze(RMSE);
    CRB_Result(doa_idx)    = squeeze(CRB_angular);
  end
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig4.7';
  
  figure(1); clf
  plot(DOArange,RMSE_Result,'+','LineWidth',2);
  hold on
  plot(DOArange,(sqrt(CRB_Result)*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('DOA vs RMSE')
  xlim([min(DOArange) max(DOArange)])
  
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% Theresa, Figure 4-7 (modified): DOA vs RMSE
%
% Parameters: Nsrc =1, Nsens =3, Nsnaps =10, SNR =[20], DoAs =variable,
%             fc =325MHz, BW =250MHz, W =3 , Nbands =11 .
% =======================================================================
if 0
  %% Setup simulation parameters
  DOArange = [0:5:45]; % the left side is the mirror of the right side
  RMSE_Result = [];
  for doa_idx = 1:length(DOArange)
    param = [];
    
    % Source parameters
    param.src.f0                      = 200e6;
    param.src.f1                      = 450e6;
    param.src.fc                      = (param.src.f1+param.src.f0)/2;
    param.src.ft_wind                 = @boxcar;
    param.src.lever_arm.fh            = @sim.lever_arm_example;
    param.src.lever_arm.args          = {[],1,[1:3],[0; c/param.src.fc/2; 0]};%0.24
    % DOA method parameters
    param.method.list                   = [2 7 8 9];
    param.method.Nsv                    = 24;
    param.method.OneD_Nsv               = 128;
    param.method.src_limits             = {[DOArange(doa_idx)-15 DOArange(doa_idx)+15]/180*pi};%{[min(DOArange) max(DOArange)]/180*pi};
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
    param.monte.SNR   = repmat(linspace(0,20,1).', [1 1]);
    num_tests = size(param.monte.SNR,1);
    param.monte.DOA   = repmat([DOArange(doa_idx)],[num_tests 1]);
    param.monte.Nsnap = repmat(100,[num_tests 1]);
    param.monte.runs  = 100;
    param.monte.random_seed_offset = 0;
    
    %% Run the simulation
    results = sim.doa(param);
    
    %% Process and save the outputs
    
    RMSE = sim.doa_rmse(param,results);  % 4X1X1
    
    [CRB_angular, CRB_spatial] = CRB(param);
    
    RMSE_Result(doa_idx,:) = squeeze(RMSE);
    CRB_Result(doa_idx)    = squeeze(CRB_angular);
  end
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig4.7';
  
  figure(1); clf
  plot(DOArange,RMSE_Result,'+','LineWidth',2);
  hold on
  plot(DOArange,(sqrt(CRB_Result)*180/pi),'+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('DOA vs RMSE')
  xlim([min(DOArange) max(DOArange)])
  
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% Theresa, Figure 4-8: Nsnaps vs RMSE
%
% Parameters: Nsrc =1, Nsens =8, Nsnaps =variable, SNR =[-5], DoAs =[25],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =11 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [2 7 8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[10 40]/180*pi};
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
  param.monte.Nsnap = round(logspace(log10(50),log10(1000),11).');
  num_tests = size(param.monte.Nsnap,1);
  param.monte.SNR   = repmat(linspace(-5,-5,1).', [num_tests 1]);
  param.monte.DOA   = repmat([25],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  semilogx(param.monte.Nsnap,RMSE(:,:,1).','+','LineWidth',2);
  hold on
  semilogx(param.monte.Nsnap,(sqrt(CRB_angular)*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('RMSE ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% Theresa, Figure 4-9: Nsnaps vs RMSE
%
% Parameters: Nsrc =1, Nsens =8, Nsnaps =variable, SNR =[20], DoAs =[25],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =11 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [2 7 8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[10 40]/180*pi};
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
  param.monte.Nsnap = round(logspace(log10(50),log10(1000),11).');
  num_tests = size(param.monte.Nsnap,1);
  param.monte.SNR   = repmat(linspace(20,20,1).', [num_tests 1]);
  param.monte.DOA   = repmat([25],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  semilogx(param.monte.Nsnap,RMSE(:,:,1).','+','LineWidth',2);
  hold on
  semilogx(param.monte.Nsnap,(sqrt(CRB_angular)*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('RMSE ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% Theresa, Figure 4-11: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =40, SNR =variable], DoAs =[25 60],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =3 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[20 30]/180*pi, [55 65]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 3;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(-5,35,11).',linspace(10,10,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(40,[num_tests 1]);
  param.monte.DOA   = repmat([25 60],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,1),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% Theresa, Figure 4-12: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =40, SNR =variable], DoAs =[25 60],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =3 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[20 30]/180*pi, [55 65]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 3;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(10,10,11).',linspace(-5,35,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(40,[num_tests 1]);
  param.monte.DOA   = repmat([25 60],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,2),RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,2),(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src2 ')
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% Theresa, Figure 4-13: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =40, SNR =variable], DoAs =[25 60],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =31 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[20 30]/180*pi, [55 65]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 31;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(-5,35,11).',linspace(10,10,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(40,[num_tests 1]);
  param.monte.DOA   = repmat([25 60],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,1),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% Theresa, Figure 4-14: Nsnaps vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =40, SNR =variable], DoAs =[25 60],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =31 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[20 30]/180*pi, [55 65]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 31;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(10,10,11).',linspace(-5,35,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(40,[num_tests 1]);
  param.monte.DOA   = repmat([25 60],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,2),RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,2),(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src2 ')
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  
  return;
end


%% Theresa, Figure 4-15/16: Nsnaps vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =variable, SNR =[5 10], DoAs =[25 60],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =31 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[10 40]/180*pi, [35 80]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 31;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.Nsnap   = round(logspace(log10(40),log10(1000),11).');
  num_tests = size(param.monte.Nsnap,1);
  param.monte.SNR   = repmat([5 10],[num_tests 1]);
  param.monte.DOA   = repmat([25 60],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.Nsnap,RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.Nsnap,(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 ')
  
  figure(2); clf
  plot(param.monte.Nsnap,RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.Nsnap,(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src2 ')
  
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

%% Theresa, Figure 4-17/18: Nsnaps vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =variable, SNR =[25 10], DoAs =[25 60],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =31 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[10 40]/180*pi, [35 80]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 31;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.Nsnap   = round(logspace(log10(40),log10(1000),11).');
  num_tests = size(param.monte.Nsnap,1);
  param.monte.SNR   = repmat([20 10],[num_tests 1]);
  param.monte.DOA   = repmat([25 60],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.Nsnap,RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.Nsnap,(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 ')
  
  figure(2); clf
  plot(param.monte.Nsnap,RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.Nsnap,(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Snapshots')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src2 ')
  
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

%% Theresa, Figure 4-20: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =40, SNR =variable], DoAs =[-8 +8],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =31 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[-20 5]/180*pi, [-5 20]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 31;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(-5,35,11).',linspace(-5,35,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(40,[num_tests 1]);
  param.monte.DOA   = repmat([-8 +8],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,1),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

% =========================================================================
% ==================== Part 2: PAST2016 presentation results ==============
% =========================================================================

%% PAST2016 Presentation, Slide 19: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =100, SNR =variable], DoAs =[25 60],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =4 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[20 30]/180*pi, [55 65]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 4;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(-5,35,11).',linspace(-5,35,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(100,[num_tests 1]);
  param.monte.DOA   = repmat([25 60],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,1),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 ')
  
  figure(2); clf
  plot(param.monte.SNR(:,2),RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,2),(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src2 ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% PAST2016 Presentation, Slide 20: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =100, SNR =variable], DoAs =[25 60],
%             fc =325MHz, BW =250MHz, W =1 , Nbands =4 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[20 30]/180*pi, [55 65]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 1;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 4;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(-5,35,11).',linspace(-5,35,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(100,[num_tests 1]);
  param.monte.DOA   = repmat([25 60],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,1),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 ')
  
  figure(2); clf
  plot(param.monte.SNR(:,2),RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,2),(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src2 ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% PAST2016 Presentation, Slide 21: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =100, SNR =variable], DoAs =[25 60],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =20 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[-90 90]/180*pi, [-90 90]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 20;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(-5,35,11).',linspace(-5,35,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(100,[num_tests 1]);
  param.monte.DOA   = repmat([25 60],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,1),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 ')
  
  figure(2); clf
  plot(param.monte.SNR(:,2),RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,2),(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src2 ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% PAST2016 Presentation, Slide 22-src1: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =100, SNR =variable], DoAs =[25 60],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =20 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[20 30]/180*pi, [55 65]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 20;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(-5,35,11).',linspace(10,10,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(100,[num_tests 1]);
  param.monte.DOA   = repmat([25 60],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,1),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% PAST2016 Presentation, Slide 22-src2: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =100, SNR =variable], DoAs =[25 60],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =20 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[20 30]/180*pi, [55 65]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 20;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(10,10,11).',linspace(-5,35,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(100,[num_tests 1]);
  param.monte.DOA   = repmat([25 60],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,2),RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,2),(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src2 ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% PAST2016 Presentation, Slide 23: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =20, SNR =variable], DoAs =[-5 5],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =4 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[-20 10]/180*pi, [-10 20]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 4;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(10,35,11).',linspace(10,35,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(20,[num_tests 1]);
  param.monte.DOA   = repmat([-5 5],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,1),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 ')
  
  figure(2); clf
  plot(param.monte.SNR(:,2),RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,2),(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src2 ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% PAST2016 Presentation, Slide 24: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =20, SNR =variable], DoAs =[-10 10],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =4 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[-20 10]/180*pi, [-10 20]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 4;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(10,35,11).',linspace(10,35,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(20,[num_tests 1]);
  param.monte.DOA   = repmat([-10 10],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,1),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 ')
  
  figure(2); clf
  plot(param.monte.SNR(:,2),RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,2),(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src2 ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% PAST2016 Presentation, Slide 25: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =20, SNR =variable], DoAs =[-10 10],
%             fc =325MHz, BW =250MHz, W =3 , Nbands =4 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[-25 -15]/180*pi, [15 25]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 4;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(10,35,11).',linspace(10,35,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(20,[num_tests 1]);
  param.monte.DOA   = repmat([-20 20],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,1),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 ')
  
  figure(2); clf
  plot(param.monte.SNR(:,2),RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,2),(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src2 ')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% PAST2016 Presentation, Slide 26: SNR vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =100, SNR =variable], DoAs =[-60 -10 10 60]==>[-60 60] are clutter,
%             fc =325MHz, BW =250MHz, W =3 , Nbands =4 .
% =======================================================================
if 0
  %% Setup simulation parameters
  param = [];
  
  % Source parameters
  param.src.f0                      = 200e6;
  param.src.f1                      = 450e6;
  param.src.fc                      = (param.src.f1+param.src.f0)/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [8 9];
  param.method.Nsv                    = 24;
  param.method.OneD_Nsv               = 128;
  param.method.src_limits             = {[-75 -45]/180*pi, [-40 20]/180*pi,[-20 40]/180*pi, [45 75]/180*pi};
  param.method.theta_guard            = 1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 3;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 4;
  
  % DOA monte carlo setup
  % Three equal power sources, SNR swept from 0 to 20 dB
  % Fixed at -5, 0, 5 deg
  % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
  param.monte.SNR   = repmat([linspace(10,35,11).',linspace(10,35,11).',linspace(10,35,11).',linspace(10,35,11).'],[1 1]);
  num_tests = size(param.monte.SNR,1);
  param.monte.Nsnap = repmat(100,[num_tests 1]);
  param.monte.DOA   = repmat([-60 -10 10 60],[num_tests 1]);
  param.monte.runs  = 100;
  param.monte.random_seed_offset = 0;
  
  %% Run the simulation
  results = sim.doa(param);
  
  %% Process and save the outputs
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig1';
  
  RMSE = sim.doa_rmse(param,results);
  [CRB_angular, CRB_spatial] = CRB(param);
  
  figure(1); clf
  plot(param.monte.SNR(:,1),RMSE(:,:,1).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,1),(sqrt(CRB_angular(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src1 (clutter) ')
  
  figure(2); clf
  plot(param.monte.SNR(:,2),RMSE(:,:,2).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,2),(sqrt(CRB_angular(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src2 ')
  
  figure(3); clf
  plot(param.monte.SNR(:,3),RMSE(:,:,3).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,3),(sqrt(CRB_angular(3,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src3 ')
  
  figure(4); clf
  plot(param.monte.SNR(:,4),RMSE(:,:,4).','+','LineWidth',2);
  hold on
  plot(param.monte.SNR(:,4),(sqrt(CRB_angular(4,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('SNR [dB]')
  ylabel('RMS error ( \circ )')
  legend({'WB','WBMLE','CRLB'})
  title('RMSE-src4 (clutter)')
  
  
  % Save Outputs
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

% =========================================================================
% ==================== Part 3: Extra (non-standard) results ===============
% =========================================================================

%% Result-1: Number of subbands vs RMSE
%
% Parameters: Nsrc =1, Nsens =3 and 8, Nsnaps =100, SNR =[10], DoAs =single src in [0 20 60],
%             fc =325MHz, BW =250MHz, W =ceil(tau*BW), Nbands =variable .
% =======================================================================
if 0
  %% Setup simulation parameters
  num_filter_banks = 1:2:30;
  DOA = [0 20 60];
  RMSE_Result = [];
  for doa_idx = 1:length(DOA)
    for idx = 1:length(num_filter_banks)
      param = [];
      
      % Source parameters
      param.src.f0                      = 200e6;
      param.src.f1                      = 450e6;
      param.src.fc                      = (param.src.f1+param.src.f0)/2;
      param.src.ft_wind                 = @boxcar;
      param.src.lever_arm.fh            = @sim.lever_arm_example;
      param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};%0.24
      % DOA method parameters
      param.method.list                   = [2 7 8 9];
      param.method.Nsv                    = 24;
      param.method.OneD_Nsv               = 128;
      param.method.src_limits             = {[DOA(doa_idx)-10 DOA(doa_idx)+10]/180*pi};
      param.method.theta_guard            = 1.5/180*pi;
      param.method.nb_nd.init             = 'ap';
      param.method.wb_td.init             = 'ap';
      param.method.wb_fd.init             = 'ap';
      %       param.method.wb_td.widening_factor = 1;
      param.method.wb_fd.filter_banks     = num_filter_banks(idx);
      
      % DOA monte carlo setup
      % Three equal power sources, SNR swept from 0 to 20 dB
      % Fixed at -5, 0, 5 deg
      % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
      param.monte.SNR   = repmat(linspace(10,10,1).', [1 1]);
      num_tests = size(param.monte.SNR,1);
      param.monte.DOA   = repmat([DOA(doa_idx)],[num_tests 1]);
      param.monte.Nsnap = repmat(100,[num_tests 1]);
      param.monte.runs  = 500;
      param.monte.random_seed_offset = 0;
      
      N     = size(param.src.lever_arm.args{3},2);
      fc    = param.src.fc;
      B     = param.src.f1-param.src.f0;
      theta = param.monte.DOA;
      param.method.wb_td.widening_factor = ceil((N-1)/2*sin(theta*pi/180)*(B/fc));
      
      if mod(param.method.wb_td.widening_factor,2) == 0
        param.method.wb_td.widening_factor = param.method.wb_td.widening_factor + 1;
      end
      
      %% Run the simulation
      results = sim.doa(param);
      
      %% Process and save the outputs
      
      RMSE = sim.doa_rmse(param,results);  % 4X1X1
      
      [CRB_angular, CRB_spatial] = CRB(param);
      
      RMSE_Result(doa_idx,idx,:) = squeeze(RMSE);
      CRB_Result(doa_idx,idx)    = squeeze(CRB_angular);
    end
  end
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig4.7';
  
  
  figure(1); clf
  subplot(311)
  plot(num_filter_banks,squeeze(RMSE_Result(1,:,:)),'+','LineWidth',2);
  hold on
  plot(num_filter_banks,(sqrt(CRB_Result(1,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Number of sub-bands')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('DoA = 0^o')
  
  figure(1);
  subplot(312)
  plot(num_filter_banks,squeeze(RMSE_Result(2,:,:)),'+','LineWidth',2);
  hold on
  plot(num_filter_banks,(sqrt(CRB_Result(2,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Number of sub-bands')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('DoA = 20^o')
  
  figure(1);
  subplot(313)
  plot(num_filter_banks,squeeze(RMSE_Result(3,:,:)),'+','LineWidth',2);
  hold on
  plot(num_filter_banks,(sqrt(CRB_Result(3,:))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Number of sub-bands')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('DoA = 60^o')
  
  %   title('Nb vs RMSE')
  
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
% Comment on Result 1: When N_sensors =3, MUSIC, MLE, and WBDOA, showed same performance
%because eventhough the BW is large, but the array length is small. Thus,
%regardless of the DoA, the system is narrowband because ceil(tau*BW)=1 for
%all DoAs. Also, the DCM of these three methods is not a function of Nb,
%thus it will not change with Nb.
%However, the WBMLE deviates from the other three methods when DoA
%increased (w.r.t the other three methods) it doesnot depend on the idea of
%the widenning factor; meaning the DCM depends on the number of subbands,
%not on the widenning factor.  Also, after some value of Nb, the RMSE of
%WBMLE stays at the same level and doesn't improve: has the same distance
%to CRB and to the other three methods. Thus, it is not always good to
%increase Nb.

%% Result-2: Number of sources vs RMSE [problem in the code]
%
% Parameters: Nsrc =variable, Nsens =15, Nsnaps =100, SNR =[10], DoAs = variable,
%             fc =325MHz, BW =250MHz, W =ceil(tau*BW), Nbands =W .
% =======================================================================
if 0
  %% Setup simulation parameters
  DOA = [0 -10 10 -20 20];% -45 45 -60 60];%[-30:10:30]
  limits1 = (DOA-10)*pi/180;
  limits2 = (DOA+10)*pi/180;
  RMSE_Result = [];
  
  for doa_idx = 1:length(DOA)
    param = [];
    
    % Source parameters
    param.src.f0                      = 200e6;
    param.src.f1                      = 450e6;
    param.src.fc                      = (param.src.f1+param.src.f0)/2;
    param.src.ft_wind                 = @boxcar;
    param.src.lever_arm.fh            = @sim.lever_arm_example;
    param.src.lever_arm.args          = {[],1,[1:10],[0; c/param.src.fc/2; 0]};%0.24
    % DOA method parameters
    param.method.list                   = [2 7 8 9];
    param.method.Nsv                    = 24;
    param.method.OneD_Nsv               = 128;
    %       param.method.src_limits             = {[DOA-10 DOA+10]/180*pi};
    param.method.theta_guard            = 1.5/180*pi;
    param.method.nb_nd.init             = 'ap';
    param.method.wb_td.init             = 'ap';
    param.method.wb_fd.init             = 'ap';
    param.method.wb_td.widening_factor = 3;
    param.method.wb_fd.filter_banks     = 3;
    
    % DOA monte carlo setup
    % Three equal power sources, SNR swept from 0 to 20 dB
    % Fixed at -5, 0, 5 deg
    % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
    param.monte.DOA   = repmat([DOA(1:doa_idx)],[1 1]);
    num_tests         = size(param.monte.DOA,2);
    param.monte.SNR   = repmat(linspace(10,10,1).', [num_tests 1])';
    param.monte.Nsnap = repmat(100,[num_tests 1]);
    param.monte.runs  = 20;
    param.monte.random_seed_offset = 0;
    
    %     N     = size(param.src.lever_arm.args{3},2);
    %     fc    = param.src.fc;
    %     B     = param.src.f1-param.src.f0;
    %     theta = param.monte.DOA(end);
    %     param.method.wb_td.widening_factor = ceil(abs((N-1)/2*sin(theta*pi/180)*(B/fc)));
    %
    %     if mod(param.method.wb_td.widening_factor,2) == 0
    %       param.method.wb_td.widening_factor = param.method.wb_td.widening_factor + 1;
    %     end
    %
    %     param.method.wb_fd.filter_banks      = param.method.wb_td.widening_factor;
    
    for idx = 1:numel(param.monte.DOA)
      param.method.src_limits{idx}       = [limits1(idx) limits2(idx)];
    end
    %% Run the simulation
    results = sim.doa(param);
    
    %% Process and save the outputs
    
    RMSE = sim.doa_rmse(param,results);  % 4X1X1
    
    [CRB_angular, CRB_spatial] = CRB(param,1);
    
    RMSE_Result(doa_idx,:) = squeeze(RMSE(:,:,1));
    CRB_Result(doa_idx)    = squeeze(CRB_angular);
  end
  
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig4.7';
  
  
  figure(1); clf
  plot(1:length(DOA),RMSE_Result,'+','LineWidth',2);
  hold on
  plot(1:length(DOA),(sqrt(CRB_Result)*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Number of sourcess')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('src1 [DoA = 0^o]')
  
  
  %   title('Nb vs RMSE')
  
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end


%% Result-3: BW vs RMSE
%
% Parameters: Nsrc =1, Nsens =8 and 3, Nsnaps =100, SNR =[10], DoAs = [10]
%             fc =1000MHz, BW =variable, W =ceil(tau*BW), Nbands =3 and 21 .
% =======================================================================
if 0
  %% Setup simulation parameters
  RMSE_Result = [];
  fc = 1000e6;
  BW = round(logspace(log10(1e6),log10(2*fc),11).');
  for bw_idx = 1:length(BW)
    param = [];
    
    % Source parameters
    
    param.src.fc                      = fc;
    param.src.f0                      = param.src.fc - BW(bw_idx)/2;
    param.src.f1                      = param.src.fc + BW(bw_idx)/2;
    param.src.ft_wind                 = @boxcar;
    param.src.lever_arm.fh            = @sim.lever_arm_example;
    param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};%0.24
    % DOA method parameters
    param.method.list                   = [2 7 8 9];
    param.method.Nsv                    = 24;
    param.method.OneD_Nsv               = 128;
    param.method.src_limits             = {[0 20]/180*pi};
    param.method.theta_guard            = 1.5/180*pi;
    param.method.nb_nd.init             = 'ap';
    param.method.wb_td.init             = 'ap';
    param.method.wb_fd.init             = 'ap';
    %       param.method.wb_td.widening_factor = 1;
    param.method.wb_fd.filter_banks     = 3;
    
    % DOA monte carlo setup
    % Three equal power sources, SNR swept from 0 to 20 dB
    % Fixed at -5, 0, 5 deg
    % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
    param.monte.SNR   = repmat(linspace(10,10,1).', [1 1]);
    num_tests = size(param.monte.SNR,1);
    param.monte.DOA   = repmat([10],[num_tests 1]);
    param.monte.Nsnap = repmat(100,[num_tests 1]);
    param.monte.runs  = 100;
    param.monte.random_seed_offset = 0;
    
    N     = size(param.src.lever_arm.args{3},2);
    fc    = param.src.fc;
    B     = param.src.f1-param.src.f0;
    theta = param.monte.DOA;
    param.method.wb_td.widening_factor = ceil((N-1)/2*sin(theta*pi/180)*(B/fc));
    
    if mod(param.method.wb_td.widening_factor,2) == 0
      param.method.wb_td.widening_factor = param.method.wb_td.widening_factor + 1;
    end
    
    %% Run the simulation
    results = sim.doa(param);
    
    %% Process and save the outputs
    
    RMSE = sim.doa_rmse(param,results);  % 4X1X1
    
    [CRB_angular, CRB_spatial] = CRB(param);
    
    RMSE_Result(bw_idx,:) = squeeze(RMSE);
    CRB_Result(bw_idx)    = squeeze(CRB_angular);
  end
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig4.7';
  
  
  figure(1); clf
  semilogx(BW,squeeze(RMSE_Result),'+','LineWidth',2);
  hold on
  semilogx(BW,(sqrt(CRB_Result)*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Signal BW [Hz]')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('BW vs RMSE')
  
  
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end


%% Result-4: Fractional BW (i.e. BW/fc) vs RMSE
% Note: BW/fc better explains the effect of the BW than BW only.
% Parameters: Nsrc =1, Nsens =8, Nsnaps =100, SNR =[10], DoAs = [60]
%             fc =1000MHz, BW =variable, W =ceil(tau*BW), Nbands =3 and 21 .
% =======================================================================
if 0
  %% Setup simulation parameters
  RMSE_Result = [];
  fc = 1000e6;
  BW = fc*linspace(0.01,2,11)';
  for bw_idx = 1:length(BW)
    param = [];
    
    % Source parameters
    param.src.fc                      = fc;
    param.src.f0                      = param.src.fc - BW(bw_idx)/2;
    param.src.f1                      = param.src.fc + BW(bw_idx)/2;
    param.src.ft_wind                 = @boxcar;
    param.src.lever_arm.fh            = @sim.lever_arm_example;
    param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};%0.24
    % DOA method parameters
    param.method.list                   = [2 7 8 9];
    param.method.Nsv                    = 24;
    param.method.OneD_Nsv               = 128;
    param.method.src_limits             = {[50 70]/180*pi};
    param.method.theta_guard            = 1.5/180*pi;
    param.method.nb_nd.init             = 'ap';
    param.method.wb_td.init             = 'ap';
    param.method.wb_fd.init             = 'ap';
    %       param.method.wb_td.widening_factor = 1;
    param.method.wb_fd.filter_banks     = 21;
    
    % DOA monte carlo setup
    % Three equal power sources, SNR swept from 0 to 20 dB
    % Fixed at -5, 0, 5 deg
    % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
    param.monte.SNR   = repmat(linspace(10,10,1).', [1 1]);
    num_tests = size(param.monte.SNR,1);
    param.monte.DOA   = repmat([60],[num_tests 1]);
    param.monte.Nsnap = repmat(100,[num_tests 1]);
    param.monte.runs  = 500;
    param.monte.random_seed_offset = 0;
    
    N     = size(param.src.lever_arm.args{3},2);
    fc    = param.src.fc;
    B     = param.src.f1-param.src.f0;
    theta = param.monte.DOA;
    param.method.wb_td.widening_factor = ceil((N-1)/2*sin(theta*pi/180)*(B/fc));
    
    if mod(param.method.wb_td.widening_factor,2) == 0
      param.method.wb_td.widening_factor = param.method.wb_td.widening_factor + 1;
    end
    
    %% Run the simulation
    results = sim.doa(param);
    
    %% Process and save the outputs
    
    RMSE = sim.doa_rmse(param,results);  % 4X1X1
    
    [CRB_angular, CRB_spatial] = CRB(param);
    
    RMSE_Result(bw_idx,:) = squeeze(RMSE);
    CRB_Result(bw_idx)    = squeeze(CRB_angular);
  end
  %   end
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig4.7';
  
  
  figure(1); clf
  %   subplot(311)
  plot(BW/fc,squeeze(RMSE_Result),'+','LineWidth',2);
  hold on
  plot(BW/fc,(sqrt(CRB_Result)*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Fractional BW [unitless]')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Fractional BW vs RMSE')
  
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
% Note: As W increases, WBDOA becomes more resistent to the increase in the
% BW (in WBMLE, the effect of Nb is equivalent to the effect of W in WBDOA).
% Thus, as DOA increases (goes away from nadir direction), which leads to
% increasing W, the WBDOA approaches CRB. Thus, numerically we conclude:
% 1) WBDOA achieves CRB for large values of W even for small number of
% antennas, small Nsnaps, or low SNR.
% 2) WBMLE approaches CRB in the large Nb regime.

%% Result-5: W vs RMSE [one source]
%
% Parameters: Nsrc =1, Nsens =15, Nsnaps =100, SNR =[10], DoAs = [60]
%             fc =100MHz, BW =11.5*fc, W =variable, Nbands =3 .
% =======================================================================
if 0
  %% Setup simulation parameters
  RMSE_Result = [];
  fc = 100e6;
  BW = 1.5*fc;
  % W_true is the calculated value of W when Nsens=15, DOA=60, BW/fc=1.5;
  W_true = 7;
  W = [1:2:W_true+5];
  for w_idx = 1:length(W)
    param = [];
    
    % Source parameters
    param.src.fc                      = fc;
    param.src.f0                      = param.src.fc - BW/2;
    param.src.f1                      = param.src.fc + BW/2;
    param.src.ft_wind                 = @boxcar;
    param.src.lever_arm.fh            = @sim.lever_arm_example;
    param.src.lever_arm.args          = {[],1,[1:15],[0; c/param.src.fc/2; 0]};%0.24
    % DOA method parameters
    param.method.list                   = [2 7 8 9];
    param.method.Nsv                    = 24;
    param.method.OneD_Nsv               = 128;
    param.method.src_limits             = {[50 70]/180*pi};
    param.method.theta_guard            = 1.5/180*pi;
    param.method.nb_nd.init             = 'ap';
    param.method.wb_td.init             = 'ap';
    param.method.wb_fd.init             = 'ap';
    param.method.wb_td.widening_factor = W(w_idx);
    param.method.wb_fd.filter_banks     = 3;
    
    % DOA monte carlo setup
    % Three equal power sources, SNR swept from 0 to 20 dB
    % Fixed at -5, 0, 5 deg
    % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
    param.monte.SNR   = repmat(linspace(10,10,1).', [1 1]);
    num_tests = size(param.monte.SNR,1);
    param.monte.DOA   = repmat([60],[num_tests 1]);
    param.monte.Nsnap = repmat(100,[num_tests 1]);
    param.monte.runs  = 500;
    param.monte.random_seed_offset = 0;
    
    %% Run the simulation
    results = sim.doa(param);
    
    %% Process and save the outputs
    
    RMSE = sim.doa_rmse(param,results);  % 4X1X1
    
    [CRB_angular, CRB_spatial] = CRB(param);
    
    RMSE_Result(w_idx,:) = squeeze(RMSE);
    CRB_Result(w_idx)    = squeeze(CRB_angular);
  end
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig4.7';
  
  
  figure(1); clf
  plot(W,squeeze(RMSE_Result),'+','LineWidth',2);
  hold on
  plot(W,(sqrt(CRB_Result)*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Fractional BW [unitless]')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Widenning factor vs RMSE')
  
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
% Note: When W=W_true, WBDOA is very close to CRB. Thus, if we choose W to
% be the next (or the prevous) odd number greater (smaller) than W_true
% (i.e. W_true-2<=W=<W_true+2), then it will be almost identical to CRB.
% [NEEDS TO BE PROVED MATHEMATICALLY].

%% Result-6: W vs RMSE [two sources]
%
% Parameters: Nsrc =2, Nsens =15, Nsnaps =100, SNR =[10], DoAs = [50 60]
%             fc =100MHz, BW =1.5*fc, W =variable, Nbands =3 .
% =======================================================================
if 0
  %% Setup simulation parameters
  RMSE_Result = [];
  fc = 100e6;
  BW = 1.5*fc;
  % W_true is the calculated value of W when Nsens=15, DOA=60, BW/fc=1.5;
  W_true = 7;
  W = [1:2:W_true+5];
  for w_idx = 1:length(W)
    param = [];
    
    % Source parameters
    param.src.fc                      = fc;
    param.src.f0                      = param.src.fc - BW/2;
    param.src.f1                      = param.src.fc + BW/2;
    param.src.ft_wind                 = @boxcar;
    param.src.lever_arm.fh            = @sim.lever_arm_example;
    param.src.lever_arm.args          = {[],1,[1:15],[0; c/param.src.fc/2; 0]};%0.24
    % DOA method parameters
    param.method.list                   = [2 7 8 9];
    param.method.Nsv                    = 24;
    param.method.OneD_Nsv               = 128;
    param.method.src_limits             = {[20 40]/180*pi, [50 70]/180*pi};
    param.method.theta_guard            = 1.5/180*pi;
    param.method.nb_nd.init             = 'ap';
    param.method.wb_td.init             = 'ap';
    param.method.wb_fd.init             = 'ap';
    param.method.wb_td.widening_factor = W(w_idx);
    param.method.wb_fd.filter_banks     = 3;
    
    % DOA monte carlo setup
    % Three equal power sources, SNR swept from 0 to 20 dB
    % Fixed at -5, 0, 5 deg
    % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
    param.monte.SNR   = repmat(linspace(10,10,1), [1 2]);
    num_tests = size(param.monte.SNR,1);
    param.monte.DOA   = repmat([30 60],[num_tests 1]);
    param.monte.Nsnap = repmat(100,[num_tests 1]);
    param.monte.runs  = 500;
    param.monte.random_seed_offset = 0;
    
    %% Run the simulation
    results = sim.doa(param);
    
    %% Process and save the outputs
    
    RMSE = sim.doa_rmse(param,results);  % 4X1X1
    
    [CRB_angular, CRB_spatial] = CRB(param);
    
    RMSE_Result(w_idx,:,:) = squeeze(RMSE);
    CRB_Result(w_idx,:)    = squeeze(CRB_angular);
  end
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig4.7';
  
  
  figure(1); clf
  plot(W,squeeze(RMSE_Result(:,:,1)),'+','LineWidth',2);
  hold on
  plot(W,(sqrt(CRB_Result(:,1))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Fractional BW [unitless]')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Widenning factor vs RMSE')
  
  figure(2); clf
  plot(W,squeeze(RMSE_Result(:,:,2)),'+','LineWidth',2);
  hold on
  plot(W,(sqrt(CRB_Result(:,2))*180/pi).','+','LineWidth',2);
  grid on
  xlabel('Fractional BW [unitless]')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Widenning factor vs RMSE')
  
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
% Note: Same conclusion as Result-5
%% Result-7: Sector size vs RMSE (one source)
%
% Parameters: Nsrc =1, Nsens =8, Nsnaps =100, SNR =[10], DoAs = [10]
%             fc =500MHz, BW =1*fc, W =3, Nbands =3 .
% =======================================================================
if 0
  %% Setup simulation parameters
  RMSE_Result = [];
  fc = 500e6;
  BW = 1.5*fc;
  DOA = 10;
  for idx = 1:8
    param = [];
    
    % Source parameters
    param.src.fc                      = fc;
    param.src.f0                      = param.src.fc - BW/2;
    param.src.f1                      = param.src.fc + BW/2;
    param.src.ft_wind                 = @boxcar;
    param.src.lever_arm.fh            = @sim.lever_arm_example;
    param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};%0.24
    % DOA method parameters
    param.method.list                   = [2 7 8 9];
    param.method.Nsv                    = 24;
    param.method.OneD_Nsv               = 128;
    %     param.method.src_limits             = {[DOA-idx*10 DOA+idx*10]/180*pi};
    param.method.src_limits             = {[-90+10*(idx-1) 90-10*(idx-1)]/180*pi};
    param.method.theta_guard            = 1.5/180*pi;
    param.method.nb_nd.init             = 'ap';
    param.method.wb_td.init             = 'ap';
    param.method.wb_fd.init             = 'ap';
    param.method.wb_td.widening_factor = 3;
    param.method.wb_fd.filter_banks     = 3;
    
    % DOA monte carlo setup
    % Three equal power sources, SNR swept from 0 to 20 dB
    % Fixed at -5, 0, 5 deg
    % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
    param.monte.SNR   = repmat(linspace(10,10,1).', [1 1]);
    num_tests = size(param.monte.SNR,1);
    param.monte.DOA   = repmat([DOA],[num_tests 1]);
    param.monte.Nsnap = repmat(100,[num_tests 1]);
    param.monte.runs  = 100;
    param.monte.random_seed_offset = 0;
    
    %% Run the simulation
    results = sim.doa(param);
    
    %% Process and save the outputs
    
    RMSE = sim.doa_rmse(param,results);  % 4X1X1
    
    RMSE_Result(idx,:) = squeeze(RMSE);
  end
  %   end
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig4.7';
  
  figure(1); clf
%   plot([10:10:idx*10],squeeze(RMSE_Result),'+','LineWidth',2);
  plot([90:-10:20],squeeze(RMSE_Result),'+','LineWidth',2);
  grid on
  xlabel('Sector size [DOA \pm x]')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Sector size vs RMSE')
  
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end


%% Result-8: Sector size (i.e. limits) vs RMSE (two sources) [has a problem]
% Note: there is a bug in the initialization of methods 7 and 8 due to the
% overlap of the search ranges of the two sources, so fix the bug then run it.
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =100, SNR =[10], DoAs = [10 20]
%             fc =500MHz, BW =1*fc, W =3, Nbands =3 .
% =======================================================================
if 0
  %% Setup simulation parameters
  RMSE_Result = [];
  fc = 500e6;
  BW = 0.5*fc;
  DOA = [-8 8];
  for idx = 1:7
    param = [];
    
    % Source parameters
    param.src.fc                      = fc;
    param.src.f0                      = param.src.fc - BW/2;
    param.src.f1                      = param.src.fc + BW/2;
    param.src.ft_wind                 = @boxcar;
    param.src.lever_arm.fh            = @sim.lever_arm_example;
    param.src.lever_arm.args          = {[],1,[1:3],[0; c/param.src.fc/2; 0]};%0.24
    % DOA method parameters
    param.method.list                   = [2 7 8 9];
    param.method.Nsv                    = 24;
    param.method.OneD_Nsv               = 128;
%     param.method.src_limits             = {[DOA(1)-idx*10 DOA(1)+idx*10]/180*pi, [DOA(2)-idx*10 DOA(2)+idx*10]/180*pi};
    param.method.src_limits             = {[-40+5*(idx-1) 40-5*(idx-1)]/180*pi, [-40+5*(idx-1) 40-5*(idx-1)]/180*pi};
    param.method.theta_guard            = 1.5/180*pi;
    param.method.nb_nd.init             = 'ap';
    param.method.wb_td.init             = 'ap';
    param.method.wb_fd.init             = 'ap';
    param.method.wb_td.widening_factor = 1;
    param.method.wb_fd.filter_banks     = 3;
    
    % DOA monte carlo setup
    % Three equal power sources, SNR swept from 0 to 20 dB
    % Fixed at -5, 0, 5 deg
    % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
    param.monte.SNR   = repmat(linspace(10,10,1), [1 length(DOA)]);
    num_tests = size(param.monte.SNR,1);
    param.monte.DOA   = repmat(DOA,[num_tests 1]);
    param.monte.Nsnap = repmat(100,[num_tests 1]);
    param.monte.runs  = 100;
    param.monte.random_seed_offset = 0;
    
    %% Run the simulation
    results = sim.doa(param);
    
    %% Process and save the outputs
    
    RMSE = sim.doa_rmse(param,results);  % 4X1X1
    
    RMSE_Result(idx,:,:) = squeeze(RMSE);
  end
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig4.7';
  
  
  figure(1); clf
  plot([-40:5:-10],squeeze(RMSE_Result(:,:,1)),'+','LineWidth',2);
  grid on
  xlabel('Sector size [DOA \pm x]')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Sector size vs RMSE')
  
  figure(2); clf
  plot([-40:5:-10],squeeze(RMSE_Result(:,:,2)),'+','LineWidth',2);
  grid on
  xlabel('Sector size [DOA \pm x]')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Sector size vs RMSE')
  
%   figure(3); clf
%   plot([-40:5:-10],squeeze(RMSE_Result(:,:,3)),'+','LineWidth',2);
%   grid on
%   xlabel('Sector size [DOA \pm x]')
%   ylabel('RMS error ( \circ )')
%   legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
%   title('Sector size vs RMSE')
  
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig'])
  save([out_fn '.mat'],'param','results')
  
  return;
end

%% Result-9: Source separation vs RMSE
%
% Parameters: Nsrc =2, Nsens =8, Nsnaps =100, SNR =[10], DoAs = [0 variable]
%             fc =500MHz, BW =1*fc, W =3, Nbands =3 .
% =======================================================================
if 0
  %% Setup simulation parameters
  RMSE_Result = [];
  fc = 500e6;
  BW = 1.5*fc;
  DOA1 = [0];
  DOA2 = DOA1+[10:10:70];
  for idx = 1:length(DOA2)
    param = [];
    
    % Source parameters
    param.src.fc                      = fc;
    param.src.f0                      = param.src.fc - BW/2;
    param.src.f1                      = param.src.fc + BW/2;
    param.src.ft_wind                 = @boxcar;
    param.src.lever_arm.fh            = @sim.lever_arm_example;
    param.src.lever_arm.args          = {[],1,[1:8],[0; c/param.src.fc/2; 0]};%0.24
    % DOA method parameters
    param.method.list                   = [2 7 8 9];
    param.method.Nsv                    = 24;
    param.method.OneD_Nsv               = 128;
    param.method.src_limits             = {[DOA1-5 DOA1+4]/180*pi [DOA2(idx)-5 DOA2(idx)+5]/180*pi};
    param.method.theta_guard            = 1.5/180*pi;
    param.method.nb_nd.init             = 'ap';
    param.method.wb_td.init             = 'ap';
    param.method.wb_fd.init             = 'ap';
    param.method.wb_td.widening_factor = 3;
    param.method.wb_fd.filter_banks     = 3;
    
    % DOA monte carlo setup
    % Three equal power sources, SNR swept from 0 to 20 dB
    % Fixed at -5, 0, 5 deg
    % Snapshots fixed at 5 fast-time and 11 slow-time samples (i.e. 55)
    param.monte.SNR   = repmat(linspace(10,10,1), [1 2]);
    num_tests = size(param.monte.SNR,1);
    param.monte.DOA   = repmat([DOA1 DOA2(idx)],[num_tests 1]);
    param.monte.Nsnap = repmat(100,[num_tests 1]);
    param.monte.runs  = 500;
    param.monte.random_seed_offset = 0;
    
    %% Run the simulation
    results = sim.doa(param);
    
    %% Process and save the outputs
    
    RMSE = sim.doa_rmse(param,results);  % 4X1X1
    
    RMSE_Result(idx,:,:) = squeeze(RMSE);
  end
  out_fn_dir = 'E:\tmp\TSP_DOA';
  out_fn_name = 'tsp_doa_fig4.7';
  
  
  figure(1); clf
  plot([10:10:idx*10],squeeze(RMSE_Result(:,:,1)),'+','LineWidth',2);
  grid on
  xlabel('Source separation [deg]')
  ylabel('RMS error ( \circ )')
  legend({'MUSIC','MLE','WB','WBMLE','CRLB'})
  title('Source separation vs RMSE')
  
  
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  saveas(1,[out_fn '_src1.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end









