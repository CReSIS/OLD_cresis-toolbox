% script sim.doa_example
%
% Example setup scripts for the sim.doa function. To run, enable one
% of the if {0|1} blocks below and then run the script.
%
% Author: John Paden, Theresa Stumpf , Sravya Athinarapu
%
% See also: doa1.m

tic

clear, close;
physical_constants;

if 1
  % =======================================================================
  % Wax and Ziskind 1988 Fig 2
  % =======================================================================
  %% Setup simulation parameters
  
  
  param = [];
  Nc = 7;
  % Source parameters
  fc = 195e6;
  BW = 0.01e6;
  M = Nc-1; % max number of sources we are trying to estimate (it can go max upto Nc-1)
  
  param.Nc = Nc;
  param.M = M;
  param.src.f0                      = fc-BW/2;
  param.src.f1                      = fc+BW/2;
  param.src.ft_wind                 = @boxcar;
  param.src.lever_arm.fh            = @sim.lever_arm_example;
  param.src.lever_arm.args          = {[],1,[1:Nc],[0; c/fc/2; 0]};
  % DOA method parameters
  param.method.list                   = [7];
  param.method.Nsv                    = 3*Nc; % sampling grid for mle cost
  param.method.OneD_Nsv               = 128;
  %param.method.src_limits             = {[-20 40]/180*pi,[-20 40]/180*pi};
  %param.method.src_limits             =  repmat({[-20 40]/180*pi} , [1,Nc-1])
  param.method.src_limits             =  repmat({[-60 60]/180*pi} , [1,Nc-1])
  
  
  % since it Avoid evaluating cost function around previous sources and
  % hence search range decreases and can't be able to evaluate doa of more
  % sources case.
  
  param.method.theta_guard            = 0.5/180*pi; %1.5/180*pi;
  param.method.nb_nd.init             = 'ap';
  param.method.wb_td.init             = 'ap';
  param.method.wb_td.widening_factor  = 1;
  param.method.wb_fd.init             = 'ap';
  param.method.wb_fd.filter_banks     = 1;
  
  %param.monte.SNR   = repmat(linspace(10,25,16).' - 10*log10(3), [1 2]);
  
  
  %   DOA_true = [-34.8499045790465,-16.6015495990202,0,16.6015495990202];
  
  DOA_true =[0];  % SRAVYA
  
  if length(DOA_true) == 0
    param.monte.SNR   = -inf ;
    num_tests = 1;
    param.monte.DOA   = [];
  else
    param.monte.SNR   = repmat(10,1, length(DOA_true));  %% SNR 10dB used in the CHEN paper
    num_tests = size(param.monte.SNR,1);
    param.monte.DOA   = repmat(DOA_true,[num_tests 1]);
    
    
  end
  
  
  param.monte.Nsnap = repmat(100,[num_tests 1]);   %%% sample size 100 used in the CHEN paper
  param.monte.runs  = 10;      %% RUNS 10000 used in the CHEN paper
  param.monte.random_seed_offset = 0;
  
  % M : Maximum number of signals
  
  for Nsig_tmp = 1:M  %% possible range of k (number of signals)
    param.Nsig_tmp = Nsig_tmp;
    
    
    %% Run the simulation
    
    % SRAVYA
    if Nsig_tmp == 1
      [results{Nsig_tmp}, DCM_runs] = sim.doa1(param,[]);
    elseif Nsig_tmp > 1
      doa_prev = squeeze(results{1,Nsig_tmp-1}.theta_est{1,param.method.list});
      
      if param.monte.runs == 1
        [results{Nsig_tmp}, DCM_runs] = sim.doa1(param, doa_prev(1:Nsig_tmp-1));
        
      else
        [results{Nsig_tmp}, DCM_runs] = sim.doa1(param, doa_prev(1:param.monte.runs,1:Nsig_tmp-1));
      end
      
      
    end
    
    %% Process and save the outputs {method}(run_idx,test_idx,1:Nsig_tmp)
    out_fn_dir = 'D:\tmp\TSP_DOA';
    out_fn_name = 'wax_ziskind_fig2';
    
    RMSE{Nsig_tmp} = sim.doa_rmse(param,results{Nsig_tmp});
    
    % %   figure(1); clf
    % %   plot(param.monte.SNR(:,1)+10*log10(3),RMSE(:,:,1).','.','LineWidth',2);
    % %   grid on
    % %   xlabel('Source SNR (dB)')
    % %   ylabel('RMS error ( \circ )')
    % %   legend({'MUSIC','MLE','WB','WBMLE'})
    % %   title('Figure 2 from Wax and Ziskind 1988')
    % %
    % %   figure(2); clf
    % %   plot(param.monte.SNR(:,2)+10*log10(3),RMSE(:,:,2).','.','LineWidth',2);
    % %   grid on
    % %   xlabel('Source SNR (dB)')
    % %   ylabel('RMS error ( \circ )')
    % %   legend({'MUSIC','MLE','WB','WBMLE'})
    % %   title('Complement to figure 2 from Wax and Ziskind (source 2)')
    % %
    % %   % Save Outputs
    % %   if ~exist(out_fn_dir,'dir')
    % %     mkdir(out_fn_dir);
    % %   end
    % %   out_fn = fullfile(out_fn_dir,out_fn_name);
    % %   saveas(1,[out_fn '_src1.fig']);
    % %   saveas(2,[out_fn '_src2.fig']);
    % %   save([out_fn '.mat'],'param','results')
    
    doa_tmp(1:param.monte.runs,:) = squeeze(results{1, Nsig_tmp}.theta_est{1, param.method.list});
    
    doa_mle_all{Nsig_tmp} = doa_tmp(:,1:Nsig_tmp) *180/pi; % each cell contain DOA for all runs. row indicates the run number.
    
    
  end
  
  
  bin = param.monte.runs;  % for now
  lineIdx = 1;
  
  
  for binIdx = 1:bin;
    
    %arranging all doa for possible k for a single run.
    for Nsig_tmp = 1:M  %%%%%CHANGE
      
      doa_mle{Nsig_tmp} = doa_mle_all{Nsig_tmp}(binIdx,:) ;
      
    end
    
    
    Rxx = DCM_runs{binIdx};
    
    
    % l: eigen values
    [l,index] = sort(eig(Rxx),'descend');
    
    
    
    % DCM is symmetric and always we get real eigen values. due to some
    % rounding errors in the data generated by matlab we got complex eigen
    % values
    
    l = real(l);
    l_all(bin,:) = l;
    [V,D] = eig(Rxx);
    %u: eigen vectors
    u = V(:,index);
    %normalizing the eigen values
    l_norm = l./min(l);               % MIGHT CHANGE THE WAY OF NORMALIZING
    l_norm_all(bin,:) = l_norm;
    
    
    
    
    %% Model Order Estimators
    
    %Optimal Methods
    
    for model_order_method = 0:7
      
      clear sources_number doa
      [sources_number,doa] = sim.model_order_optimal(doa_mle,Nc,Rxx,param.monte.Nsnap,binIdx, model_order_method);
      
      [doa,sort_idxs] = sort(doa);
      
      switch model_order_method
        case 0
          
          model_order_results.MDL_sravya.doa(binIdx,:,lineIdx)  = doa;
          model_order_results.MDL_sravya.Nest(binIdx,lineIdx)  = sources_number;
          
          
        case 1
          
          model_order_results.MDL.doa(binIdx,:,lineIdx)  = doa;
          model_order_results.MDL.Nest(binIdx,lineIdx)  = sources_number;
          
        case 2
          
          model_order_results.AIC.doa(binIdx,:,lineIdx)  = doa;
          model_order_results.AIC.Nest(binIdx,lineIdx)  = sources_number;
          
        case 3
          model_order_results.BIC.doa(binIdx,:,lineIdx)  = doa;
          model_order_results.BIC.Nest(binIdx,lineIdx)  = sources_number;
          
        case 4
          
          model_order_results.HQ.doa(binIdx,:,lineIdx)  = doa;
          model_order_results.HQ.Nest(binIdx,lineIdx)  = sources_number;
          
        case 5
          model_order_results.AICc.doa(binIdx,:,lineIdx)  = doa;
          model_order_results.AICc.Nest(binIdx,lineIdx)  = sources_number;
          
        case 6
          
          model_order_results.KICvc.doa(binIdx,:,lineIdx)  = doa;
          model_order_results.KICvc.Nest(binIdx,lineIdx)  = sources_number;
          
        case 7
          model_order_results.WIC.doa(binIdx,:,lineIdx)  = doa;
          model_order_results.WIC.Nest(binIdx,lineIdx)  = sources_number;
          
          
          
          
        otherwise
          error('Not supported')
      end
      
    end
    %% Suboptimal Methods
    
    for model_order_method = 0:7
      
      clear sources_number doa
      sources_number = sim.model_order_suboptimal(Nc,Rxx,param.monte.Nsnap,binIdx, model_order_method);
      
      doa = NaN *ones(1,Nc-1);
      
      %
      if sources_number > 0
        doa(1:sources_number) = doa_mle{sources_number};
      end
      
      [doa,sort_idxs] = sort(doa);
      
      
      switch model_order_method
        case 0
          
          model_order_results_suboptimal.MDL_sravya.doa(binIdx,:,lineIdx)  = doa;
          model_order_results_suboptimal.MDL_sravya.Nest(binIdx,lineIdx)  = sources_number;
          
          
        case 1
          
          model_order_results_suboptimal.MDL.doa(binIdx,:,lineIdx)  = doa;
          model_order_results_suboptimal.MDL.Nest(binIdx,lineIdx)  = sources_number;
          
        case 2
          
          model_order_results_suboptimal.AIC.doa(binIdx,:,lineIdx)  = doa;
          model_order_results_suboptimal.AIC.Nest(binIdx,lineIdx)  = sources_number;
          
        case 3
          
          model_order_results_suboptimal.BIC.doa(binIdx,:,lineIdx)  = doa;
          model_order_results_suboptimal.BIC.Nest(binIdx,lineIdx)  = sources_number;
        case 4
          
          model_order_results_suboptimal.HQ.doa(binIdx,:,lineIdx)  = doa;
          model_order_results_suboptimal.HQ.Nest(binIdx,lineIdx)  = sources_number;
          
        case 5
          model_order_results_suboptimal.AICc.doa(binIdx,:,lineIdx)  = doa;
          model_order_results_suboptimal.AICc.Nest(binIdx,lineIdx)  = sources_number;
          
        case 6
          
          
          model_order_results_suboptimal.KICvc.doa(binIdx,:,lineIdx)  = doa;
          model_order_results_suboptimal.KICvc.Nest(binIdx,lineIdx)  = sources_number;
          
        case 7
          
          model_order_results_suboptimal.WIC.doa(binIdx,:,lineIdx)  = doa;
          model_order_results_suboptimal.WIC.Nest(binIdx,lineIdx)  = sources_number;
          
          
          
          
        otherwise
          error('Not supported')
      end
      
      
    end
    
  end
  
  
  %% Saving results from all model order estimators for comparision
  
  Nest_results(:,1) = length(DOA_true).*ones(bin,1);
  Nest_results(:,2) = model_order_results_suboptimal.MDL.Nest;
  Nest_results(:,3) = model_order_results_suboptimal.AIC.Nest;
  Nest_results(:,4) = model_order_results_suboptimal.BIC.Nest;
  Nest_results(:,5) = model_order_results_suboptimal.HQ.Nest;
  Nest_results(:,6) = model_order_results_suboptimal.AICc.Nest;
  Nest_results(:,7) = model_order_results_suboptimal.KICvc.Nest;
  Nest_results(:,8) = model_order_results_suboptimal.WIC.Nest;
  Nest_results(:,9) = model_order_results_suboptimal.MDL_sravya.Nest;
  
  Nest_results(:,10) = model_order_results.MDL.Nest;
  Nest_results(:,11) = model_order_results.AIC.Nest;
  Nest_results(:,12) = model_order_results.BIC.Nest;
  Nest_results(:,13) = model_order_results.HQ.Nest;
  Nest_results(:,14) = model_order_results.AICc.Nest;
  Nest_results(:,15) = model_order_results.KICvc.Nest;
  Nest_results(:,16) = model_order_results.WIC.Nest;
  Nest_results(:,17) = model_order_results.MDL_sravya.Nest;
  
  
  for i= 1:17
    correct_estimation(i) =  length(find(Nest_results(:,i)==length(DOA_true))) ;
  end
  
  percentage_correct_q =  (correct_estimation(2:17)/correct_estimation(1))*100;
  method_title =   {'MDL  ', 'AIC  ', 'BIC  ', 'HQ   ', 'AICc ', 'KICvc', 'WIC  ','MDL S'}
  
  figure(8),clf;
  suptitle (['Comparision for number of sources =' num2str(length(DOA_true))'']);
  subplot(1,2,1)
  stem(percentage_correct_q(1:8),'LineWidth',10,'Marker','.')
  set(gca,'xticklabel',method_title)
  title('Suboptimal Model Order Estimators')
  ylim([0 100])
  grid on
  xlabel('model order estimators')
  ylabel('percentage of estimating the number of sources correctly')
  
  
  subplot(1,2,2)
  stem(percentage_correct_q(9:16),'LineWidth',10,'Marker','.')
  set(gca,'xticklabel',method_title)
  title('Optimal Model Order Estimators')
  
  xlabel('model order estimators')
  ylabel('percentage of estimating the number of sources correctly')
  
  
  ylim([0 100])
  grid on
  
  
  % COMPARISION TABLE
  comparision = strcat( ['# q   =';'MDL   ='; 'AIC   ='; 'BIC   ='; 'HQ    ='; 'AICc  ='; 'KICvc ='; 'WIC   =';'MDLs  =';'MDL   ='; 'AIC   ='; 'BIC   ='; 'HQ    ='; 'AICc  ='; 'KICvc ='; 'WIC   =';'MDLs  ='] ,num2str(Nest_results.')) ;
  
  
  
  toc
  
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
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
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
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

if 0
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
  param.src.lever_arm.args          = {[],1,[1:7],[0; c/fc/2; 0]};
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
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  saveas(3,[out_fn '_src3.fig']);
  saveas(4,[out_fn '_src4.fig']);
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
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  saveas(3,[out_fn '_src3.fig']);
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
  param.src.lever_arm.args          = {[],1,[1:8],[0; 0.24; 0]};
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
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
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
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end

if 0
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
  saveas(1,[out_fn '_src1.fig']);
  saveas(2,[out_fn '_src2.fig']);
  saveas(3,[out_fn '_src3.fig']);
  saveas(4,[out_fn '_src4.fig']);
  save([out_fn '.mat'],'param','results')
  
  return;
end
