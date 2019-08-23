% script sim.doa_example_moe
%
% Example setup scripts for the sim.doa function. This includes tutorial
% examples that illustrate how to use the model order estimation.
%
% Author: John Paden, Theresa Stumpf
%
% See also: doa.m

%%  Model Order Estimation
% =========================================================================

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
opt_norm = ~norm_saved;

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
          [results{Nsig_tmp}, DCM_runs] = sim.doa(param);
          
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
      % for each k we should run sim.doa.m whereas only suboptimal methods
      % are required just k= 1 is enough to generate data and as we do not
      % need doa estimation(no need to run sim.doa.m multiple times)
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
      [results{Nsig_tmp}, DCM_runs] = sim.doa(param);
      
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
      % Eigenvalues for all the SNRs for a given Nsrc
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
        % Eigenvalues for all the SNRs for a given Nsrc
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
