function [results, DCM_runs] = doa(param)
% results = doa(param)
%
% Function for simulating direction of arrival algorithms
%
% Author: John Paden, Theresa Stumpf

param.src.fc = (param.src.f0 + param.src.f1)/2;
param.src.fs = param.src.f1 - param.src.f0;

% Get phase center positions of antenna array
if isfield(param.src,'phase_center')
  param.src.y_pc   = param.src.phase_center(2,:).';
  param.src.z_pc   = param.src.phase_center(3,:).';
  
elseif isfield(param.src,'lever_arm') && isfield(param.src.lever_arm,'fh')
  [phase_center] = param.src.lever_arm.fh(param.src.lever_arm.args{:});
  param.src.y_pc   = phase_center(2,:).';
  param.src.z_pc   = phase_center(3,:).';
elseif isfield(param.src,'y_pc') && isfield(param.src,'z_pc')
  % Do nothing. Already defined 
else
  warning('Phase center information are not defined')
  keyboard
end

% Setup theta grid and steering vectors
[param.method.theta, param.method.SV] ...
  = array_proc_sv(param.method.Nsv,param.src.fc, param.src.y_pc, param.src.z_pc);
param.method.theta             = fftshift(param.method.theta);
param.method.SV                = fftshift(param.method.SV,2);

[param.method.OneD_theta, param.method.OneD_SV] ...
  = array_proc_sv(param.method.OneD_Nsv,param.src.fc, param.src.y_pc, param.src.z_pc);
param.method.OneD_theta             = fftshift(param.method.OneD_theta);
param.method.OneD_SV                = fftshift(param.method.OneD_SV,2);

% Setup parameterization structure for narrowband 1-D search
% methods (e.g. MUSIC)
doa_nb_1d_param = [];
doa_nb_1d_param.fs          = param.src.fs;
doa_nb_1d_param.fc          = param.src.fc;
doa_nb_1d_param.src_limits  = param.method.src_limits;
doa_nb_1d_param.theta_guard = param.method.theta_guard;
doa_nb_1d_param.y_pc        = param.src.y_pc;
doa_nb_1d_param.z_pc        = param.src.z_pc;
doa_nb_1d_param.theta       = param.method.OneD_theta;
doa_nb_1d_param.SV          = param.method.OneD_SV;
doa_nb_1d_param.options     = optimoptions(@fmincon,'Display','off','Algorithm','sqp','TolX',1e-3);

% Setup parameterization structure for narrowband N-D search
% methods (e.g. MLE)
doa_nb_nd_param = [];
doa_nb_nd_param.fs          = param.src.fs;
doa_nb_nd_param.fc          = param.src.fc;
doa_nb_nd_param.src_limits  = param.method.src_limits;
doa_nb_nd_param.theta_guard = param.method.theta_guard;
doa_nb_nd_param.y_pc        = param.src.y_pc;
doa_nb_nd_param.z_pc        = param.src.z_pc;
doa_nb_nd_param.theta       = param.method.theta;
doa_nb_nd_param.SV          = param.method.SV;
doa_nb_nd_param.options     = optimoptions(@fmincon,'Display','off','Algorithm','sqp','TolX',1e-3);
doa_nb_nd_param.search_type = param.method.nb_nd.init;

% Setup wideband parameterization structure for wideband time domain
% search methods (e.g. WB)
doa_wb_td_param= [];
doa_wb_td_param.fs          = param.src.fs;
doa_wb_td_param.fc          = param.src.fc;
doa_wb_td_param.src_limits  = param.method.src_limits;
doa_wb_td_param.theta_guard = param.method.theta_guard;
doa_wb_td_param.y_pc        = param.src.y_pc;
doa_wb_td_param.z_pc        = param.src.z_pc;
doa_wb_td_param.theta       = param.method.theta;
doa_wb_td_param.SV          = param.method.SV;
doa_wb_td_param.options     = optimoptions(@fmincon,'Display','off','Algorithm','sqp','TolX',1e-3);
doa_wb_td_param.search_type = param.method.wb_td.init;

% Setup wideband parameterization structure for wideband frequency domain
% search methods (e.g. MLEWB)
doa_wb_fd_param= [];
doa_wb_fd_param.fs          = param.src.fs;
doa_wb_fd_param.fc          = param.src.fc;
doa_wb_fd_param.src_limits  = param.method.src_limits;
doa_wb_fd_param.theta_guard = param.method.theta_guard;
doa_wb_fd_param.y_pc        = param.src.y_pc;
doa_wb_fd_param.z_pc        = param.src.z_pc;
doa_wb_fd_param.theta       = param.method.theta;
doa_wb_fd_param.SV          = param.method.SV;
doa_wb_fd_param.options     = optimoptions(@fmincon,'Display','off','Algorithm','sqp','TolX',1e-3);
doa_wb_fd_param.search_type = param.method.wb_fd.init;
doa_wb_fd_param.nb_filter_banks = param.method.wb_fd.filter_banks;

%% Monte Carlo Trials
%==========================================================================

theta_est = [];
if isfield(param,'M') && ~isempty(param.M)
    % Model order estimation simulation. M=Nc-1 usually.
    for method = param.method.list
        theta_est{method} = zeros(param.monte.runs,size(param.monte.SNR,1),param.Nc-1);
        hessian_est{method} = zeros(param.monte.runs,size(param.monte.SNR,1),param.Nc-1);
        amp_est{method} = zeros(param.monte.runs,size(param.monte.SNR,1),param.Nc-1);
        cost_func{method} = zeros(param.monte.runs,size(param.monte.SNR,1));
    end
else
    for method = param.method.list
        theta_est{method} = zeros(param.monte.runs,size(param.monte.SNR,1),size(param.monte.SNR,2));
        hessian_est{method} = zeros(param.monte.runs,size(param.monte.SNR,1),size(param.monte.SNR,2));
        amp_est{method} = zeros(param.monte.runs,size(param.monte.SNR,1),size(param.monte.SNR,2));
        cost_func{method} = zeros(param.monte.runs,size(param.monte.SNR,1));
    end
end

rng_args = zeros(param.monte.runs);

doa_nonlcon_fh = eval(sprintf('@(x) doa_nonlcon(x,%f);', param.method.theta_guard));

start_time = tic;
for test_idx = 1:size(param.monte.SNR,1)
  
    if isempty(param.monte.DOA)   % (for k = 0 model order estimation)  
        param.src.DOAs   = []; 
    else
        param.src.DOAs   = param.monte.DOA(test_idx,:);
    end
  param.src.SNR    = param.monte.SNR(test_idx,:);
  param.src.Nsnap  = param.monte.Nsnap(test_idx);
  
  % Set number of signals, Nsig, field
  if isfield(param,'Nsig_tmp') && ~isempty(param.Nsig_tmp)
      % For model order estimation simulation.
      doa_nb_1d_param.Nsig = param.Nsig_tmp;
      doa_nb_nd_param.Nsig = param.Nsig_tmp;
      doa_wb_td_param.Nsig = param.Nsig_tmp;
      doa_wb_fd_param.Nsig = param.Nsig_tmp;
      
      LB = zeros(param.Nsig_tmp,1);
      UB = zeros(param.Nsig_tmp,1);
  else
      doa_nb_1d_param.Nsig = size(param.src.SNR,2);
      doa_nb_nd_param.Nsig = size(param.src.SNR,2);
      doa_wb_td_param.Nsig = size(param.src.SNR,2);
      doa_wb_fd_param.Nsig = size(param.src.SNR,2);
      
      LB = zeros(length(param.src.DOAs),1);
      UB = zeros(length(param.src.DOAs),1);
  end
  
  % Set source limits for N-dimensional constrained optimization
  for src_idx = 1:length(LB)
    LB(src_idx) = param.method.src_limits{src_idx}(1);
    UB(src_idx) = param.method.src_limits{src_idx}(2);
  end
  
  for run_idx = 1:param.monte.runs
      if isfield(param,'Nsig_tmp') && ~isempty(param.Nsig_tmp)
          % For model order estimation simulation.
          fprintf('q:%2d  k:%2d of %d / run %4d of %d (%.1f sec)\n', ...
              length(param.monte.DOA),param.Nsig_tmp, param.M, run_idx, param.monte.runs, toc(start_time));
      else
          if ~mod(run_idx-1,50)
              fprintf('test: %2d of %d / run %4d of %d (%.1f sec)\n', ...
                  test_idx, size(param.monte.SNR,1), run_idx, param.monte.runs, toc(start_time));
          end
      end
    % Setup random number generator
    rng_args(run_idx) = param.monte.random_seed_offset + run_idx;
    rng(rng_args(run_idx));
    
    % Simulate array data
    [Data,DCM,imp_response,DCM_fd] = sim.doa_wideband_data(param);
    
    if isfield(param,'Nsig_tmp') && ~isempty(param.Nsig_tmp)
        % For model order estimation simulation.
        DCM_runs{run_idx} = DCM;
    end
  
    % Set up estimation parameters for each method
    doa_wb_td_param.h    = conv(imp_response.vals,imp_response.vals,'same');
    doa_wb_td_param.h    = doa_wb_td_param.h ./ max(abs(doa_wb_td_param.h));
    doa_wb_td_param.t0   = imp_response.time_vec(1);
    doa_wb_td_param.dt   = mean(diff(imp_response.time_vec));
    doa_wb_td_param.DCM  = DCM;
    
    DCM_nb_idxs = (param.method.wb_td.widening_factor-1)/2*length(param.src.y_pc) + (1:length(param.src.y_pc));
    DCM_nb = DCM(DCM_nb_idxs,DCM_nb_idxs);
    doa_nb_1d_param.Rxx = DCM_nb;
    doa_nb_nd_param.Rxx = DCM_nb;
    
    doa_wb_fd_param.Rxx = DCM_fd;
    
    if 0
      figure;plot(imp_response.time_vec,imp_response.vals);
      hold on
      plot(imp_response.time_vec,doa_param.h,'r')
      keyboard
    end 
    
    for method = param.method.list
      %fprintf('test: %2d / run: %2d / method: %2d \n', test_idx, run_idx, method);
      switch method
        case 2
          % MUSIC method: init and minimization
          doa0 = sort(music_initialization(DCM_nb,doa_nb_1d_param));
          
          [doa,Jval,exitflag,OUTPUT,~,~,HESSIAN] = ...
            fmincon(@(theta_hat) music_cost_function(theta_hat,doa_nb_1d_param), doa0,[],[],[],[],LB,UB,doa_nonlcon_fh,doa_nb_1d_param.options);
          
        case 7
          % MLE method: init and minimization
          doa0 = sort(mle_initialization(DCM_nb,doa_nb_nd_param));
          
          [doa,Jval,exitflag,OUTPUT,~,~,HESSIAN] = ...
            fmincon(@(theta_hat) mle_cost_function(theta_hat,doa_nb_nd_param), doa0,[],[],[],[],LB,UB,doa_nonlcon_fh,doa_nb_nd_param.options);
          
        case 8 % WB method
          % WB method: init and minimization
          doa0 = sort(wb_initialization(DCM,doa_wb_td_param));
          
          [doa,Jval,exitflag,OUTPUT,~,~,HESSIAN] = ...
            fmincon(@(theta_hat) wb_cost_function(theta_hat,doa_wb_td_param), doa0,[],[],[],[],LB,UB,doa_nonlcon_fh,doa_wb_td_param.options);
          
        case 9
          % WBMLE method: init and minimization
          doa0 = sort(wbmle_initialization(DCM_fd,doa_wb_fd_param));
          
          [doa,Jval,exitflag,OUTPUT,~,~,HESSIAN] = ...
            fmincon(@(theta_hat) wbmle_cost_function(theta_hat,doa_wb_fd_param), doa0,[],[],[],[],LB,UB,doa_nonlcon_fh,doa_wb_fd_param.options);
      end
      
      % Store outputs into variables
      if isfield(param,'Nsig_tmp') && ~isempty(param.Nsig_tmp)
          % For model order estimation simulation.
          [theta_est{method}(run_idx,test_idx,1:param.Nsig_tmp),sort_idxs] = sort(doa);
          HESSIAN = diag(HESSIAN);
          
          % for only suboptimal methods (param.subopt_only). DOa estimation not required.
          %In that case run using doa_example_suboptimal.m and the following
          %section cannot be evaluated as we do not have doa.
          
          if param.doa_example == 1
              hessian_est{method}(run_idx,test_idx,1:param.Nsig_tmp) = HESSIAN(sort_idxs);
              cost_func{method}(run_idx,test_idx) = Jval;
          end
      else
          [theta_est{method}(run_idx,test_idx,:),sort_idxs] = sort(doa);
          HESSIAN = diag(HESSIAN);
          hessian_est{method}(run_idx,test_idx,:) = HESSIAN(sort_idxs);
          cost_func{method}(run_idx,test_idx) = Jval;
      end
    end
  end
end

% Copy outputs into output argument structure
if isfield(param,'Nsig_tmp') && ~isempty(param.Nsig_tmp)
    % For model order estimation simulation.
    results.theta_est = theta_est;
    results.rng_args = rng_args;
    
    if param.doa_example == 1
        results.hessian_est = hessian_est;
        results.cost_func = cost_func;
    end
    results.DCM = DCM_runs;
else
    results.theta_est = theta_est;
    results.rng_args = rng_args;
    results.hessian_est = hessian_est;
    results.cost_func = cost_func;
end

return;
