%% Data loading and setup
% nb = load('Y:\mohanad\Results\MOE\DCMandEigenvalues\1D_NB.mat');
data_training = load('/users/mohanad/IceSheetProject/MOE work/DCMandEigenvalues/1D_WB_LL_eigenval_10000runs.mat');
eigenvalues_training = data_training.eigenvalues_runs_snr; % Each cell is Nruns by Ntargets by Nb * Nc
LL_training         = data_training.loglikelihood_runs_snr;
LL_training_opt     = LL_training.opt;    % Each cell is Nruns by Ntargets by Ntargets
LL_training_subopt  = LL_training.subopt; % Each cell is Nruns by Ntargets by Ntargets

% m: number of training examples. Don't use the same training dataset for
% testing. Split the training examples into 70 percent for training and 30
% percent for testing.
m = size(eigenvalues_training{1},1);
SampleSize = 5000;%round(logspace(log10(10),log10(m),5));
for loop_i = 1:length(SampleSize)
  clearvars -except SampleSize m loop_i percentage_correct_all data_training eigenvalues_training J_testing J_training Jval
  h1 = sprintf('\nLoop # %d of %d...\n',loop_i,length(SampleSize));
  disp(h1)
  % rng(loop_i)
  
  m_training = SampleSize(loop_i);
  m_testing  = floor(m - m_training);
  
  Nc = 7;
  
  % Choos input type: 'loglikelihoods', 'geometric/arithmetic means'
  % (geometric and arithmetic means of eigenvalues), or 'eigenvalues'.
  % input_type = 'eigenvalues';
  % input_type = 'loglikelihoods';
  input_type = 'geometric/arithmetic means';
  
  % training_flag/testing_flag: boolean variable to turn on/off training or testing parts/modes
  training_flag  = 1;
  testing_flag   = 1;
  optimizer_flag = 1;
  
  %% Training (apply logistic regression to the training dataset)
  if training_flag
    % Geometric and arithmetic mean of the eigenvalues
    if strcmp(input_type,'geometric/arithmetic means')
      clear geom_mean arithm_mean
      for snr_i = 1:length(data_training.eigenvalues_runs_snr)
        for run_i = 1:m_training;%size(data_training.eigenvalues_runs_snr{1},1)
          for band_i = 1:size(data_training.eigenvalues_runs_snr{1},3)
            for q_i = 1:size(data_training.eigenvalues_runs_snr{1},2)
              eigval = squeeze(data_training.eigenvalues_runs_snr{snr_i}(run_i,q_i,band_i,:));
              for k_i = 0:size(data_training.eigenvalues_runs_snr{1},4)-1
                geom_mean_training{snr_i}(run_i,q_i,band_i,k_i+1)   = sum(log(eigval(k_i+1:end)));
                arithm_mean_training{snr_i}(run_i,q_i,band_i,k_i+1) = (Nc-k_i)*(log(sum((eigval(k_i+1:end))/(Nc-k_i))));
                
                %     log_lq(nb_idx,k+1) = Nsnap*(geom_mean- arithm_mean);
              end
            end
          end
        end
      end
    end
    
    % max_Ntargets: maximum expected/required number of targets
    % (q=0:max_Ntargets)
    max_Ntargets = 6;
    
    % Nexamples_training_q: number of training examples per specific number of targets (Q)
    Nexamples_training_q = m_training*length(eigenvalues_training);
    Nexamples = Nexamples_training_q*(max_Ntargets+1);
    training_input = [];
    for q_i = 1 + [0:max_Ntargets]
      % Concatinate the input data from all SNRs for this number of sources
      for snr_i = 1:length(eigenvalues_training)
        if strcmp(input_type,'eigenvalues')
          input_snr = eigenvalues_training{snr_i};
        elseif strcmp(input_type,'loglikelihoods')
          % Loglikelihoods were stored as -2L. So, here we divide by -2.
          input_snr = LL_training_opt{snr_i}*(-1/2);
          %       input_snr = exp(LL_training_subopt{snr_i}*(-1/2));
        elseif strcmp(input_type,'geometric/arithmetic means')
          input_snr1 = geom_mean_training{snr_i}(1:m_training,q_i,:,:);
          input_snr2 = arithm_mean_training{snr_i}(1:m_training,q_i,:,:);
          %       input_snr3 = eigenvalues_training{snr_i}(1:m_training,:,:,:);
          %       input_snr3 = input_snr1./input_snr2;
          %       input_snr = cat(4,input_snr1,input_snr2,input_snr3);
          input_snr = cat(4,input_snr1,input_snr2);
        else
          warning('Unsupported training-input data type.')
          keyboard
        end
        %     training_input(end+1:end+m_training,:) = mean(real(squeeze(input_snr(1:m_training,q_i,:,:))),2);
%         training_input(end+1:end+m_training,:) = real(squeeze(input_snr(1:m_training,q_i,:,:)));
        training_input(end+1:end+m_training,:) = real(squeeze(input_snr));
      end
      %       for q_i = 1 + [0:max_Ntargets]
      %         training_input(end+1:end+m_training,:) = real(squeeze(input_snr(1:m_training,q_i,:,:)));
    end
    
    % Handle the case of bad examples (i.e. examples that have NaN or Inf entries)
    bad_examples_idxs = [];
    for bad_example_i = 1:size(training_input,1)
      if any(isnan(training_input(bad_example_i,:))) || any(isinf(training_input(bad_example_i,:)))
        bad_examples_idxs(end+1) = bad_example_i;
      end
    end
    training_input(bad_examples_idxs,:) = min(training_input(:));% [];
    
    % Normalize the input data so that it is insensitive to scaling. Choosing
    % the right option is input value-dependent (trial and error type of scenario).
    if 0
      % Normalize wrt the maximum value of each example
      %  --- works best for eigenvalues ---
      training_input = training_input./repmat(max(abs(training_input),[],2),[1 size(training_input,2)]);
    elseif 0
      % Subtract the mean and normalize wrt the maximum value of each example
      training_input = (training_input-repmat(mean(training_input,2),[1 size(training_input,2)]))./repmat(max(training_input,[],2),[1 size(training_input,2)]);
    elseif 1
      % Subtract the mean and normalize wrt the std value of each example
      training_input = (training_input-repmat(mean(training_input,2),[1 size(training_input,2)]))./repmat(std(training_input,[],2),[1 size(training_input,2)]);
    elseif 0
      % Subtract the mean and normalize wrt the std value of each example
      %  --- works for arithmetic and geometric means ---
      training_input = (training_input-repmat(mean(training_input,1),[Nexamples 1]))./repmat(std(training_input,[],1),[Nexamples 1]);
      %   training_input = training_input./repmat(max(abs(training_input),[],2),[1 Nc]);
    end
    
    % Make training_input an Nc by Nexamples matrix
    training_input = training_input.';
    % training_input = [ones(1,size(training_input,2)); training_input];
  end
  %% Model and initialization
  % h: Model of the in-out relationship. It is the most important piece.
  % Watch the size of the parameter theta0 (initial value of theta).
  if 0
    % This model works best for NB eigenvalues
    h = @(theta,x) (1./(1+exp(-(theta(1:size(x,1)).'*x.^(-1/2) + theta(size(x,1)+1:2*size(x,1)).'*x.^(-1/3)))));
    theta0 = zeros(2*size(training_input,1),1);
    %   theta0 = rand(2*size(training_input,1),1)*(2*1e-2)-1e-2;
  elseif 1
    % Best for WB geometric and arithmetic mean (concatenated)
    h = @(theta,x) (1./(1+exp(-(theta(1:Nc).'*x(1:Nc,:).^(1) + theta(Nc+1:2*Nc).'*x(Nc+1:end,:).^(1)))));
%     theta0 = zeros(2*Nc,1);
    %   d_theta = 1e0;
    %   theta0 = 2*d_theta*rand(2*Nc,1)-d_theta;
    %   theta0_all(:,loop_i) = theta0;
  elseif 0
    h = @(theta,x) (1./(1+exp(-(theta(1)*x(1,:) + theta(2:Nc+1).'*x(2:Nc+1,:).^(1) + theta(Nc+2:2*Nc+1).'*x(Nc+2:end,:).^(1)))));
    theta0 = zeros(2*Nc+1,1);
  elseif 0
    h = @(theta,x) (1./(1+exp(-(theta(1:Nc).'*x(1:Nc,:).^(1) + theta(Nc+1:2*Nc).'*x(Nc+1:2*Nc,:).^(1)))) + ...
      theta(2*Nc+1:3*Nc).'*x(2*Nc+1:3*Nc,:).^(1));
    theta0 = zeros(3*Nc,1);
  elseif 0
    h = @(theta,x) (1./(1+exp(-(theta(1:Nc).'*x(1:Nc,:).^(1) + theta(Nc+1:2*Nc).'*x(Nc+1:2*Nc,:).^(1)))) + ...
      theta(2*Nc+1:3*Nc).'*x(2*Nc+1:3*Nc,:).^(-1/2) );
    theta0 = zeros(3*Nc,1);
  elseif 0
    h = @(theta,x) (1./(1+exp(-(theta(1:Nc).'*x(1:Nc,:).^(1) + theta(Nc+1:2*Nc).'*x(Nc+1:2*Nc,:).^(1)))) + ...
      theta(2*Nc+1:3*Nc).'*x(1:Nc,:).^(-1/2) + theta(3*Nc+1:4*Nc).'*x(Nc+1:2*Nc,:).^(-1/2));
    theta0 = zeros(4*Nc,1);
  elseif 0
    h = @(theta,x) (1./(1+exp(-(theta(1:size(x,1)).'*x.^(-1/2) + theta(size(x,1)+1:2*size(x,1)).'*x.^(-1/3) + ...
      theta(2*size(x,1)+1:3*size(x,1)).'*x.^(-1/4)))));
    theta0 = zeros(3*size(training_input,1),1);
  elseif 0
    h = @(theta,x) (1./(1+exp(-(theta(1:size(x,1)).'*x + theta(size(x,1)+1:2*size(x,1)).'*x.^(1/2)))));
    theta0 = zeros(2*size(training_input,1),1);
  elseif 0
    h = @(theta,x) (1./(1+exp(-(theta.'*x.^(-1/2)))));
    theta0 = zeros(size(training_input,1),1);
  elseif 0
    h = @(theta,x) (1./(1+exp(-(theta(1:size(x,1)).'*x + theta(size(x,1)+1:2*size(x,1)).'*x.^(1/2) + ...
      theta(2*size(x,1)+1:3*size(x,1)).'*x.^(3)))));
    theta0 = zeros(3*size(training_input,1),1);
  end
  
  if 0
    % In case needed
    training_label_all = zeros(size(training_input,2),1);
    for class_i = 1:7
      training_label_all((class_i-1)*Nexamples_training_q+1:class_i*Nexamples_training_q) = class_i;
    end
  end
  
  %% Optimizer setup
  if optimizer_flag
%     reg_param = [0 0 100 0  0  0   1400];
    reg_param = [0:100:5000];
%     reg_param = [0     0     0     0   400   300     0];
%     reg_param = [0 100 900 500 200 700 600];
    
    lb = -inf(size(training_input,1),1);
    ub = +inf(size(training_input,1),1);
%     lb = -inf(length(theta0),1);
%     ub = +inf(length(theta0),1);
    
    options =  optimoptions(@fmincon,'TolX',1e-6,'TolFun',1e-6,'MaxIter',10^7,'MaxFunEvals',10^7,'Display','off');%, ...
%       'PlotFcn',{@optimplotfval,@optimplotstepsize});%,'Algorithm','sqp');
    for reg_param_i = 1:length(reg_param)
      h2 = sprintf('\n   lambda # %d of %d ...\n',reg_param_i,length(reg_param));
      disp(h2)
%     clear theta Jval exitflag
    training_label_all = [];
    for class_i = 1+[0:max_Ntargets]
      % Output training data (or class ID): Each time we set all outputs to 0
      % except for one class, which we set to class_i.
      training_label = zeros(size(training_input,2),1);
      training_label((class_i-1)*Nexamples_training_q+1:class_i*Nexamples_training_q) = class_i;
      
      % Cost function parameters (LR stands for Logistic Regression)
      LR_cost_fn_param.training_input  = training_input;
      LR_cost_fn_param.training_output = training_label;
      LR_cost_fn_param.model           = h;
%       LR_cost_fn_param.reg_param       = reg_param(class_i);
      LR_cost_fn_param.reg_param       = reg_param(reg_param_i);
      
      % Random search
      d_theta = 2e0;
      for rand_init_idx = 1:25
        rng(rand_init_idx)
         theta0 = 2*d_theta*rand(2*Nc,1)-d_theta;
%          theta0_all(:,loop_i) = theta0;
%          theta_q: estimated parameters of class q (or class_i)
        [theta_q,Jval_q,exitflag_q] = fmincon(@(theta_q) moe_logistic_regression_cost_fn(theta_q,LR_cost_fn_param), ...
          theta0,[],[],[],[],lb,ub,[],options);
      
        theta_q_rand(:,rand_init_idx) = theta_q;
        Jval_q_rand(rand_init_idx) = Jval_q;
      end
      
      [min_Jval, min_Jval_idx] = min(Jval_q_rand);
      theta_q = theta_q_rand(:,min_Jval_idx);
      
      theta(:,class_i)     = theta_q;
%       Jval(loop_i,class_i) = min_Jval; 
      Jval(class_i,reg_param_i) = min_Jval;
      opt_rand_idx(class_i,reg_param_i) = min_Jval_idx;
%       exitflag(class_i)    = exitflag_q;
    end
    
    end
    
    if 1
      % Choose the best regularization parameter for each class
      [min_cost,best_idx] = min(Jval,[],2);
      best_reg_param = reg_param(best_idx);
      % The best regularization coefficient are (for q=0 to q=6 respectively):
      % 0   100   200   300   600   700     0
    end
    if 1
      % Debug: plot regularization parameter vs cost
      figure;
      plot(reg_param,Jval)
      xlabel('Regularization parameter value')
      ylabel('Cost')
      title('Learning curve: training')
      grid on
      legend('q=0','q=1','q=2','q=3','q=4','q=5','q=6','Location','best')
    end
%
  else
    % Use theta if you have it estimated already
    % theta =
  end
  %% Testing
  if testing_flag
    % Load testing data
    % ------------------
    % m_testing = size(eigenvalues_testing{1},1);
    
    if 0 && m_testing > 0
      % Use portion of the training dataset for testing
      % -----------------------------------------------
      % testing_input = [];
      % for q_i = 1 + [0:max_Ntargets]
      %   training_input_reshaped = reshape(training_input(q_i,:),[length(training_input)/(max_Ntargets+1), max_Ntargets+1]);
      %   testing_input(:,end+1:end + m_training) = training_input_reshaped(m_testing+1:end,:);
      % end
      
      % eigenvalues_testing = data_testing.eigenvalues_runs_snr; % Each cell is Nruns by Ntargets by Nb * Nc
      % LL_testing         = data_testing.loglikelihood_runs_snr;
      % LL_testing_opt     = LL_testing.opt;    % Each cell is Nruns by Ntargets by Ntargets
      % LL_testing_subopt  = LL_testing.subopt; % Each cell is Nruns by Ntargets by Ntargets
      
      % Geometric and arithmetic mean of the eigenvalues
      if strcmp(input_type,'geometric/arithmetic means')
        clear geom_mean arithm_mean
        for snr_i = 1:length(data_training.eigenvalues_runs_snr)
          for run_i = 1:m_testing;%size(data_testing.eigenvalues_runs_snr{1},1)
            for band_i = 1:size(data_training.eigenvalues_runs_snr{1},3)
              for q_i = 1:size(data_training.eigenvalues_runs_snr{1},2)
                eigval = squeeze(data_training.eigenvalues_runs_snr{snr_i}(run_i+m_training,q_i,band_i,:));
                for k_i = 0:size(data_training.eigenvalues_runs_snr{1},4)-1
                  geom_mean_testing{snr_i}(run_i,q_i,band_i,k_i+1)   = sum(log(eigval(k_i+1:end)));
                  arithm_mean_testing{snr_i}(run_i,q_i,band_i,k_i+1) = (Nc-k_i)*(log(sum((eigval(k_i+1:end))/(Nc-k_i))));
                  
                  %     log_lq(nb_idx,k+1) = Nsnap*(geom_mean- arithm_mean);
                end
              end
            end
          end
        end
      end
      
      Nexamples_testing_q = m_testing*length(eigenvalues_training);
      Nexamples = Nexamples_testing_q*(max_Ntargets+1);
      testing_input = [];
      for q_i = 1 + [0:max_Ntargets]
        % Concatinate the input data from all SNRs for this number of sources
        for snr_i = 1:length(eigenvalues_training)
          if strcmp(input_type,'eigenvalues')
            input_snr = eigenvalues_training{snr_i};
          elseif strcmp(input_type,'loglikelihoods')
            % Loglikelihoods were stored as -2L. So, here we divide by -2.
            input_snr = LL_training_opt{snr_i}*(-1/2);
            %       input_snr = exp(LL_training_subopt{snr_i}*(-1/2));
          elseif strcmp(input_type,'geometric/arithmetic means')
            input_snr1 = geom_mean_testing{snr_i}(1:m_testing,q_i,:,:);
            input_snr2 = arithm_mean_testing{snr_i}(1:m_testing,q_i,:,:);
            %       input_snr3 = eigenvalues_training{snr_i}(m_training+1:end,:,:,:);
            %       input_snr3 = input_snr1./input_snr2;
            %       input_snr = cat(4,input_snr1,input_snr2,input_snr3);
            input_snr = cat(4,input_snr1,input_snr2);
          else
            warning('Unsupported testing-input data type.')
            keyboard
          end
          %     testing_input(end+1:end+m_testing,:) = mean(real(squeeze(input_snr(1:m_testing,q_i,:,:))),2);
%           testing_input(end+1:end+m_testing,:) = real(squeeze(input_snr(1:m_testing,q_i,:,:)));
          testing_input(end+1:end+m_testing,:) = real(squeeze(input_snr));
        end
      end
      
      % Handle the case of bad examples (i.e. examples that have NaN or Inf entries)
      bad_examples_idxs = [];
      for bad_example_i = 1:size(testing_input,1)
        if any(isnan(testing_input(bad_example_i,:))) || any(isinf(testing_input(bad_example_i,:)))
          bad_examples_idxs(end+1) = bad_example_i;
        end
      end
      testing_input(bad_examples_idxs,:) = min(testing_input(:));% [];
      
      % Normalize the input data so that it is insensitive to scaling. Choosing
      % the right option is input value-dependent (trial and error type of scenario).
      if 0
        % Normalize wrt the maximum value of each example (works best for
        % eigenvalues)
        testing_input = testing_input./repmat(max(abs(testing_input),[],2),[1 size(testing_input,2)]);
      elseif 0
        % Subtract the mean and normalize wrt the maximum value of each example
        testing_input = (testing_input-repmat(mean(testing_input,2),[1 size(testing_input,2)]))./repmat(max(testing_input,[],2),[1 size(testing_input,2)]);
      elseif 1
        % Subtract the mean and normalize wrt the std value of each example
        testing_input = (testing_input-repmat(mean(testing_input,2),[1 size(testing_input,2)]))./repmat(std(testing_input,[],2),[1 size(testing_input,2)]);
      elseif 0
        % Subtract the mean and normalize wrt the std value of each example
        testing_input = (testing_input-repmat(mean(testing_input,1),[Nexamples 1]))./repmat(std(testing_input,[],1),[Nexamples 1]);
        %   testing_input = testing_input./repmat(max(abs(testing_input),[],2),[1 Nc]);
      end
      
      % Make testing_input an Nc by Nexamples matrix
      testing_input = testing_input.';
      % testing_input = [ones(1,size(testing_input,2)); testing_input];
      
    else
      % Use different dataset
      % ----------------------
      data_testing = load('/users/mohanad/IceSheetProject/MOE work/DCMandEigenvalues/1D_WB_LL_eigenval_5000runs.mat');
      eigenvalues_testing = data_testing.eigenvalues_runs_snr; % Each cell is Nruns by Ntargets by Nb * Nc
      LL_testing         = data_testing.loglikelihood_runs_snr;
      LL_testing_opt     = LL_testing.opt;    % Each cell is Nruns by Ntargets by Ntargets
      LL_testing_subopt  = LL_testing.subopt; % Each cell is Nruns by Ntargets by Ntargets
      
      m_testing = size(data_testing.eigenvalues_runs_snr{1},1);
      
      % Geometric and arithmetic mean of the eigenvalues
      if strcmp(input_type,'geometric/arithmetic means')
        clear geom_mean arithm_mean
        for snr_i = 1:length(data_testing.eigenvalues_runs_snr)
          for run_i = 1:m_testing;%size(data_testing.eigenvalues_runs_snr{1},1)
            for band_i = 1:size(data_testing.eigenvalues_runs_snr{1},3)
              for q_i = 1:size(data_testing.eigenvalues_runs_snr{1},2)
                eigval = squeeze(data_testing.eigenvalues_runs_snr{snr_i}(run_i,q_i,band_i,:));
                for k_i = 0:size(data_testing.eigenvalues_runs_snr{1},4)-1
                  geom_mean_testing{snr_i}(run_i,q_i,band_i,k_i+1)   = sum(log(eigval(k_i+1:end)));
                  arithm_mean_testing{snr_i}(run_i,q_i,band_i,k_i+1) = (Nc-k_i)*(log(sum((eigval(k_i+1:end))/(Nc-k_i))));
                  
                  %     log_lq(nb_idx,k+1) = Nsnap*(geom_mean- arithm_mean);
                end
              end
            end
          end
        end
      end
      
      Nexamples_testing_q = m_testing*length(eigenvalues_testing);
      Nexamples = Nexamples_testing_q*(max_Ntargets+1);
      testing_input = [];
      for q_i = 1 + [0:max_Ntargets]
        % Concatinate the input data from all SNRs for this number of sources
        for snr_i = 1:length(eigenvalues_testing)
          if strcmp(input_type,'eigenvalues')
            input_snr = eigenvalues_testing{snr_i};
          elseif strcmp(input_type,'loglikelihoods')
            % Loglikelihoods were stored as -2L. So, here we divide by -2.
            input_snr = LL_testing_opt{snr_i}*(-1/2);
            %       input_snr = exp(LL_testing_subopt{snr_i}*(-1/2));
          elseif strcmp(input_type,'geometric/arithmetic means')
            input_snr1 = geom_mean_testing{snr_i}(1:m_testing,q_i,:,:);
            input_snr2 = arithm_mean_testing{snr_i}(1:m_testing,q_i,:,:);
            %       input_snr3 = eigenvalues_training{snr_i}(m_training+1:end,:,:,:);
            %       input_snr = input_snr1./input_snr2;
            %       input_snr = cat(4,input_snr1,input_snr2,input_snr3);
            input_snr = cat(4,input_snr1,input_snr2);
            
          else
            warning('Unsupported testing-input data type.')
            keyboard
          end
          
          %     testing_input(end+1:end+m_testing,:) = mean(real(squeeze(input_snr(1:m_testing,q_i,:,:))),2);
%           testing_input(end+1:end+m_testing,:) = real(squeeze(input_snr(1:m_testing,q_i,:,:)));
          testing_input(end+1:end+m_testing,:) = real(squeeze(input_snr));
        end
      end
      
      % Handle the case of bad examples (i.e. examples that have NaN or Inf entries)
      bad_examples_idxs = [];
      for bad_example_i = 1:size(testing_input,1)
        if any(isnan(testing_input(bad_example_i,:))) || any(isinf(testing_input(bad_example_i,:)))
          bad_examples_idxs(end+1) = bad_example_i;
        end
      end
      testing_input(bad_examples_idxs,:) = min(testing_input(:));% [];
      
      % Normalize the input data so that it is insensitive to scaling. Choosing
      % the right option is input value-dependent (trial and error type of scenario).
      if 0
        % Normalize wrt the maximum value of each example (works best for
        % eigenvalues)
        testing_input = testing_input./repmat(max(abs(testing_input),[],2),[1 size(testing_input,2)]);
      elseif 0
        % Subtract the mean and normalize wrt the maximum value of each example
        testing_input = (testing_input-repmat(mean(testing_input,2),[1 size(testing_input,2)]))./repmat(max(testing_input,[],2),[1 size(testing_input,2)]);
      elseif 1
        % Subtract the mean and normalize wrt the std value of each example
        testing_input = (testing_input-repmat(mean(testing_input,2),[1 size(testing_input,2)]))./repmat(std(testing_input,[],2),[1 size(testing_input,2)]);
      elseif 0
        % Subtract the mean and normalize wrt the std value of each example
        testing_input = (testing_input-repmat(mean(testing_input,1),[Nexamples 1]))./repmat(std(testing_input,[],1),[Nexamples 1]);
        %   testing_input = testing_input./repmat(max(abs(testing_input),[],2),[1 Nc]);
      end
      
      % Make testing_input an Nc by Nexamples matrix
      testing_input = testing_input.';
      % testing_input = [ones(1,size(testing_input,2)); testing_input];
      
    end
    
    % Classify
    % ----------
    h = LR_cost_fn_param.model;
    L_total = Nexamples_testing_q;
    clear thr_q x_test theta_q percentage_correct
    for test_i = 1+[0:max_Ntargets]
%       if test_i == 6
%         keyboard;
%       end
      % x_test: test set corresponding to this class
      x_test = testing_input(:,(test_i-1)*Nexamples_testing_q+1:test_i*Nexamples_testing_q);
      
      % Now apply the max_Ntargets+1 models to this dataset
      for q_i = 1+[0:max_Ntargets]
        theta_q = theta(:,q_i);
        %    z_q = theta_q.'*x_test;
        thr_q(:,q_i) = h(theta_q,x_test);
      end
      [max_prob, class_i] = max(thr_q,[],2);
      pred_class = class_i;
      
      % percentage_correct: Nc by Nc matrix, where each column corresponds to
      % one class (or a specific number of possible targets) and each row corrsponds
      % to a specific dataset of one class. The entries represent the percentage
      % of the number of times the class is estimated correctly. Ideally, we
      % should get all zeros in each column, except for one value that is 100.
      % In other wors, percentage_correct is a diagonal matrix in the ideal case.
      for q_i = 1+[0:max_Ntargets]
        actual_class = q_i;
        u = find(pred_class == actual_class);
        percentage_correct(test_i,q_i) = length(u)/L_total *100;
      end
    end
    
    % format shortG
    percentage_correct = round(percentage_correct,2)
    percentage_correct_all(:,loop_i) = diag(percentage_correct);
    
    % Calculate the cost
    % ------------------
    for q_i=1:max_Ntargets+1
      % Calculate the cost using the testing dataset
      testing_label = zeros(size(testing_input,2),1);
      testing_label((q_i-1)*Nexamples_testing_q+1:q_i*Nexamples_testing_q) = q_i;
      
      LR_cost_fn_param.training_input  = testing_input;
      LR_cost_fn_param.training_output = testing_label;
      LR_cost_fn_param.model           = h;
      
      theta_q = theta(:,q_i);
%       J_testing(q_i,loop_i) = 1/(2*length(testing_label))*sum((h(theta_q,testing_input).' - testing_label).^2);
      J_testing(q_i,loop_i) = moe_logistic_regression_cost_fn(theta_q,LR_cost_fn_param);
      
      training_label = zeros(size(training_input,2),1);
      training_label((q_i-1)*Nexamples_training_q+1:q_i*Nexamples_training_q) = q_i;
      
      LR_cost_fn_param.training_input  = training_input;
      LR_cost_fn_param.training_output = training_label;
      LR_cost_fn_param.model           = h;
      
%       J_training(q_i,loop_i) = 1/(2*length(training_label))*sum((h(theta_q,training_input).' - training_label).^2);
      J_training(q_i,loop_i) = moe_logistic_regression_cost_fn(theta_q,LR_cost_fn_param);
      
    end
  end
  
  if 0
    % Debug: plot the logistic function
    figure(2000);clf;
    for q_i = 1+[0:max_Ntargets]
      theta_q = theta(:,q_i);
      % z_q should match the model (h)
      z_q = theta_q(1:size(testing_input,1)).'*testing_input.^(-1/2)+theta_q(size(testing_input,1)+1:2*size(testing_input,1)).'*testing_input.^(-1/3);
      subplot(2,ceil((max_Ntargets+1)/2),q_i)
      plot(z_q,h(theta_q,testing_input))
      xlim([min(z_q)  max(z_q)])
      xlabel('\theta^T x')
      ylabel('Cost')
      title(sprintf('pr(y=%1.0d|x;theta)',q_i))
      grid on
    end
  end
  
end

if 0
  % Save the results
  results_WB.percentage_correct = percentage_correct_all;
  results_WB.Jtraining=J_training;
  results_WB.Jtesting=J_testing;
  results_WB.m = perc/100*m;
  results_WB.note='Each row of Jtraining and Jtesting belongs to one class and each column belongs to the cost associated with a specific number of training examples, m';
end

 if 0
    % Plot the cost function
    figure(2001);clf
    figure(2002);clf
    for q_i=1:max_Ntargets+1
      figure(2001)
      subplot(2,4,q_i)
      plot(SampleSize*3,J_testing(q_i,:),'b','LineWidth',2)
%       semilogx(perc/100*m*3,J_testing(q_i,:),'b','LineWidth',2)
%       hold on
      
      xlim([min(SampleSize*3)  max(SampleSize*3)])
      grid on
      
      if q_i == 5 || q_i == 6 || q_i == 7
%         xlabel('m (*3000)')
        xlabel('m')
      end
      if q_i == 1 || q_i == 5
        ylabel('Cost')
      end
      title(sprintf('q=%1d',q_i-1))
 %     
      figure(2002)
      subplot(2,4,q_i)
      plot(SampleSize*3,J_training(q_i,:),'-r','LineWidth',2)
%       semilogx(perc/100*m*3,J_training(q_i,:),'-r','LineWidth',2)
      xlim([min(SampleSize*3)  max(SampleSize*3)])
      grid on
      
      if q_i == 5 || q_i == 6 || q_i == 7
%         xlabel('m (*3000)')
        xlabel('m')
      end
      if q_i == 1 || q_i == 5
        ylabel('Cost')
      end
      title(sprintf('q=%1d',q_i-1))
    end
    figure(2001)
    suptitle('Learning curves: testing case')
    figure(2002)
    suptitle('Learning curves: training case')
 end
  
if 0
  % Debug: Plot learning curves
  figure(103);clf
  for q_i = 1:7
    percentage_correct_q = percentage_correct_all(q_i,:);
    
    subplot(2,4,q_i)
    plot(perc/100*m/1000,percentage_correct_q)
    grid on
    xlim([min(perc/100*m/1000)  max(perc/100*m/1000)])
    if q_i == 5 || q_i == 6 || q_i == 7
      xlabel('Sample size per SNR (*1000)')
    end
    if q_i == 1 || q_i == 5
      ylabel('Percentage correct')
    end
    title(sprintf('q=%1d',q_i-1))
  end
  
  subplot(248)
  plot(percentage_correct_all.')
  xlim([min(perc/100*m/1000)  max(perc/100*m/1000)])
  grid on
  xlabel('Sample size per SNR (*1000)')
  title('All in one plot')
  legend('q=0','q=1','q=2','q=3','q=4','q=5','q=6','Location','best')
  
  suptitle('Learning curves for NB case: same testing dataset')
end

if 0
  % Debug: plot effect of initialization
  figure(104);clf
  for q_i = 1:7
    percentage_correct_q = percentage_correct_all(q_i,:);
    
    subplot(2,4,q_i)
    plot(percentage_correct_q)
    grid on
    xlim([1  loop_i])
    if q_i == 5 || q_i == 6 || q_i == 7
      xlabel('Initial point index')
    end
    if q_i == 1 || q_i == 5
      ylabel('Percentage correct')
    end
    title(sprintf('q=%1d',q_i-1))
  end
  
  subplot(248)
  plot(percentage_correct_all.')
  xlim([1  loop_i])
  grid on
  xlabel('Initial point index')
  title('All in one plot')
  legend('q=0','q=1','q=2','q=3','q=4','q=5','q=6','Location','best')
  
  suptitle('Effect of initial point for NB case(with reg.): same testing dataset')
end

%% Plot input data
if 0
  % Plot the loglikelihoods
  clear mean_LL
  LL_idxs = repmat([1:7]',[1 m]);
  figure(100);clf
  for q = 0:Nc-1
    subplot(2,4,q+1)
    LL = squeeze(LL_training_opt{2}(1:m_training,q+1,:));
    LL = LL./repmat(nanmax(abs(LL),[],2),[1 Nc]);
    plot(LL_idxs,LL.','b*');
    hold on
    plot(nanmean(LL,1),'+r')
    xlim([1 7])
    grid on
    if q == 4 || q == 5 || q == 6
      xlabel('q')
    end
    if q == 0 || q == 4
      ylabel('Log-likelihood')
    end
    mean_LL(:,q+1) = nanmean(LL,1);
  end
  
  figure(100);subplot(248);
  plot(mean_LL,'*')
  xlim([1 7])
  grid on
  xlabel('q')
  title('Mean Log-likelihood')
  legend('q=0','q=1','q=2','q=3','q=4','q=5','q=6','Location','best')
end


if 0
  % Plot the eigenvalues
  eigval_idx=repmat([1:7]',[1 m_training]);
  figure(1);clf
  figure(2);clf
  x_eigenvalue = [];
  
  for q = 0:Nc-1
    clear y
    sample_eigenvalues = data_training.eigenvalues_runs_snr{2}; % Nruns by Ntargets by Nb * Nc
    x_eigenvalue = squeeze(sample_eigenvalues(1:m_training,q+1,:,:))';
    %   y = squeeze(arithm_mean{2}(1:m,q+1,:,:))';
    y_norm=x_eigenvalue./repmat(max(abs(x_eigenvalue),[],1),[7 1]);
    
    %     figure(q+1);clf
    figure(1); subplot(2,4,q+1)
    plot(eigval_idx,y_norm,'*b')
    hold on
    plot(mean(y_norm,2),'+r')
    mean_eigval(:,q+1) = mean(y_norm,2);
    
    xlim([1 7])
    if q == 4 || q == 5 || q == 6
      xlabel('Eigenvalue index')
    end
    if q == 0 || q == 4
      ylabel('Norm. eigenvalue')
    end
    grid on
    title(sprintf('q = %1d',q))
    
    y_diff = flipud(diff(flipud(y_norm),1));
    mean_diff_eigval(:,q+1) = mean(y_diff,2);
    figure(2);subplot(2,4,q+1)
    plot(eigval_idx(1:end-1,:),y_diff,'*b')
    hold on
    plot(mean(y_diff,2),'+r')
    xlim([1 6])
    if q == 4 || q == 5 || q == 6
      xlabel('Eigenvalue index')
    end
    if q == 0 || q == 4
      ylabel('Relative eigenvalue')
    end
    grid on
    title(sprintf('q = %1d',q))
  end
  figure(1);subplot(248);
  plot(mean_eigval,'*')
  xlim([1 7])
  grid on
  xlabel('Eigenvalue index')
  title('Mean eigenvalues')
  legend('q=0','q=1','q=2','q=3','q=4','q=5','q=6')
  
  figure(2);subplot(248);
  plot(mean_diff_eigval,'*')
  xlim([1 6])
  grid on
  xlabel('Eigenvalue index')
  title('Mean difference-eigenvalues')
  legend('q=0','q=1','q=2','q=3','q=4','q=5','q=6')
  
  figure(1); suptitle('1D sim-Norm. eigenvalues: Narrowband, 1000 runs, 21 snapshots')
  figure(2); suptitle('1D sim-Relative eigenvalues: Narrowband, 1000 runs, 21 snapshots')
end


