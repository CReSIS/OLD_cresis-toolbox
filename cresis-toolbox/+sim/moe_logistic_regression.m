%% Data loading and setup
% data_training = load('/users/mohanad/IceSheetProject/MOE work/ML MOE/training  and testing data/1D_WB_10000runs.mat');
% data_training = load('/users/mohanad/IceSheetProject/MOE work/ML MOE/training  and testing data/1D_NB_5000runs.mat');
data_training = load('/users/mohanad/IceSheetProject/MOE work/ML MOE/training  and testing data/1D_WB_5000runs_edgeSNRs.mat');
% data_training = load('/users/mohanad/IceSheetProject/MOE work/ML MOE/training  and testing data/1D_NB_5000runs_edgeSNR.mat');
% data_training = load('/users/mohanad/IceSheetProject/MOE work/NT results/1D MOE_WB_2 DOAs_63snaps_10000runs/1D.mat');

% data_training.eigenvalues_runs_snr = data_training.eigenvalues_runs_snr_decimated;
eigenvalues_training = data_training.eigenvalues_runs_snr; % Each cell is Nruns by Ntargets by Nb * Nc

% This is used to train on all SNRs but test on one specific SNR.
single_snr_test_flag = 0;

if single_snr_test_flag
  reqd_snr_idx = 3;
  data_training = rmfield(data_training,'eigenvalues_runs_snr');
  data_training.eigenvalues_runs_snr{1} = eigenvalues_training{reqd_snr_idx};
  eigenvalues_training = data_training.eigenvalues_runs_snr;
end

% LL_training         = data_training.loglikelihood_runs_snr;
% LL_training_opt     = LL_training.opt;    % Each cell is Nruns by Ntargets by Ntargets
% LL_training_subopt  = LL_training.subopt; % Each cell is Nruns by Ntargets by Ntargets
% m: number of training examples. Don't use the same training dataset for
% testing. Split the training examples into 70 percent for training and 30
% percent for testing.
m = size(eigenvalues_training{1},1);
% SampleSize = round(logspace(log10(10),log10(m),5));
SampleSize = 3000;% m % this is the number of monte carlo runs. the actual number of examples is: m*7*n_snr
for loop_i = 1:length(SampleSize)
  clearvars -except single_snr_test_flag reqd_snr_idx SampleSize m loop_i percentage_correct_all training_input data_training eigenvalues_training J_testing J_training Jval
  h1 = sprintf('\nLoop # %d of %d...\n',loop_i,length(SampleSize));
  disp(h1)
  % rng(loop_i)
  
  m_training = SampleSize(loop_i);
  m_testing  = floor(m - m_training);
  
  Nc = 7;
  
  % max_Ntargets: maximum expected/required number of targets
  % (q=0:max_Ntargets)
  max_Ntargets = 6;
  
  % Nexamples_training_q: number of training examples per specific number of targets (Q)
    Nexamples_training_q = m_training*length(eigenvalues_training);
    Nexamples = Nexamples_training_q*(max_Ntargets+1);
    
   h = @(theta,x) (1./(1+exp(-(theta(1:Nc).'*x(1:Nc,:).^(1) + theta(Nc+1:2*Nc).'*x(Nc+1:end,:).^(1)))));
   LR_cost_fn_param.model = h;
   
  % Choos input type: 'loglikelihoods', 'geometric/arithmetic means'
  % (geometric and arithmetic means of eigenvalues), or 'eigenvalues'.
  % input_type = 'eigenvalues';
  % input_type = 'loglikelihoods';
  input_type = 'geometric/arithmetic means';
  
  % training_flag/testing_flag: boolean variable to turn on/off training or testing parts/modes
  training_flag  = 1;
  testing_flag   = 1;
  optimizer_flag = 0;
  
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
    
    if 0
      % Stor training data in a different format (for Victor)
      m_class = size(training_input,2)/(max_Ntargets+1);
      m_snr = m_class/length(data_training.eigenvalues_runs_snr);
      for class_i=1:7
        % Training data per class
        training_input_class = training_input(:,(class_i-1)*m_class+1:class_i*m_class);
        for snr_i=1:3
          training_input_snr{snr_i,class_i} = training_input_class(:,(snr_i-1)*m_snr+1:snr_i*m_snr).';
        end
      end
      
      
      
    end
  
  end
  
  %% Model and initialization
  % h: Model of the in-out relationship. It is the most important piece.
  % Watch the size of the parameter theta0 (initial value of theta).
  if 0
    % This model works best for NB eigenvalues
    h = @(theta,x) (1./(1+exp(-(theta(1:size(x,1)).'*x.^(-1/2) + theta(size(x,1)+1:2*size(x,1)).'*x.^(-1/3)))));
    theta0 = zeros(2*size(training_input,1),1);
    %   theta0 = rand(2*size(training_input,1),1)*(2*1e-2)-1e-2;
  elseif 0
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
%     reg_param =  [1e-5  1e-4  1e-3  1e-2  1e-1  0  1:50:500];
    reg_param = 1e-5*ones(1,max_Ntargets+1);
%     reg_param = [0 100 900 500 200 700 600];
%     reg_param = [0 1e-4 1 1e-3 1e-1 0 1];
    
    lb = -inf(size(training_input,1),1);
    ub = +inf(size(training_input,1),1);
%     lb = -inf(length(theta0),1);
%     ub = +inf(length(theta0),1);
    
    options =  optimoptions(@fmincon,'TolX',1e-6,'TolFun',1e-6,'MaxIter',10^7,'MaxFunEvals',10^7,'Display','off');%, ...
%       'PlotFcn',{@optimplotfval,@optimplotstepsize});%,'Algorithm','sqp');
    for reg_param_i = 1;%1:length(reg_param)
      h2 = sprintf('\n   lambda # %d of %d ...\n',reg_param_i,length(reg_param));
      disp(h2)
%     clear theta Jval exitflag
    training_label_all = [];
    for class_i = 1+[0:max_Ntargets]
      % Output training data (or class ID): Each time we set all outputs to 0
      % except for one class, which we set to class_i.
      training_label = zeros(size(training_input,2),1);
      training_label((class_i-1)*Nexamples_training_q+1:class_i*Nexamples_training_q) = 1;%class_i;
      
      % Cost function parameters (LR stands for Logistic Regression)
      LR_cost_fn_param.training_input  = training_input;
      LR_cost_fn_param.training_output = training_label;
      LR_cost_fn_param.model           = h;
%       LR_cost_fn_param.reg_param       = reg_param(class_i);
      LR_cost_fn_param.reg_param       = reg_param(reg_param_i);
      
      % Random search
      d_theta = 10;
      for rand_init_idx = 1:2
        rng(rand_init_idx)
        theta0 = 2*d_theta*rand(2*Nc,1)-d_theta;
%         theta0_all(:,loop_i) = theta0;
%         theta_q: estimated parameters of class q (or class_i)
        if 1
          % Mohanad's
          [theta_q,Jval_q,exitflag_q] = fmincon(@(theta_q) moe_logistic_regression_cost_fn(theta_q,LR_cost_fn_param), ...
            theta0,[],[],[],[],lb,ub,[],options);
        else
          % Victor's
          options = optimset('MaxIter', 1e3, 'Display', 'off');
          [theta_q, Jval_q] = fmincg(@(theta_q) moe_logistic_regression_cost_fn(theta_q, LR_cost_fn_param),theta0, options);
        end
        theta_q_rand(:,rand_init_idx) = theta_q;
        Jval_q_rand(rand_init_idx) = Jval_q;
      end
      
      [min_Jval, min_Jval_idx] = min(Jval_q_rand);
      theta_q = theta_q_rand(:,min_Jval_idx);
      
      theta(:,class_i)     = theta_q;
      Jval(loop_i,class_i) = min_Jval; 
%       Jval(class_i,reg_param_i) = min_Jval;
      opt_rand_idx(class_i,reg_param_i) = min_Jval_idx;
%       exitflag(class_i)    = exitflag_q;
    end
    
    end
    
    if 0
      % Choose the best regularization parameter for each class
      [min_cost,best_idx] = min(Jval,[],2);
      best_reg_param = reg_param(best_idx);
      % The best regularization coefficient are (for q=0 to q=6 respectively):
      % 0   100   200   300   600   700     0
    end
    if 0
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
    % This theta is from Victor -- fantastic results
    if 1
      % WB: 1D simulator (lambda=1e-5)
      theta = [ ...
        115.3932  -41.8957   22.4887    8.4901  -22.1140  -21.7548   -3.1442
   21.8544   88.5577 -117.6403   12.4482   85.6742   -7.8794  -22.6042
  -65.2461  -57.7029  146.3293 -152.8717  -48.1071  113.1228    9.3234
  -88.2990   -3.3980  -74.0628  156.1764  -88.2640  -52.3726  100.3692
 -103.4253   -9.1910   10.8360  -34.6267   79.2251 -152.2212   18.7016
  -46.6482   20.4726   21.7595   29.1974   -5.0215  150.9595 -196.2829
   27.9845  -11.9824   -7.3009  -14.2238  -11.0077  -14.6233   22.2914
  -94.8233   16.1435   -7.3837   -5.8193   -0.4774   -0.6628   -4.1323
    6.6343  -49.2732   41.8384   -6.8791  -32.3504    5.9261    5.0152
  -15.5534   23.0922  -71.9905   58.6820   13.0195  -45.3194   -2.1512
   44.9007   -0.2421   37.8710  -71.9921   37.2458    6.8307  -51.2647
   34.5614    1.9546   -3.7568   18.5648  -47.3845   72.2483  -35.4074
   87.5448  -11.6786  -14.7524  -21.3240    3.3978  -78.1073   98.3919
   41.8129    1.7448   -6.7097   -8.0586    2.8344   -8.4311   28.4972];
    elseif 0
      % WB: 2D simulator .. 11 snapshots(lambda=1e-5):
      % Without decimation
      theta = [ ...
        -167.5825  123.8552  252.5915   30.8863  -29.1894  -17.3484  111.7195
        149.4602 -160.2116 -130.2510  179.1966   60.5440    7.7267  -31.3331
        35.8196   26.1451 -266.4613 -120.4866  148.5158   55.6375  -46.1517
        8.3693    8.9919  109.7341 -203.7672  -87.2943  136.3176   23.7713
        2.0931   14.4602   46.0689   99.6758 -188.5086  -96.8944   76.3352
        9.7902   11.3471   46.1115   49.1115  104.4304 -186.9385 -266.4674
        -7.3053    2.5993    2.0793    5.9393   14.3033   65.0841   63.9010
        23.5777  -14.5448  -59.3954   -7.1159    4.4380    2.2324  -34.8270
        -3.7391    2.5444   16.5945  -52.1294   -7.0552    2.5055   -4.4683
        -6.1513    2.7365   18.9980   16.2458  -52.0041  -12.1967    4.8376
        -5.1806    2.0853    8.6026   19.6762    9.8124  -60.8024  -15.6005
        -8.3819   -2.7354   -3.3521    4.5793   23.3928    9.0526  -49.6594
        -13.0411   -9.4493  -32.9770  -17.3280   -5.2655   40.9630  114.4649
        -11.2726   -1.3680   -1.8880    1.9720   10.3360   61.1168   59.9337];

      % WB: 2D simulator .. 21 snapshots(lambda=0.01):
         % Without decimation
%          theta = [ ...
%           -167.1903   15.0358  163.7966  149.0409   92.1379   32.9896   77.4535
%           124.5746    8.0789  -53.5423   12.5203   30.3144   17.5658  -14.1088
%           32.0806   15.9001  -98.2231  -53.3720   34.1875   38.0647  -20.4454
%           8.9547    7.2154   24.7653  -97.8965  -50.0509   62.8524    2.8378
%           0.0724    7.4588   21.1291   27.5027 -111.7288  -60.2421   53.2927
%           -3.0316   11.4625   21.0483   21.4199    3.9848 -183.4078 -165.7021
%           -16.9002   -4.1467   20.9568   30.3485   43.5045   68.0626   23.4884
%           13.3513   31.8381  -34.0264  -26.5855  -20.8924   -6.6623  -16.5832
%           -10.9763  -41.7811   34.8269  -18.8794  -13.8978   -5.3108   -2.5348
%           -14.7861  -25.8492   -2.4145   36.5712  -26.7474  -18.6122    0.7053
%           -20.4518   -3.1208   14.8071    1.7301   34.4575  -54.6090   -1.6361
%           -23.0940    5.6205   12.2739   16.9992    1.7156   37.0959  -66.5195
%           -21.9626    2.3116    7.3765   14.1959   29.8283   39.6167  159.0462
%           -17.5653   -7.4557   20.9028   30.4019   43.5113   67.7705   22.4239];
      
         % With decimation (decimation factor = 2)
%         theta = [ ...
%           -235.3037   43.0651  216.2298  190.0582   78.3813   17.7600   82.8526
%           200.2188   16.7953  -74.8196   25.8918   56.3869   19.0569  -11.2783
%           42.2876   39.3491 -152.9607  -81.7143   69.1676   61.8253  -29.3375
%           0.5292   26.0391   48.6358 -139.5413  -78.0503  115.2109   18.6032
%           5.0719   35.7559   22.4831   54.3935 -154.6644  -89.6614   82.7067
%           -7.2194   23.2655   18.1564   22.7139   13.9458 -258.1179 -248.2455
%           -20.2993  -19.4456   24.1408   37.7607   48.2334   85.8171   36.5587
%           23.3105   58.5045  -44.2185  -36.1748  -17.5235   -3.3469  -20.2353
%           -15.9600 -111.6534   31.0386  -22.4826  -26.2548   -3.8799   -2.7354
%           -20.2822  -39.2974    4.4196   41.4231  -36.2745  -33.1677   -0.6380
%           -27.3657  -30.9680   13.9578   10.6027   34.7168  -64.4500   -8.2711
%           -24.3350  -23.9020   10.4226   18.9229    6.9449   31.8059  -86.0773
%           -22.1291   17.1821    6.4859   16.9398   32.7460   60.4692  184.1125
%           -16.6049   -8.8025   26.4013   37.7631   48.3011   89.6365   43.6042];
     
       
    elseif 0
      % NB: 1D simulator (lambda=1e-5)
      theta = [ ...
        21.4380  -34.1855  -18.8216  -88.0541   10.1419  -55.3386    6.1199
        -5.8818    4.9128  -18.4896   65.8542  -18.2141   28.3606   32.8638
        -11.3701   16.1711    2.8871  -70.1359   50.2207   53.7158   22.9762
        -13.0747   11.9325   48.0741   59.4007 -199.0351   25.6690   11.0571
        -11.6139    5.4855   12.5482  114.0562  153.0354 -129.8525   21.6085
        -4.7807    4.5220   11.7470   10.5063   95.5616  138.0173   40.0402
        14.7222   -8.1734  -14.5300  -28.7718  -40.8947  -14.1418  -88.1731
        -12.5413   18.4482   -1.2671    6.5262   -5.4486    5.7333  -51.0520
        2.4435  -33.6139   47.9696   -4.0999   -1.1316    1.3580   39.7649
        2.9760    7.4611 -119.1599   97.4575  -12.1052   -9.8669    1.6307
        7.0628   -1.2289   81.6667 -203.7587   82.5779  -17.2509  -32.0586
        9.3928  -10.9943  -14.9405   62.5042  -88.5795   59.1064  -42.2346
        7.7229   -6.5373  -21.9593  -63.4295  -19.5545 -107.2161  132.0497
        14.7213   -8.2035  -14.1102  -25.6917  -27.3239   -0.1498  -88.1735];
    end
  
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
      data_testing = load('/users/mohanad/IceSheetProject/MOE work/ML MOE/training  and testing data/1D_WB_5000runs.mat');
%      data_testing.eigenvalues_runs_snr = data_testing.eigenvalues_runs_snr_decimated;
      eigenvalues_testing = data_testing.eigenvalues_runs_snr; % Each cell is Nruns by Ntargets by Nb * Nc
      
      if single_snr_test_flag
        data_testing = rmfield(data_testing,'eigenvalues_runs_snr');
        data_testing.eigenvalues_runs_snr{1} = eigenvalues_testing{reqd_snr_idx};
        eigenvalues_testing = data_testing.eigenvalues_runs_snr;
      end
%       LL_testing         = data_testing.loglikelihood_runs_snr;
%       LL_testing_opt     = LL_testing.opt;    % Each cell is Nruns by Ntargets by Ntargets
%       LL_testing_subopt  = LL_testing.subopt; % Each cell is Nruns by Ntargets by Ntargets
      
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
    
    %% Classify
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
    if 0
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

if 0
  % Plot percentage corrct
  plot(diag(percentage_correct),'b-*','LineWidth',2)
  grid on
  grid minor
%   pbaspect([1 1/2 1])
  legend('1D: train=WB .. test=NB')
end
return

%% To make 2D dataset match 1D
data_training = load('/users/mohanad/IceSheetProject/MOE work/DCMandEigenvalues/2D_WB_1000runs.mat');

eigval = data_training.eigenvalues_all;
num_targets = data_training.actual_num_targets{1}{1}(1:end-5);
class_length = length(find(num_targets==1));

% Each run has class_length examples for each class (so 45 examples)
eigenvalues_snr = [];
label_y = [];
y = [];
for run_i = 1:1000
  eigval_run = eigval{run_i};
  for snr_i=1:3
    eigval_snr = eigval_run{snr_i}; % Nt-by-Nc
    for class_i = 1:7
      true_class = find(num_targets == class_i-1);
      eigval_class = eigval_snr(true_class,:); % class_length-by-Nc
      
      eigenvalues_snr{snr_i}((run_i-1)*class_length+1:run_i*class_length,class_i,1,:) = eigval_class;
    end
  end
end
    
eigenvalues_runs_snr=eigenvalues_snr;
data_training.eigenvalues_runs_snr=eigenvalues_runs_snr;

    
  
