function results = crosstrack(param)
% results = sim.crosstrack(param)
%
% Function for simulating test data and testing direction of arrival
% algorithms on the test data.
%
% param: structure describing the simulation and testing
%  .src: structure describing signal source characteristics
%    .f0: low frequency (Hz)
%    .f1: high frequency (Hz)
%    .lever_arm.fh: function handle to a lever arm function that returns
%      the locations of the sensor phase centers. Phase centers returned
%      should be in a 3 by Nc array. First row is "x" or along-track
%      and is generally not used, second row is "y" , and third row is "z".
%      Nc is the number of sensors/channels of data.
%    .lever_arm.fh_args: arguments to be passed to the lever arm function
%      handle
%    FIELDS REQUIRED BY sim.crossover_data.m
%  .monte: structure describing the monte carlo simulation setup
%    .runs: positive integer scalar representing number of simulation runs
%    .random_seed_offset: scalar double usually set to zero to allow
%                         repeatability
%
% results: structure containing results from simulation runs
%
% Author: Sean Holloway, John Paden, Theresa Stumpf
%
% See also: sim.crosstrack.m, sim.crosstrack_data.m, sim.crosstrack_example.m

%% Setup

% Load standard physical constants like c = speed of light
physical_constants;
array_proc_methods;
% fc: center frequency
param.src.fc = (param.src.f0 + param.src.f1)/2;
% fs: sampling frequency
param.src.fs = abs(param.src.f1 - param.src.f0);
% dt: Fast time sample spacing (sec)
dt = 1/param.src.fs;
% Nt: number of range bins (fast time samples)
Nt = floor((param.src.t1-param.src.t0)/dt);

% phase_center: get phase center positions of antenna array.
% Convert phase_center's from aircraft coordinate system into flight coordinate system
% y_pc: cross track with positive pointing to the left
% z_pc: elevation with positive pointing up

if isfield(param.src,'phase_center')
  param.src.y_pc   = -param.src.phase_center(2,:).';
  param.src.z_pc   = -param.src.phase_center(3,:).'; 
elseif isfield(param.src,'lever_arm') && isfield(param.src.lever_arm,'fh')
  [phase_center] = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
  param.src.y_pc   = -phase_center(2,:).';
  param.src.z_pc   = -phase_center(3,:).';
elseif isfield(param.src,'y_pc') && isfield(param.src,'z_pc')
  % Do nothing. Already defined 
else
  warning('Phase center information are not defined')
  keyboard
end

%% Monte Carlo Trials
%==========================================================================
if ~isfield(param,'SNR_db') || isempty(param.SNR_db)
  % SNR_db is used with MOE only.
  % If param.SNR_db is not used, then its value doesn't matter, just
  % its length. param.SNR_db is used in model order estimation only.
  param.SNR_db = -999;
end
start_time = tic;
for run_idx = 1:param.monte.runs

  %% Setup random number generator
  if isfield(param.monte,'rng_seed') && ~isempty(param.monte.rng_seed)
    rng(param.monte.rng_seed + run_idx)
  else
    rng_args = param.monte.random_seed_offset + run_idx;
    rng(rng_args);
  end
%   rng(param.monte.rng_seed)
  %% Create layered surface model
  fprintf('  Create simulated data\n');
  for test_idx = 1:length(param.SNR_db)
    fprintf('\nTest %2d of %2d ... Run %4d of %d (%.1f sec)\n',test_idx,length(param.SNR_db), run_idx, param.monte.runs, toc(start_time));
  
    surf_model.x = [];
    surf_model.y = [];
    surf_model.z = [];
    surf_model.rcs = [];
    surf_model.layer = [];
    if 0 && param.debug_level >= 3
      % Plot surface model
      figure(11); clf;
    end
  
    for surf_idx = 1:length(param.monte.target_param)
%       if isfield(param,'SNR_db') || ~isempty(param.SNR_db)
      if isfield(param.monte,'rng_seed') && ~isempty(param.monte.rng_seed)
        param.monte.target_param{surf_idx}.rng_seed = param.monte.rng_seed + run_idx + 1;
      else
        param.monte.target_param{surf_idx}.rng_seed = param.monte.random_seed_offset + run_idx+1;
      end
    
      if isfield(param,'optimal_test') || isfield(param,'suboptimal_test')
        % SNR_db is used with MOE only.
        source_power = param.SNR_db(test_idx) + param.src.noise_power(1);
        param.monte.target_param{surf_idx}.rcs_in.var = 10.^(source_power/10) ;% LINEAR
      end
      
      new_surf = param.monte.target_func(param.monte.target_param{surf_idx});
      
      if surf_idx ==1
        surf_model.x = [surf_model.x; new_surf.x];
      end
%       surf_model.x = [surf_model.x; new_surf.x];
      surf_model.y = [surf_model.y; new_surf.y];
      surf_model.z = [surf_model.z; new_surf.z];
%       if surf_idx == 1
%         surf_model.z = [surf_model.z; new_surf.z];
%       else
        % Do not allow negative layer thickness
%         new_surf.z(new_surf.z<0) = 0;
%         surf_model.z = [surf_model.z; surf_model.z + new_surf.z];
%       end
      surf_model.rcs = [surf_model.rcs; new_surf.rcs];
      surf_model.layer = [surf_model.layer; surf_idx*ones(size(new_surf.rcs(:,1)))];
      
      if 0 & param.debug_level >= 3
        % Plot surface model
        figure(11);
        plot_style = {'b-','r:'};
        plot(new_surf.y, new_surf.z(:,1), plot_style{1+mod(surf_idx-1,length(plot_style))});
        hold on;
      end
    end
    
    %% Create simulated data
%     if isfield(param.monte.target_param{surf_idx},'MP_weights')
%       surf_model.MP_weights = param.monte.target_param{surf_idx}.MP_weights;
%     end
    [sim_data,sim_time] = sim.crosstrack_data(param, surf_model);
    %   param.src.y_pc = param.src.y_pc;
    %   param.src.z_pc = param.src.z_pc;
    
    if 0
      imagesc(squeeze(lp(sim_data(:,1,:))))
    end
    
    % =======================================================================
    %% Process simulated data
    % =======================================================================
%     param.src.y_pc   = phase_center(2,:).';
%     param.src.z_pc   = phase_center(3,:).';
    %% Setup array processing input arguments
    sim_data = {permute(sim_data,[1 3 4 5 2])};
    param.array.tomo_en = true;
    for wf_adc_idx = 1:length(param.src.y_pc)
      param.array.fcs{1}{wf_adc_idx}.pos(2,:) = repmat(param.src.y_pc(wf_adc_idx), [1 size(sim_data{1},2)]);
      param.array.fcs{1}{wf_adc_idx}.pos(3,:) = repmat(param.src.z_pc(wf_adc_idx), [1 size(sim_data{1},2)]);
     
      if isfield(param.array,'surface') && ~isempty(param.array.surface)
        % For S-MLE DoA estimation
        param.array.fcs{1}{wf_adc_idx}.surface = param.array.surface;
      end
    end
    
    param.array_proc.wfs.fc = (param.src.f1 + param.src.f0)/2;
    param.array_proc.wfs.time = sim_time;
    
    if isfield(param,'optimal_test') || isfield(param,'suboptimal_test')
      % For model order estimation (MOE)
      sim_data_SNR{test_idx}= sim_data;
    end
    
    %% Determine the actual DoAs and their number per pixel
    if isfield(param,'optimal_test') || isfield(param,'suboptimal_test')
      % For model order estimation (MOE)
      actual_doa_param.rlines       = 1:size(surf_model.z,2);
      actual_doa_param.rbins        = 1:length(param.array_proc.wfs.time)-1;
      actual_doa_param.BW           = param.src.fs;
%       actual_doa_param.phase_center = phase_center;
      actual_doa_param.src.y_pc     = param.src.y_pc;
      actual_doa_param.src.z_pc     = param.src.z_pc;
      actual_doa_param.fs           = param.src.fs;
      actual_doa_param.surf_model   = surf_model;
      actual_doa_param.t0           = param.src.t0;
      actual_doa_param.fc           = param.src.fc;
      actual_doa_param.line_rng    = param.array.line_rng;
      actual_doa_param.dline        = param.array.dline;
      
      if isfield(param,'SS') && ~isempty(param.SS) && param.SS==0
        actual_doa_param.SS = param.SS;
      end
      
      if isfield(param,'dist_target_factor') && ~isempty(param.dist_target_factor)
        actual_doa_param.dist_target_factor = param.dist_target_factor;
      end
      
      if isfield(param,'optimal_test')
        actual_doa_param.optimal_test = param.optimal_test;
      end
      
      if isfield(param,'suboptimal_test')
        actual_doa_param.suboptimal_test = param.suboptimal_test;
      end
      
      [actual_doa, actual_doa_len] = sim.actual_number_of_targets(actual_doa_param);
      
      actual_doa_all_SNR{test_idx}   = actual_doa;
      sources_true_all_SNR{test_idx} = actual_doa_len;
      clear sim_data actual_doa_len  actual_doa
    else
      % All simulations that don't use MOE
      doa_param.rlines = abs(param.array.line_rng(1))+1:param.array.dline:length(surf_model.x)-abs(param.array.line_rng(end));
      doa_param.rbins = abs(param.array.bin_rng(1))+1:param.array.dbin:Nt-abs(param.array.bin_rng(end));
      doa_param.phase_center(2,:) = -param.src.y_pc;
      doa_param.phase_center(3,:) = -param.src.z_pc;
      %     doa_param.phase_center = -param.phase_center;
      doa_param.src.y_pc     = param.src.y_pc;
      doa_param.src.z_pc     = param.src.z_pc;
      doa_param.BW = (param.src.f1-param.src.f0);
      doa_param.fc = param.src.fc;
      doa_param.fs = param.src.fs;
      doa_param.t0 = param.src.t0;
      doa_param.surf_model.y = surf_model.y;
      doa_param.surf_model.z = surf_model.z;
      doa_param.plot_doa_lims = [-pi/2 pi/2];%param.plot_doa_lims;
      
      [actual_doa,Nsrc_true] = sim.actual_number_of_targets(doa_param);
      if 0
        % Remove DoAs outside 3dB beamwidth (or pre-specified limits in
        % general)
        for line_idx = 1:size(actual_doa,2)
          for bin_idx = 1:size(actual_doa,1)
            good_doa = actual_doa{bin_idx,line_idx} ...
              (actual_doa{bin_idx,line_idx}>doa_param.plot_doa_lims(1) & actual_doa{bin_idx,line_idx}<doa_param.plot_doa_lims(2));
            
            actual_doa{bin_idx,line_idx} = good_doa;
            good_doa_idx = numel(good_doa);
            Nsrc_true(bin_idx,line_idx) = good_doa_idx;
          end
        end
      end
    end
   
  end
  
  if isfield(param,'optimal_test') || isfield(param,'suboptimal_test')
    % For model order estimation (MOE).
    % Save actual number of targets
    actual_num_targets{run_idx} = sources_true_all_SNR;
    actual_doa_targets{run_idx} = actual_doa_all_SNR;
  end
   
    %% Loop through and run each array processing method
    for method_idx = 1:length(param.method.list)
      
      fprintf('  Array processing method %s/%d\n', array_proc_method_str(param.method.list(method_idx)), param.method.list(method_idx));
      
      %% Run array processing
      % Some notes about the outputs
      %  param.array.theta: radians, zero points toward -z_pc, increases toward positive y_pc
      param.array.method = param.method.list(method_idx);
      if isfield(param.method,'method_mode') && ~isempty(param.method.method_mode)
        param.array.method_mode = param.method.method_mode;
      end
      if isfield(param,'optimal_test') || isfield(param,'suboptimal_test')
        %% For model order estimation (MOE)
        if isfield(param,'testing') && ~isempty(param.testing) && param.testing == 0
          % This section is used in the training phase of model order estimation
          param.array.optimal_test    = param.optimal_test;
          param.array.suboptimal_test = param.suboptimal_test;
          
          if ~isfield(param,'opt_norm')
            param.array.opt_norm = 0;
          else
            param.array.opt_norm = param.opt_norm;
          end
          
          if ~isfield(param,'optimizer')
            param.array.optimizer = 0;
          else
            param.array.optimizer = param.optimizer;
          end
          
          if ~isfield(param,'testing')
            param.array.testing = 0;
          else
            param.array.testing = param.testing;
          end
          
          param.array.SNR_db = param.SNR_db;
          param.array.norm_allign_zero = param.norm_allign_zero;
          
          % Determine the log-likelihoods for all SNRs and runs
          loglikelihood_2D_param.array_param          = param.array;
          loglikelihood_2D_param.sim_data_SNR         = sim_data_SNR;
          loglikelihood_2D_param.y_pc                 = param.src.y_pc;
          loglikelihood_2D_param.z_pc                 = param.src.z_pc;
          loglikelihood_2D_param.norm_allign_zero     = param.array.norm_allign_zero;
          
          if ~isfield(param,'opt_norm_term') || isempty(param.opt_norm_term)
            param.opt_norm_term = zeros(1,max(param.array.Nsrc)+1);
          end
          if ~isfield(param,'norm_term_suboptimal') || isempty(param.norm_term_suboptimal)
            param.norm_term_suboptimal = zeros(1,max(param.array.Nsrc)+1);
          end
          if ~isfield(param,'norm_term_optimal') || isempty(param.norm_term_optimal)
            param.norm_term_optimal = zeros(1,max(param.array.Nsrc)+1);
          end
          
          loglikelihood_2D_param.opt_norm_term        = param.opt_norm_term;
          loglikelihood_2D_param.norm_term_suboptimal = param.norm_term_suboptimal;
          loglikelihood_2D_param.norm_term_optimal    = param.norm_term_optimal;
          
          % LL_subopt{run_idx} and LL_opt{run_idx} are cell arrarys, where
          % each cell contains the log-likelihoods for a single SNR. Also,
          % LL_subopt{run_idx}{SNR_idx} and LL_opt{run_idx}{SNR_idx} are
          % matrices of size (Nx*Nt)-by-M
          [LL_subopt{run_idx}, LL_opt{run_idx},~,eigenvalues_all{run_idx}] = loglikelihood_fun_2D(loglikelihood_2D_param);
          %         actual_num_targets{run_idx}          = sources_true_all_SNR;
          
          
          if 0 && exist('sources_true_all_SNR','var') && ~isempty(sources_true_all_SNR) && length(eigenvalues_all)==3
            % Debug: plot eigenvalues
            figure(999);clf
            clear rbin_idxs
            Nloops = length(param.param.array.Nsrc)+1;
            for q_idx = 1:Nloops %size(eivenvalues_all,3)
              bin_Idx = find(sources_true_all_SNR{1} == q_idx-1,1);
              rbin_idxs(q_idx) = bin_Idx;
              
              for i = 1:3
                 sample_eigenvalues(i,:) = eigenvalues_all{1}{i}(bin_Idx,:); 
              end
              eigenvalue_q = squeeze(sample_eigenvalues).';
              eigenvalue_q = 10*log10(eigenvalue_q);%./repmat(max(eigenvalue_q,[],1),[size(eigenvalue_q,1) 1]));
              
              if length(param.param.array.Nsrc)+1 > 3
                subplot(ceil(Nloops/2),2,q_idx)
              else
                 subplot(Nloops,1,q_idx)
              end
              
              plot(repmat([0:1:6].',[1 3]),eigenvalue_q,'-*')
              title(sprintf('q = %d',q_idx-1))
              
              if q_idx == Nloops % || q_idx == size(eivenvalues_all,3)
                xlabel('Eigenvalue index')
              end
              
              if q_idx == 1 %mod(q_idx,2) == 1
                ylabel('Eigenvalue (dB)')
              end
              grid on
              if q_idx == 1
                legend('10dB','20dB','30dB','Location','southwest')
              end
              xlim([0 6])
              if min(eigenvalue_q(:)) < max(eigenvalue_q(:))
                ylim([min(eigenvalue_q(:))  max(eigenvalue_q(:))])
              end
            end
            suptitle(sprintf('2D simulator: Bf = %2.2f percent',param.src.fs/param.src.fc *100))
          end
          
        else
          for test_idx = 1:size(param.SNR_db,2)
            if param.suboptimal_test==1
              param.array.NT = param.NT;   % suboptimal
            end
            
            if param.optimal_test==1
              param.array.NT_opt = param.NT_opt;   % Optimal
            end
            
            sim_data = sim_data_SNR{test_idx};
            param.array.Nsrc = max(param.array.Nsrc);
            param.array.testing = param.testing;
            param.array.suboptimal_test = param.suboptimal_test;
            param.array.optimal_test = param.optimal_test;
            if param.suboptimal_test
              param.array.penalty_NT = param.NT;
              %             else
              %               param.array.penalty_NT = zeros(max(param.array.Nsrc),1);
            end
            if param.optimal_test
              param.array.penalty_NT_opt = param.NT_opt;
              %             else
              %               param.array.penalty_NT_opt = zeros(max(param.array.Nsrc),1);
            end
            param.array.norm_allign_zero = param.norm_allign_zero;
            
            param.array.opt_norm_term        = param.opt_norm_term;
            param.array.norm_term_optimal    = param.norm_term_optimal;
            param.array.norm_term_suboptimal = param.norm_term_suboptimal;
            param.array.moe_methods      = param.moe_methods;
            [array_param_tmp{run_idx}{test_idx},tomo_tmp{run_idx}{test_idx}] = array_proc(param.array, sim_data);
          end
        end
        
      else
        %% All simulations that don't use MOE
        param.array_proc.bin0 = param.array_proc.wfs.time(1)/dt;
        [param,tomo] = array_proc(param, sim_data);
      end
      %% Debug Plots
      if exist('tomo','var')
        if param.debug_level >= 2
          range = param.array_proc.wfs.time(param.array_proc.bins)*c/2;
          
          if param.array.method < DOA_METHOD_THRESHOLD
            [theta,theta_idxs] = sort(param.array.theta);
            
            z = bsxfun(@times,-range,cosd(theta(theta_idxs)));
            y = bsxfun(@times,range,sind(theta(theta_idxs)));
            y_axis = (min(y(:)):min(diff(sort(y(1,:)))):max(y(:))).';
            dz_axis = min(diff(range));
            z_axis = min(z(:))-5*dz_axis:dz_axis:max(z(:))+5*dz_axis;
            
            slice = 1;
            tomo.img_rect = griddata(y,z,double(tomo.tomo.img(:,:,slice)),y_axis,z_axis);
            
            
            slice = 11;
            figure(1); clf;
            subplot(2,1,1);
            plot_style = {'b-','r:'};
            for layer = unique(surf_model.layer(:)).'
              if layer == 1
                layer_z = surf_model.z(surf_model.layer == layer,slice);
              else
                layer_z = surf_model.z(surf_model.layer == layer,slice);
              end
              plot(surf_model.y(surf_model.layer == layer,1), ...
                layer_z-param.monte.target_param{1}.z.mean, ...
                plot_style{1+mod(layer-1,length(plot_style))});
              hold on;
            end
            hold off;
            ylabel('Elevation (m)');
            h_axes = gca;
            grid on;
            
            subplot(2,1,2);
            imagesc(y_axis,z_axis-param.monte.target_param{1}.z.mean,lp(tomo.img_rect));
            set(gca,'YDir','normal');
            xlabel('Cross-track (m)');
            ylabel('Elevation (m)');
            h_axes(2) = gca;
            
            linkaxes(h_axes,'xy');
            ylim(h_axes(1),ylim(h_axes(1)) + dz_axis*[-5 5]);
            
          else
            z = bsxfun(@times,-range,cos(tomo.tomo.theta(:,:,1)));
            y = bsxfun(@times,range,sin(tomo.tomo.theta(:,:,1)));
            
            if 0
              figure(1); clf;
              subplot(2,1,1);
              plot_style = {'b-','r:'};
              for layer = unique(surf_model.layer).'
                plot(surf_model.y(surf_model.layer == layer), ...
                  surf_model.z(surf_model.layer == layer,1)-param.monte.target_param{1}.z.mean, ...
                  plot_style{1+mod(layer-1,length(plot_style))});
                hold on;
              end
              hold off;
              ylabel('Elevation (m)');
              h_axes = gca;
              grid on;
              
              subplot(2,1,2);
              plot(y,z-param.monte.target_param{1}.z.mean,'.');
              set(gca,'YDir','normal');
              xlabel('Cross-track (m)');
              ylabel('Elevation (m)');
              h_axes(2) = gca;
              grid on;
              linkaxes(h_axes,'xy');
            end
            
            if 0
              figure(2); clf;
              imagesc(surf_model.x, surf_model.y, surf_model.z-param.monte.target_param{1}.z.mean);
              title('Ground truth');
              h_axes = gca;
              xlabel('Along-track (m)');
              ylabel('Cross-track (m)');
              cc = caxis;
              h_cb = colorbar;
              set(get(h_cb,'YLabel'),'String','WGS-84 elevation (m)');
            end
            
            z = bsxfun(@times,-range,cos(tomo.tomo.theta));
            y = bsxfun(@times,range,sin(tomo.tomo.theta));
            x = repmat(permute(surf_model.x(param.array_proc.lines),[1 3 2]),[size(y,1) size(y,2) 1]);
            
            z(isnan(z)) = -9999;
            y(isnan(y)) = -9999;
            x(isnan(x)) = -9999;
            
            good_mask = zeros(size(tomo.tomo.theta));
            good_mask = good_mask | db(tomo.tomo.img) > 10;
            good_mask = good_mask | repmat(permute(tomo.tomo.cost,[1 3 2]),[1 size(good_mask,2) 1]) < -25;
            
            % Choose layer to isolate in z_grid
            layer = 1;
            min_z = min(min(surf_model.z(surf_model.layer == layer,:)));
            max_z = max(max(surf_model.z(surf_model.layer == layer,:)));
            good_mask = good_mask & z >= min_z & z <= max_z;
            
            z_grid = griddata(double(x(good_mask)),double(y(good_mask)),double(z(good_mask)), ...
              surf_model.x,surf_model.y);
            if 0
              figure(3); clf;
              
              imagesc(surf_model.x, surf_model.y, z_grid-param.monte.target_param{1}.z.mean);
              h_axes(2) = gca;
              title('Array processing with basic surface extraction');
              xlabel('Along-track (m)');
              ylabel('Cross-track (m)');
              caxis(cc);
              h_cb = colorbar;
              set(get(h_cb,'YLabel'),'String','WGS-84 elevation (m)');
              
              linkaxes(h_axes,'xy');
            end
          end
        end
      end
      
      if ~isfield(param,'optimal_test') || ~isfield(param,'suboptimal_test')
        % Non MOE simulations
        if param.debug_level >= 3
          % End early for debug testing of outputs
          if param.array.method ~= PF_METHOD
            results.param_crosstrack = param;
            results.tomo = tomo;
            results.sim_data = sim_data;
            results.surf_model = surf_model;
            results.actual_doa = actual_doa;
            results.Nsrc_true = Nsrc_true;
            
            % required for slice plots in crosstrack_example.m
            if exist('z_grid','var')
              results.z_grid = z_grid ;
            end
          else
            % Particle filter simulation
            results.sim_data = sim_data;
            results.actual_doa = actual_doa;
          end
          return
        end
      end
      %% Determine number of sources errors and source DOA errors
      
      % To do...
      
    end
  
end

if isfield(param,'optimal_test') || isfield(param,'suboptimal_test')
  % For model order estimation (MOE)
  if isfield(param.array,'testing') && ~isempty(param.array.testing) && param.array.testing == 0
    % This section us used in the training phase of model order estimation
    results.LL_subopt          = LL_subopt;
    results.LL_opt             = LL_opt;
    results.eigenvalues_all    = eigenvalues_all;
    results.actual_num_targets = actual_num_targets;
  else
    if param.debug_level >= 3
      % End early for debug testing of outputs
      results.param_crosstrack = param;
      results.tomo        = tomo_tmp;
      results.surf_model  = surf_model;
      results.actual_num_targets = actual_num_targets;
    end
  end
end

if isempty(param.method.list) %~exist('results','var')
  % This is for cases where you don't need DOA estimation
  results.param_crosstrack = param;
  results.sim_data = sim_data;
  results.surf_model = surf_model;
  results.actual_doa = actual_doa;
  results.Nsrc_true = Nsrc_true;
end
