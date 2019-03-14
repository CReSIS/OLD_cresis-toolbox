% Goal: given an array geometry, can we come up with a set of radar
% parameters such that the RMSE of the DoA is minimum? So, here we seek a
% radar operation mode that minimizes the RMSE of the DoA.
%
% The main script was copied from sim.crosstrack_example.m script.
%
physical_constants;

% param = read_param_xls(ct_filename_param('rds_param_2016_Greenland_Polar6.xls'),'20160413_01');
%
% param.cmd.frms = [];

fprintf('=======================================================================\n');
fprintf('  Radar depth sounder simulation: 2016_Greenland_Polar6 \n');
fprintf('  Linear array in y-dimension\n');
fprintf('=======================================================================\n');

if 0
  %   rline_rng = round([-1*linspace(3,250,51).' , linspace(3,250,51).']);   % Ntests by 2 (start and end of range-lins range)
  %   rline_rng = rline_rng(1:2:end,:);
%   rline_rng = [-[5  10   20  30   40  50 100  150 200 250 350].' , [5  10   20  30   40  50 100  150 200 250 350].'];
  rline_rng = [-ceil(((logspace(log10(10),log10(1000),21)).')./2), ceil(((logspace(log10(10),log10(1000),21)).')./2)];
  
  Ntests = size(rline_rng,1);
  snr_dB = 20;%10^(20/10);
  test_type = 'Number of snapshots';
elseif 0
  %   Bf = [0.1:0.1:0.6 0.7 0.8:0.1:2].'; %0.7
  Bf = [0.001 0.01 0.1 0.2 0.3 0.4 0.5 1 1.5 2].';
  Ntests = length(Bf);
  test_type = 'Fractional BW';
elseif 1
  snr_dB = linspace(5,40,15);
  snr = 10.^(snr_dB./10); % This is basically the rcs.var
  Ntests = length(snr);
  test_type = 'SNR (dB)';
end

% Ntests = 1;
Ntrials = 100;

% DOA method parameters: MLE = 7, WBMLE = 9
methods_list    = [2 7];

test_param = [];
rmse_tests = [];
clear rmse_tests_mean
for method_idx = 1:length(methods_list)
  param = [];
  param.method.list = methods_list (method_idx);
  
  for test_idx = 1:Ntests
    %% Setup simulation parameters
    % Debug level of 3 causes this function to stop early and output
    % simulated data, simulation parameters, and array processing results
    % in a "results" structure.
    param.debug_level = 3;
    
    %% Source parameters
    % param.fc     = 195e6;
    param.BW     = 30e6;
    if strcmp(test_type,'Fractional BW')
      param.fc     = param.BW/Bf(test_idx);
    else
      param.fc     = 195e9;
    end
    param.fs     = param.BW;
    param.src.f0 = param.fc-param.BW/2;
    param.src.f1 = param.fc+param.BW/2;
    lambda = c/param.fc;
    
    % param.src.f0      = 150e6;
    % param.src.f1      = 520e6;
    % param.fc          = (param.src.f1+param.src.f0)/2;
    % param.BW          = abs(param.src.f1-param.src.f0);
    % param.fs          = param.BW;
    
    param.pulse_width = 3e-6;
    % flight_height is the average/approximate flight height.
    flight_height     = 1500;%param.pulse_width*c;
    % guard_bellow matters more than guard_above because as this parameter
    % increases you are zooming out the slice view while keeping the resolution
    % the same (like having a better look to the surface).
    guard_above       = 500;%flight_height/10; % in meters
    guard_bellow      = 1000;%flight_height/2; % in meters
    % t0/t1: time gate to calculate the number of range bins=(t1-t0)*fs
    param.src.t0      = 2*(flight_height-guard_above)/c;
    param.src.t1      = 2*(flight_height+guard_bellow)/c;
    
    t_ref = 2*flight_height/c;
    Nt_above_surface = floor((t_ref-param.src.t0)*param.fs);% The surface starts at this range-bin
    Nt_bellow_surface = floor((param.src.t1-t_ref)*param.fs);
    Nt = floor((param.src.t1-param.src.t0)*param.fs);
    param.actual_Nt = Nt;
    range_res = c/(2*param.fs);
    
    Nt_st = 1;
    Nx_st = 1;
    
    param.src.ft_func = @(t) tukeywin_cont(t * param.BW);
    param.src.ft_wind = @(N) hanning(N);
%     param.src.ft_wind = @(N) boxcar(N);
    
    % Find the phase center for each channel
    % param.src.lever_arm.fh = @sim.lever_arm_example;
    Nc = 3;
    % Nsep: number of lambda/4 steps between sensors
    Nsep = 1;
%     tx_weights         = ones(Nc,1);
    % tx_weights         = hanning(Nc);
    % tx_weights         = blackman(Nc);
    % tx_weights         = [hanning(4).',0 0 0 0].'; % Steering left
    % tx_weights         = [0 0 0 0 hanning(4).'].'; % Steering right
    param.src.lever_arm.fh            = @sim.lever_arm_example;
    param.src.lever_arm.fh_args       = {[], 1, 1:Nc, [0; Nsep*c/param.fc/2; 0]};
    
    if 0
      lever_arm_fh       = @lever_arm;
      param.season_name  = '2016_Greenland_Polar6';
      param.radar_name   = 'rds';
      param.gps_source   = 'AWI-'; % 'ATM-'
      
      for rxchannel_idx = 1:Nc
        % For left sub-array (chan#1 to 8)
        %         phase_center(:,rxchannel_idx) = lever_arm_fh(param,tx_weights,rxchannel_idx);
        % For center sub-array (chan#9 to 16)
        param.phase_center(:,rxchannel_idx) = lever_arm_fh(param,tx_weights,rxchannel_idx+8);
        % For right sub-array (chan#17 to 24)
        %         phase_center(:,rxchannel_idx) = lever_arm_fh(param,tx_weights,rxchannel_idx+16);
      end
    end
    % d_y = lambda/2;
    % param.phase_center =  zeros(3,Nc);
    % param.phase_center(2,:) = d_y/2*(-floor((Nc-1)/2):floor(Nc/2));
    
    phase_center = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
    param.phase_center = -phase_center;
    y_pc = -phase_center(2,:).'; %-y is the direction of increasing y (negative y is associated with + theta).
    z_pc = -phase_center(3,:).';
    
    % y_pc = phase_center(2,:).';
    % z_pc = phase_center(3,:).';
    
    % Arguments for a linear array in y-dimension
    % param.src.lever_arm.fh_args: Arguments for a linear array in y-dimension
    % param.src.lever_arm.fh_args       = {[], 1, 1:Nc, [0; Nsep*c/param.fc/2; 0]};
    param.src.noise_power = zeros(1,Nc); %in dB
    % param.src.noise_power =  10*log10(BoltzmannConst*290*param.BW) + 2*ones(1,Nc);
    
    % DOA method parameters: MLE = 7, WBMLE = 9
    % param.method.list    = [2];
    
    if param.method.list == 2
      method.name = 'MUSIC';
    elseif param.method.list == 7
      method.name = 'MLE';
    elseif param.method.list == 9
      method.name = 'WBMLE';
    end
    
    %% Place targets on the surface at specific angles
%     ang = [-1.2154 -1.0654 -0.9484 -0.8481 -0.7580 -0.6751 -0.5974 -0.5236 -0.4528 -0.3844 -0.3178 -0.2527 -0.1886 -0.1253 -0.0625...         0
%           0.0625 0.1253 0.1886 0.2527 0.3178 0.3844 0.4528 0.5236 0.5974 0.6751 0.7580 0.8481 0.9484 1.0654 1.2154]';
%     ang = ang(ang>-31*pi/180 & ang<31*pi/180);
    
    ang = [-[20 17 15 13 11 9 7 5 3], [0 3 5 7 9 11 13 15 17 20]]'*pi/180;
%     ang = [-10 10]'*pi/180;
    
    %% Beam steering and DoA limits (All DoAs should be inside the 3dB beamwidth)
    k = 2*pi/(lambda/2);
    tx_weights         = ones(Nc,1);
    % tx_weights         = hanning(Nc);
    % Beam steering angle
    steer_ang = 0*pi/180;
    
    steering_delay = (y_pc*sin(steer_ang) - z_pc*cos(steer_ang))/c;
    % steering_delay = 2*(param.phase_center(2,:)*sin(steer_ang))/c;
    w = tx_weights.*exp(-1i*4*pi*param.fc*steering_delay);
    steering_mtx_fh = @(DoA) ((1/sqrt(length(y_pc)))*exp(1i*(y_pc*k*sin(DoA) - z_pc*k*cos(DoA))));
    
    % Generate radiation pattern (RP)
    theta_RP = linspace(-90,90,2048)*pi/180;
    %   RP = abs(sum(diag(w)*steering_mtx(theta_RP),1)).^2;
    RP = abs(w.'*steering_mtx_fh(theta_RP)).^2;
    RP = RP./max(RP);
    RP_dB = 10*log10(RP);
    
    % Focus only on the 3dB beam
    threshold_dB = -3.5;
    idxs_3dB = find(RP_dB>=threshold_dB);
    RP_3dB = RP_dB(idxs_3dB);
    theta_3dB = theta_RP(idxs_3dB);
    
    min_doa_lim = min(ang); %min(theta_3dB);
    max_doa_lim = max(ang); %max(theta_3dB);
    param.plot_doa_lims = [min_doa_lim max_doa_lim];
    
    % Restrict actual DoAs to be inside the beampattern
    % ang = ang(ang>min_doa_lim & ang<max_doa_lim);
    if isempty(ang)
      warning('Actual target locations are not specified')
      keyboard
    end
    
    % Complex transmit weights (or transmit LUT)
    if 0
      % Use this LUT for normal simulations
      comp_tx_weight = w.'*steering_mtx_fh(ang.');
    else
      % Use this LUT for comapring DOA estimation methods
      comp_tx_weight = ones(size(ang));
    end
    param.src.tx_weights = comp_tx_weight;
    
    %% Simulation Runs Setup
    
    % Cross track monte carlo setup
    param.monte.target_func = @surface_gen;
    param.monte.runs = 1;
    param.monte.random_seed_offset = 0;
    est_error = [];
    % Ntrials = 50;
    tic
    for trial_idx = 1:Ntrials
      sprintf('\n Method %d of %d (%s) ... Test %u of %u... Trial %u of %u\n',method_idx,length(methods_list),method.name,test_idx, Ntests,trial_idx,Ntrials)
      % SEE sim.surfgen FOR DEFINITION OF THESE PARAMETERS
      param.monte.rng_seed = trial_idx+test_idx;
      % Target surface parameters
      surf_param = [];
      % z.mean: mean flight height (the radar is the reference and z is pointing
      % downward)
      surf_param.z.mean = -flight_height;
      % z.rms_height: how much variability around the mean. For flat surfaces,
      % set it to 0.
      surf_param.z.rms_height = 0;
      % z.corr_length_x: this number is calculated as the azimuth resolution
      % times the number of snapshots. Shoud be chosen such that the ratio of
      % z.corr_length_x to z.mean (the flight height)<<1 (e.g. 0.16).
      surf_param.z.corr_length_x = 400;%0.05*abs(surf_param.z.mean);
      % corr_length_y: similar to corr_length_x
      surf_param.z.corr_length_y = 400;%0.05*abs(surf_param.z.mean);
      % rcs_in.mean: rcs is modeled as a Gaussian random variable with mean
      % rcs_in.mean and variance rcs_in.var.
      surf_param.rcs_in.mean = 0;
      % rcs_in.var: rcs variance. Use it to control SNR
      if strcmp(test_type , 'SNR (dB)')
        surf_param.rcs_in.var = 1;
%         surf_param.rcs_in.var = 10^((snr_dB(test_idx)-10*log10(Nc))/10);
        param.snr_db = snr_dB(test_idx);%-10*log10(Nc);
      elseif exist('snr_dB','var')
        surf_param.rcs_in.var = 1;%snr_dB-10*log10(Nc);
        param.snr_db = snr_dB;%-10*log10(Nc);
      else
        surf_param.rcs_in.var = 1;%1e2;
        param.snr_db = snr_dB;%-10*log10(Nc);
      end
      % dy: this is the cross-track resolution (spatial distance between snapshots).
      surf_param.dy = 10;
      % y_range: the total range over which the cross-track extends.
      surf_param.y_range = [-2500 2500];
      % dx: this is the azimuth resolution (spatial distance between
      % snapshots). John defines it as: sample spacing in x-dimension
      surf_param.dx = 10;
      % x_range: the total range over which the along-track extends.
      surf_param.x_range = [-2500 2500];
      % x: the range over which the imaged swath extends along the flight path (<= x_range).
      %   surf_param.x = [-1500:-1500+10]; %[-1500:10:1500]; %
      % y: the range over which the imaged swath extends in the cross-track
      % dimension (<= y_range). Note here that the spacing here represents how
      % far apart are your targets from each other. [-500:10:500] means that the
      % targets are separated by 10m. To control the location of the targets, do
      % something like y=[first target, second target, ..., last target]'. You
      % can also add multiple layers by doing something like y{1} and y{2} and
      % also z_mean = [z_mean of layer 1, z_mean of layer 2]. But x doesn't
      % change (i.e. single x). Use ground range formula to place targets if you
      % want.
      % surf_param.y = [-1500:10:1500].';
      
      surf_param.y = flight_height*tan(ang);
      
      if 0
        % Draw range-bins and surface targets
        figure(100);clf
        
        time = param.src.t0 + (1/param.fs)*(0:Nt-1).';
        R = time*c/2;
        %   R = (param.src.t0:1/param.fs:param.src.t1)*c/2;
        y = flight_height*tan(ang);
        
        th = pi:pi/100:2*pi;
        
        hold on
        for i = 1:length(R)
          x_c = y(round(length(y)/2))+R(i).*cos(th);
          x_s = 0+R(i).*sin(th);
          figure(100);plot(x_c,x_s)
        end
        scatter(y,surf_param.z.mean*ones(length(y),1),'b.','LineWidth',40)
        y_min = -sqrt(max(R)^2-flight_height^2);
        y_max = sqrt(max(R)^2-flight_height^2);
        y_extra = y_max/4;
        
        xlim([y_min-y_extra, y_max+y_extra])
        xlabel('y-axis (cross-track)')
        ylabel('z-axis')
        title('Sample slice showing the traget locations within the range-bins')
        grid on
      end
      
      if 0
        % Change of DoA from range-bin to the next. Needed for PF process model
        time = param.src.t0 + (1/param.fs)*(0:Nt-1).';
        R = time*c/2;
        rng_res = c/(2*param.BW);
        h = flight_height;
        
        good_R_idx = find(R>=h);
        theta_1 = acos(h./R(good_R_idx(1:end-1)));
        theta_2 = acos(h./(R(good_R_idx(2:end))+rng_res));
        
        delta_theta = (theta_2-theta_1)*180/pi;
      end
      
      % param.snr_db = 10*log10(surf_param.rcs_in.var^2)- 10*log10(param.src.noise_power);
      
      %% Array Processing parameters
      array_param = [];
      if param.method.list == 2
        param.method.method_mode = 'estimator';
      elseif param.method.list == 7
        param.method.method_mode = 'mle';
      end
      
      if 0
        % This is the standard way of calculating SVs
        
        % k: wavenumber
        k = 2*pi/(lambda/2);
        
        % My: over sampling factor
        My = 4;
        
        % dy: phase center spacing
        dy = Nsep*lambda/4;
        
        % dky and ky: y-component of wavenumber (spacing and axis)
        dky = 2*pi / (Nc*dy) / My;
        
        % dy = mean(diff(phase_center(2,:))); % The mean spacing
        % array_len = sum(diff(phase_center(2,:)));
        % dky = 2*pi / array_len / My;
        
        ky = dky * ifftshift(-floor(My*Nc/2) : floor((My*Nc-1)/2));
        
        % theta: theta values associated with ky axis
        theta = fftshift(asin(ky./k));
        %   theta = [theta,pi/2];
        if 0
          % Exclude angles that are outside the field of view (FOV), which is
          % determined based on the range vector and the flight height.
          R = (param.src.t0:1/param.fs:param.src.t1)*c/2;
          max_ang_FOV = max(atan(sqrt(R.^2-flight_height^2)/flight_height));
          theta_min_idx = find(theta<=(-max_ang_FOV), 1, 'last');
          theta_max_idx = find(theta>=max_ang_FOV, 1, 'first');
          theta = theta(theta_min_idx:theta_max_idx);
          % FOV = abs(theta(end)-theta(1))*180/pi;
          
          % sprintf('\nFoV: %f to %f\n\n',min(theta)*180/pi,max(theta)*180/pi)
        end
        
        array_param.Nsv = {'theta', asin(ky/k)};
      elseif 0
        % SVs can be calculated using the standard 64 DOA that we used to
        % process real data (the rds_Greenland_2014 data set)
        k = 2*pi/(lambda/2);
        load('/cresis/snfs1/dataproducts/ct_data/ct_tmp/sv_calibration/rds/2014_Greenland_P3/theta_cal.mat');
        array_param.Nsv = {'theta', theta};
      else
        k = 2*pi/(lambda/2);
        theta = linspace(-30,30,128)*pi/180;
        array_param.Nsv = {'theta', theta};
      end
        
      array_param.sv_fh = @array_proc_sv;
      
      % bin_rng/rline_rng: the range of range-bins/lines around a given
      % range-bin/line where the snapshots are taken when constructing the DCM.
      % Remember that length(surf_param.x)>=length(array_param.rline_rng).
      array_param.bin_rng   = 0;%-1:1;
      if  strcmp(test_type,'Number of snapshots')
        array_param.rline_rng = round([rline_rng(test_idx,1) : rline_rng(test_idx,2)]);
      else
        array_param.rline_rng = -10:10;
      end
      Nsnaps = length(array_param.bin_rng )*length(array_param.rline_rng);
      
      % dbin/dline: number of range-bins/lines to be skipped.
      array_param.dbin  = 1;
      array_param.dline = 1;
      
      N_skipped_rlines = 1;
      N_reqd_rlines    = 1;
      surf_param.x = [-750:N_skipped_rlines:-750+N_skipped_rlines*(N_reqd_rlines+2*max(array_param.rline_rng))-1];
      param.monte.target_param{1} = surf_param;
      
      Nx = length(surf_param.x);
      
      rlines = Nx_st+max(array_param.rline_rng):array_param.dline:Nx-max(array_param.rline_rng);
      rbins = Nt_st+max(array_param.bin_rng):array_param.dbin:Nt-max(array_param.bin_rng);
      
      % Maximim number of snapshots depends on the correlation length and the
      % azimuth resolution, and can be written approximately as:
      % array_param.rline_rng =...
      %   -floor(surf_param.z.corr_length_x/surf_param.dx/2):floor(surf_param.z.corr_length_x/surf_param.dx/2);
      
      % TO DO: ADD MODEL ORDER ESTIMATION.
      %====================================
      % Nsig: number of impinging signals.
      Nsig = 2;
      array_param.Nsig = Nsig;
      
      % init: the optimization algorithm used to initialize the array processing.
      % We can use 'ap' for alternating projectiona and 'grid' for grid search.
      array_param.init = 'ap';
      array_param.theta_guard = 1.5*pi/180;%(max(theta)-min(theta))/(4*Nc);
      
      % W: the widening factor. W=1 for narrowband.
      array_param.W = 1;
      % NB: number of subbands (used in WBMLE). NB should be <=length(bin_rng).
      
      % Each group of data has represents the number of snapshots per subband. So,
      % the total number of fast-time snapshots is
      % length(1:NB:length(bin_rng)-NB+1) * NB. To compare against MLE, the
      % number of fast-time snapshots in MLE must be equal to the
      % number of snapshots per subband, which is the n umber of data
      % groups length(1:NB:length(bin_rng)-NB+1) or ceil((length(bin_rng)-NB+1)/NB).
      array_param.NB = 1;
      dt = 1/param.fs;
      
      array_param.imp_resp.time_vec = -3*array_param.W*dt : dt/8 : 3*array_param.W*dt;
      array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / param.BW);
      
      % TO DO: CONSTRAINT THE SEARCH AT THE NEXT STEP BASED ON THE PREVIOUS STEP.
      %==========================================================================
      for idx = 1:array_param.Nsig
        array_param.doa_constraints(idx).method = 'fixed';
%         array_param.doa_constraints(idx).init_src_limits = [-20 40];
%         array_param.doa_constraints(idx).src_limits = [-20 40];
            array_param.doa_constraints(idx).init_src_limits = [min(theta) max(theta)]*180/pi;
            array_param.doa_constraints(idx).src_limits = [min(theta) max(theta)]*180/pi;
      end
      
      param.array_param = array_param;
      % clear array_param;
      
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
      
      results = crosstrack(param);
      
      if 1
        Nsig_true = results.Nsig_true;
        actual_doa  = results.actual_doa ;
        est_doa   = results.tomo.doa;  % Nt by Nsig by Nx matrix
        
        if 1%trial_idx == 1
          % Save a sample of the estimated DOA for plotting the surface later
          est_doa_save{method_idx}{test_idx,trial_idx} = est_doa;
          est_pwr_save{method_idx}{test_idx,trial_idx} = results.tomo.power;
        end
        
        % Convert actual_doa from cell array to matrix
        actual_doa_mtx = NaN(size(actual_doa,1),Nsig,size(actual_doa,2)); % Nt by Nsig by Nx matrix
        for col_idx = 1:size(actual_doa,2)
          for row_idx = 1:size(actual_doa,1)
            doa_t = NaN(Nsig,1);
            if ~isempty(actual_doa{row_idx,col_idx})
              if length(actual_doa{row_idx,col_idx}) == Nsig
                doa_t = actual_doa{row_idx,col_idx};
              elseif length(actual_doa{row_idx,col_idx}) > Nsig
                doa_t = actual_doa{row_idx,col_idx}(1:Nsig);
              else
                if actual_doa{row_idx,col_idx} < 0
                  doa_t(1) = actual_doa{row_idx,col_idx};
                else
                  doa_t(2) = actual_doa{row_idx,col_idx};
                end
              end
              actual_doa_mtx(row_idx,:,col_idx) = doa_t;
            end
          end
        end
      end
      
      if 0
        sim_data = results.sim_data;
        sim_data = squeeze(sim_data{1});
        sim_data_all{trial_idx} = sim_data;
      else
        sim_data_all = [];
      end
      
      % RMSE
      results.param.actual_doa = actual_doa_mtx;
      results.array_param.Nsig = Nsig;
      results.param.plot_doa_lims = [min_doa_lim max_doa_lim];
      
      if 0
        % Remove DoAs outside 3dB beamwidth
        for line_idx = 1:size(actual_doa,2)
          for bin_idx = 1:size(actual_doa,1)
            good_doa = actual_doa{bin_idx,line_idx}(actual_doa{bin_idx,line_idx}>min_doa_lim & actual_doa{bin_idx,line_idx}<max_doa_lim);
            %       actual_doa{bin_idx,line_idx} = good_doa;
            good_doa_idx = numel(good_doa);
            Nsig_true(bin_idx,line_idx) = good_doa_idx;
          end
        end
        
        %   est_doa(est_doa < min_doa_lim) = NaN;
        %   est_doa(est_doa > max_doa_lim) = NaN;
      end
    end
    
    if strcmp(test_type,'Number of snapshots')
      test_param(test_idx) = Nsnaps;
    elseif strcmp(test_type,'SNR (dB)')
      test_param(test_idx) = snr_dB(test_idx);
    elseif strcmp(test_type,'Fractional BW')
      test_param(test_idx) = Bf(test_idx);
    end
  end
end

%% Estimated RMSE for each DOA and each method separately 
threshold_val = 3;
error_mle   = [];
error_music = [];
rmse_music  = [];
rmse_mle    = [];
rmse_tests_mean = [];

for method_idx = 1:length(methods_list)
  method_i = methods_list(method_idx);
  for test_i = 1:length(test_param)
    for run_i = 1:Ntrials
      if method_i == 2
        % MUSIC
        tmp_est_doa_music = est_doa_save{1}{test_i,run_i};
        tmp_est_doa_music = tmp_est_doa_music(~isnan(tmp_est_doa_music));
        
        est_doa_music = nan(size(ang));
        est_doa_music(1:length(tmp_est_doa_music)) = tmp_est_doa_music;
        error_music(run_i,:) = abs(est_doa_music - ang);
        
      elseif method_i == 7
        % MLE
        tmp_est_doa_mle = est_doa_save{2}{test_i,run_i};
        tmp_est_doa_mle = tmp_est_doa_mle(~isnan(tmp_est_doa_mle));
        
        est_doa_mle = nan(size(ang));
        est_doa_mle(1:length(tmp_est_doa_mle)) = tmp_est_doa_mle;
        error_mle(run_i,:) = abs(est_doa_mle - ang);
      end
    end
    
    % RMSE
    if method_i == 2
      % MUSIC
      threshold = threshold_val *std(error_music(~isnan(error_music(:))));
      
      error_music(error_music>threshold) = 0;
      rmse_music(test_i,:) = sqrt(nanmean(error_music.^2,1))*180/pi;
      
    elseif method_i == 7
      % MLE
      threshold = threshold_val *std(error_mle(~isnan(error_mle(:))));
      
      error_mle(error_mle>threshold) = 0;
      rmse_mle(test_i,:) = sqrt(nanmean(error_mle.^2,1))*180/pi;
    end
  end
  % Average RMSE over all targets of each range-bin
  if method_i == 2
  rmse_music_rbin = nanmean(rmse_music,2);
  rmse_tests_mean(:,1) = rmse_music_rbin;
  elseif method_i == 7
  rmse_mle_rbin   = nanmean(rmse_mle,2);
  rmse_tests_mean(:,2) = rmse_mle_rbin;
  end
end


if 0
  % Save
  sim_param.fc                = param.fc;
  sim_param.BW                = param.BW;
  sim_param.Nc                = Nc;
  sim_param.Nsnap             = Nsnaps;
  % sim_param.SNR             = snr_db;
  sim_param.Nruns             = Ntrials;
  sim_param.actual_doa_deg    = ang*180/pi;
  sim_param.tx_window         = 'hanning';
  sim_param.beamwidth_3dB_deg = (max(theta_3dB)-min(theta_3dB))*180/pi;
  sim_param.note = 'rmse has dimension of Nt by Nx by Ntests by Nmethods. rmse_tests_mean has dimension of Ntests by Nmethods (this is what you should plot)';
  
  out_fn_dir = '/users/mohanad/IceSheetProject/MUSIC-MLE-PF comparison/results/';
  out_fn_name = 'MLE-MISUC_WB_SNR';%sprintf('MLE-MUSIC_RMSEvsSNR_Nsnaps%d_fc%dG_BW%dM',Nsnaps,param.fc/1e9,param.BW/1e6);
  
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn = fullfile(out_fn_dir,out_fn_name);
  save([out_fn '.mat'],'error_all','rmse_tests','rmse_tests_mean','test_param','sim_param')
end


%% To plot Slice model
% The 'Ground-truth surface' is the surface that we modeled based on the
% papameters provided in surface_param (our model for the surface). The
% 'Actual surface' is the surface generated based on the estimated
% DoA and the range to each range-bin, which depends on the time gate (t0
% and t1) and the sampling frequency==> range=2*(t1:1/fs:t0)*c. Thus, the
% two surfaces may not coincide. I think better to call it the ' estimated
% surface' rather than the 'Actual surface'.
slice = 1;
if 0
  surface_z = results.z_grid-results.param.monte.target_param{1}.z.mean;
  surface_z_groung_truth =results.surf_model.z-results.param.monte.target_param{1}.z.mean;
  figure(5); clf;
  plot(results.surf_model.y, surface_z_groung_truth(:,slice),'b');
  hold on
  plot(results.surf_model.y, surface_z(:,slice),'r');
  xlim([-1500 1500])
  ylim([-200 250])
  title('Slice - surface model');
  xlabel('Cross-track (m)');
  ylabel('WGS84-Elevation (m)');
  hold off
  legend('Ground-truth surface','Actual (i.e. estimated) surface');
  grid on
end
%

if 0
  % Plot beam pattern
  figure(33);clf;
  subplot(121)
  plot(theta_RP*180/pi,RP_dB);
  grid on
  xlabel('theta [deg.]')
  ylabel('Relative power pattern [dB]')
  title('Radiation pattern')
  xlim([theta_RP(1) theta_RP(end)]*180/pi)
  ylim([-100 0])
  %   ylim([-30 0]);
  
  % Plot 3dB beampattern
  figure(33);
  subplot(122)
  plot(theta_3dB*180/pi,RP_3dB);
  grid on
  xlabel('theta [deg.]')
  title('3dB Radiation pattern')
  xlim([theta_3dB(1) theta_3dB(end)]*180/pi)
end

if 0
  RMSE_mean = RMSE.RMSE_mean;
  figure(44);clf;
  %    scatter(actual_doa{:,slice}*180/pi,RMSE_mean(:,slice)*180/pi,20,'fill')
  scatter(mean(squeeze(mean(est_doa,2)),2)*180/pi,RMSE_mean(:,slice)*180/pi,20,'fill')
  xlim([0 +50])
  xlabel('True DoA (deg.)')
  ylabel('RMSE (deg.)')
  title('2D simulator')
  grid on
end

%
% Prepare the estimated and actual DoA matrices so that the maximum Nsig = 2.
actual_doa_tmp = NaN(size(actual_doa_mtx,1),Nsig);
est_doa_tmp    = NaN(size(est_doa,1),Nsig,length(methods_list));
pwr_tmp        = NaN(size(est_doa,1),Nsig,length(methods_list));

%   reqd_test_idx = find(snr_dB==10);
[~, reqd_test_idx] = min(abs(snr_dB-10));
reqd_trial_idx = 1;
%   reqd_test_idx = 1;
for method_idx = 1:length(methods_list)
  for bin_idx = 1:size(est_doa,1)
    len_est = size(est_doa,2);
    len_act = size(actual_doa_mtx,2);
    len     = min(len_est,len_act);
    
    if len == 0
      % Do nothing
    elseif len == 1
      est_doa_tmp(bin_idx,1,method_idx) = est_doa_save{method_idx}{reqd_test_idx,reqd_trial_idx}(bin_idx,1,slice);
      est_doa_tmp(bin_idx,2,method_idx) = NaN;
      
      actual_doa_tmp(bin_idx,1)         =  actual_doa_mtx(bin_idx,1,slice);
      actual_doa_tmp(bin_idx,2)         = NaN;
      
      pwr_tmp(bin_idx,1,method_idx)     = 10*log10(est_doa_save{method_idx}{reqd_test_idx,reqd_trial_idx}(bin_idx,1,slice));
      pwr_tmp(bin_idx,1,method_idx)     = NaN;
    elseif len == 2
      est_doa_tmp(bin_idx,:,method_idx) = est_doa_save{method_idx}{reqd_test_idx,reqd_trial_idx}(bin_idx,1:2,slice);
      actual_doa_tmp(bin_idx,:)         = actual_doa_mtx(bin_idx,1:2,slice);
      pwr_tmp(bin_idx,:,method_idx)     = 10*log10(est_doa_save{method_idx}{reqd_test_idx,reqd_trial_idx}(bin_idx,1:2,slice));
    else
      est_doa_tmp(bin_idx,:,method_idx) = est_doa_save{method_idx}{reqd_test_idx,reqd_trial_idx}(bin_idx,1:2,slice);
      actual_doa_tmp(bin_idx,:)         = actual_doa_mtx(bin_idx,1:2,slice);
      pwr_tmp(bin_idx,:,method_idx)     = 10*log10(est_doa_save{method_idx}{reqd_test_idx,reqd_trial_idx}(bin_idx,1:2,slice));
    end
  end
end

% Range-bin indices for each DOA
%   range_bin_idxs_DOA = find(~isnan(doa_actual_tmp))';

% Range Bin vs DOA
idx_1_max = find(~isnan(actual_doa_tmp(:,1)),1,'last');
idx_2_max = find(~isnan(actual_doa_tmp(:,2)),1,'last');

idx_1_min = find(~isnan(actual_doa_tmp(:,1)),1,'first');
idx_2_min = find(~isnan(actual_doa_tmp(:,2)),1,'first');

if ~isempty(idx_1_max) && ~isempty(idx_2_max)
  max_y_lim = max(idx_1_max,idx_2_max);
elseif isempty(idx_1_max) || isempty(idx_2_max)
  if isempty(idx_1_max)
    max_y_lim = idx_2_max;
  else
    max_y_lim = idx_1_max;
  end
else
  max_y_lim = Nt-max(array_param.bin_rng);
end

if ~isempty(idx_1_min) && ~isempty(idx_2_min)
  min_y_lim = min(idx_1_min,idx_2_min);
elseif isempty(idx_1_min) || isempty(idx_2_min)
  if isempty(idx_1_min)
    min_y_lim = idx_2_min;
  else
    min_y_lim = idx_1_min;
  end
else
  min_y_lim = Nt_st+max(array_param.bin_rng);
end



% if isempty(idx_1)
%   max_y_lim = idx_2+array_param.dbin;
% elseif isempty(idx_2)
%   max_y_lim = idx_1+array_param.dbin;
% else
%   max_y_lim = max(idx_1+array_param.dbin,idx_2+array_param.dbin);
% end
%
% idx_1_min = find(~isnan(actual_doa_tmp(:,1)),1,'first');
% idx_2_min = find(~isnan(actual_doa_tmp(:,2)),1,'first');
% if isempty(idx_1)
%   min_y_lim = idx_2+array_param.dbin;
% elseif isempty(idx_2)
%   min_y_lim = idx_1+array_param.dbin;
% else
%   min_y_lim = max(idx_1+array_param.dbin,idx_2+array_param.dbin);
% end

lim_guard = round((max_y_lim-min_y_lim)/10);
lim_guard = max(lim_guard,5);

%   Color  = {'b','c','g','m','k',[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560]};
Color  = {'b','r','g','m','k',[0.4940, 0.1840, 0.5560]};
Marker = {'x','*','s','+','^','<','>'};
LineWidth = [3 1];

figure(6);clf
hold on
figure(7);clf
%   figure(8);clf
figure(9);clf

for method_idx = 1:length(methods_list)
  method_num = methods_list(method_idx);
  if method_num == 2
    method.name = 'MUSIC';
  elseif method_num == 7
    method.name = 'MLE';
  elseif method_num == 9
    method.name = 'WBMLE';
  end
  
  % Sample surface
  % -----------------
  figure(6);
  %     subplot(1,length(methods_list),method_idx)
  for signal_idx = 1:Nsig
    doa_tmp = [NaN(array_param.dbin,1);est_doa_tmp(:,signal_idx,method_idx);NaN(array_param.dbin,1)]*(180/pi);
    doa_actual_tmp = [NaN(array_param.dbin,1);actual_doa_tmp(:,signal_idx);NaN(array_param.dbin,1)]*(180/pi);
    if method_idx == 1
      h_actual = scatter(doa_actual_tmp,1:length(doa_actual_tmp), 40,'+','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
    end
    hold on
    if method_idx == 1
      % scatter(est_doa_tmp(:,1)*(180/pi),results.array_param.bins, 20 , pwr_tmp(:,1),'fill','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1.5);
      h1 = scatter(doa_tmp,1:length(doa_tmp), 10*LineWidth(method_idx),Marker{method_idx}','MarkerFaceColor',Color{method_idx},'MarkerEdgeColor',Color{method_idx},'LineWidth',LineWidth(method_idx));
    elseif method_idx == 2
      h2 = scatter(doa_tmp,1:length(doa_tmp), 10*LineWidth(method_idx),Marker{method_idx}','MarkerFaceColor',Color{method_idx},'MarkerEdgeColor',Color{method_idx},'LineWidth',LineWidth(method_idx));
    end
  end
  set(gca,'Ydir','reverse')
  if method_idx == length(methods_list)
    % xlim([min_doa-5*(pi/180) max_doa+5*(pi/180)]*(180/pi))
    xlim([min_doa_lim-5*pi/180 max_doa_lim+5*pi/180]*180/pi)
    ylim([min_y_lim-lim_guard max_y_lim+lim_guard])
    %       ylim([1 results.array_param.bins(end)])
    xlabel('DOA (deg.)');
    ylabel('Range bin');
    title('Sample surface')
  end
  %     title(sprintf('%s: Sample surface for slice # %u',method.name,slice))
  
  if method_idx == length(methods_list)
    if length(methods_list) == 1
      legend_est_doa = sprintf('%s',method.name);
      legend('Actual DoA',legend_est_doa ,'Location','best')
    elseif length(methods_list) == 2
      legend_est_doa_1 = sprintf('MUSIC');
      legend_est_doa_2 = sprintf('MLE');
      legend([h_actual  h1  h2],{'Actual DoA',legend_est_doa_1,legend_est_doa_2} ,'Location','best')
    end
    
  end
  grid on
  
  % True vs estimated DoA
  % ------------------------
  figure(7);
  %     subplot(1,length(methods_list),method_idx)
  if method_idx == 1
    h3 = scatter(actual_doa_tmp(:,1)*(180/pi),est_doa_tmp(:,1,method_idx)*(180/pi),40, '+', 'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',2);
    hold on
    scatter(actual_doa_tmp(:,2)*(180/pi),est_doa_tmp(:,2,method_idx)*(180/pi),40, '+', 'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',2);
  elseif method_idx == 2
    h4 = scatter(actual_doa_tmp(:,1)*(180/pi),est_doa_tmp(:,1,method_idx)*(180/pi),30, '*', 'MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',2);
    hold on
    scatter(actual_doa_tmp(:,2)*(180/pi),est_doa_tmp(:,2,method_idx)*(180/pi),30, '*', 'MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',2);
  end
  
  if method_idx == length(methods_list)
    % xlim([min_doa-5*(pi/180) max_doa+5*(pi/180)]*(180/pi))
    xlim([min_doa_lim-5*pi/180 max_doa_lim+5*pi/180]*180/pi)
    xlabel('True DoA (deg.)')
    ylabel('Estimated DoA (deg.)')
    title('Sample true vs estimated DOA')
    grid on
  end
  %     title(sprintf('%s: Sample surface for slice # %u',method.name,slice))
  
  if method_idx == length(methods_list)
    if length(methods_list) == 1
      legend_est_doa = fprintf('Est. DOA: %s',method.name);
      legend('Actual DoA',legend_est_doa ,'Location','best')
    elseif length(methods_list) == 2
      legend([h3  h4],{'MUSIC','MLE'} ,'Location','best')
    end
  end
  
  if 0
    % DoA vs RMSE
    figure(44);clf;
    scatter(theta*180/pi,rmse(rbins,slice),20,'fill');
    xlim([min_doa_lim-5*pi/180 max_doa_lim+5*pi/180]*180/pi)
    xlabel('True DoA (deg.)')
    ylabel('RMSE (deg.)')
    title(strcat('2D simulator: ','Slice# ',num2str(slice)))
    grid on
  end
  
  if 0
    test_idx = 1;
    % RMSE plots
    figure(8 + method.list);
    subplot(311)
    imagesc(rmse_tests.rmse(:,:,test_idx,method_idx)*180/pi)
    h_colorbar = colorbar;
    set(get(h_colorbar,'YLabel'),'String','RMSE (deg)');
    xlabel('Range line');
    ylabel('Range bin');
    title(fprintf('%s: RMSE over %u runs',method.name,Ntrials))
    ylim([min_y_lim-5 max_y_lim+5])
    caxis([min(rmse(:)) max(rmse(:))]*180/pi)
    
    subplot(312)
    scatter(1:length(rmse_tests.avg_rbin_rmse(:,test_idx,method_idx)),rmse_tests.avg_rbin_rmse(:,test_idx,method_idx)*180/pi,20,'fill', 'MarkerFaceColor','b');
    xlabel('Range-bin')
    ylabel('RMSE (deg)')
    title(fprintf('%s: Average RMSE (over range-lines)',method.name))
    grid on
    
    subplot(313)
    scatter(1:length(rmse_tests.avg_rline_rmse(:,test_idx,method_idx)),rmse_tests.avg_rline_rmse(:,test_idx,method_idx)*180/pi,20,'fill', 'MarkerFaceColor','b');
    xlabel('Range-line')
    ylabel('RMSE (deg)')
    title(fprintf('%s: Average RMSE (over range-bins)',method.name))
    grid on
  end
  
  % RMSE vs TEST_PARAM
  if exist('test_param','var') && ~isempty(test_param)
    figure(9);
    %       subplot(1,length(methods_list),method_idx)
    if method_idx == 1
      if strcmp(test_type,'Number of snapshots')
        h5 = semilogx(test_param.',rmse_tests_mean(:,method_idx),'+b','LineWidth',2);
      else
      h5 = scatter(test_param.',rmse_tests_mean(:,method_idx), 40,'+','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',2);
      %       h5 = plot(test_param,rmse_tests_mean(:,method_idx)*180/pi,'x','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',2);
      end
      hold on
    elseif method_idx == 2
      if strcmp(test_type,'Number of snapshots')
        h6 = semilogx(test_param.',rmse_tests_mean(:,method_idx),'*r','LineWidth',2);
      else
      h6 = scatter(test_param.',rmse_tests_mean(:,method_idx), 30,'*','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',2);
      %       h6 = plot(test_param,rmse_tests_mean(:,method_idx)*180/pi,'x','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',2);
      end
    end
    
    if length(test_param) > 1 && test_param(2) > test_param(1)
%       xlim([test_param(1)  test_param(end)])
    else
      % No xlim
    end
    
    if method_idx == length(methods_list)
      xlabel(test_type)
      ylabel('RMSE (deg.)')
      %         title(sprintf('%s: %s vs average RMSE (over range-lines/bins)',method.name,test_type))
      grid on
    end
    
    if method_idx == length(methods_list)
      if length(methods_list) == 1
        %           legend([h5  h6],{'MUSIC','MLE'} ,'Location','best')
      elseif length(methods_list) == 2
        legend([h5  h6],{'MUSIC','MLE'} ,'Location','best')
      end
    end
    
  end
end

%% CRLB
if 0
  crb_DOA = [0 30];
  crb_param.src.fc                   = param.fc;
  crb_param.src.BW                   = param.BW;
  crb_param.method.wb_td.widening_factor = 1;
  crb_param.method.wb_fd.filter_banks    = 1;
  crb_param.method.list                  = methods_list;
  crb_param.src.lever_arm.fh             = @sim.lever_arm_example;
  crb_param.src.lever_arm.args           = {[],1,[1:Nc],[0; c/crb_param.src.fc/2; 0]};
  
  Color = {'k','c'};
  Marker = {'s','^'};
  for doa_idx = 1:length(crb_DOA)
    DOA = crb_DOA(doa_idx); % In degrees
    if 1
      crb_param.monte.ASNR  = repmat(test_param.',[1 length(DOA)]); % In dB
      crb_param.monte.SNR   = crb_param.monte.ASNR; % Array SNR
      % param.monte.SNR   = param.monte.ASNR - 10*log10(Nchan); % Channel SNR
      num_tests         = size(crb_param.monte.SNR,1);
      crb_param.monte.DOA   = repmat(DOA,[num_tests 1]);
      crb_param.monte.Nsnap = repmat(Nsnaps,[num_tests 1]);
    elseif 0
      crb_param.monte.Nsnap = repmat(test_param.',[1 length(DOA)]);
      num_tests             = size(crb_param.monte.Nsnap,1);
      crb_param.monte.ASNR  = repmat(snr_db,[num_tests length(DOA)]); % In dB
      crb_param.monte.SNR   = crb_param.monte.ASNR; % Array SNR
      % param.monte.SNR   = param.monte.ASNR - 10*log10(Nchan); % Channel SNR
      
      crb_param.monte.DOA   = repmat(DOA,[num_tests 1]);
      
    end
    desired_src = 1;
    [CRB_angular , CRB_spatial] = crb(crb_param, desired_src);
    
    figure(9);
    hold on
    plot_name  = sprintf('h%d',doa_idx+6);
    plot_name = scatter(test_param,(sqrt(CRB_angular(desired_src,:))*180/pi).',20,Marker{doa_idx},'MarkerFaceColor',Color{doa_idx},'MarkerEdgeColor',Color{doa_idx},'LineWidth',2);
    xlim([test_param(1)  test_param(end)])
    
    legend({'MUSIC','MLE','CRLB: 0^\circ target','CRLB: 30^\circ target'} ,'Location','best')
    
  end
end

%%
if ~isempty(sim_data_all)
  % Estimate SNR per channel (for comparing against 1D simulator results)
  xx = cat(3,sim_data_all{:});
  n_pwr  = abs(xx(max_y_lim+1:end,:,:)).^2;
  n_pwr  = mean(n_pwr(:));
  s_pwr  = abs(xx(min_y_lim:max_y_lim,:,:)).^2;
  s_pwr  = mean(s_pwr(:));
  snr_db = 10*log10(s_pwr/n_pwr);
  
  sprintf('\nEstimated SNR is: %2.2f dB\n',snr_db)
  %     fprintf('\nFoV is: %2.0f to %2.0f\n',min(theta)*180/pi,max(theta)*180/pi)
end

toc

if param.debug_level >= 3
  return
end


%% Notes
% 1. The FoV is limited by the maximum range (i.e. the time gate) and the
% fligt height. So, to increase FoV, increase the time gate.

