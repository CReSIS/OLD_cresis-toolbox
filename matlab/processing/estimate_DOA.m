% script estimate_DOA.m
%
% Function for estimating direction of arrival

% roll_to_ant_mapping(ANTENNA) = index into surf_vals containing this
% antenna's data (transmit and receive with the same antenna)
% 
% surf_bins = The relative range bin into surf_vals that we will use for 
%             extracting values from


physical_constants


%% 2013_Antarctica_Basler Setup
if 1
  
  % fn = output from coh_noise surf tracker
  
  % 
  if 0
    fn    =  '/cresis/snfs1/dataproducts/ct_data/rds/2013_Antarctica_Basler/CSARP_noise/surf_20131216_05.mat';
    roll_to_ant_mapping = [17:24];
    surf_bins           = 6;
    bin_rng             = [-1:1];
    rline_rng           = [-5:5];
    ref_ant             = 4;
  end
  
  if 1
    fn    =  '/cresis/snfs1/dataproducts/ct_data/rds/2013_Antarctica_Basler/CSARP_noise/surf_20131216_05.mat';
    roll_to_ant_mapping = [25:32];
    surf_bin_offset     = -6;
    bin_rng             = 0;
    bin_rng             = [-1:1];
    rline_rng           = [-3:3];
    ref_ant             = 3;
  end
  
  % Multibeam Data
  if 0
    fn        = '/cresis/snfs1/dataproducts/ct_data/rds/2013_Antarctica_Basler/CSARP_noise/surf_20131216_08.mat'; 
    roll_to_ant_mapping = [1:8]; 
    surf_bins           = 17;
    bin_rng             = [-1:1];
    rline_rng           = [-3:3];
    ref_ant             = 4;
  end
  %   output_fn = '/cresis/snfs1/dataproducts/ct_data/rds/2013_Antarctica_Basler/CSARP_noise/sv_table_2013_Antarctica_Basler.mat';
  
  rlines      = [];
  % Antenna phase reference channel (usually a center element)
  fc          = (200e6 + 450e6)/2; % Center frequency
  
  Nsv         = 48; % coarse grid for initialization
  Nsig        = 1;
  
end

%% Load data
data = load(fn);
zero_surf_bin = 1-data.param_analysis.analysis.surf.bin_rng(1);
surf_bins = zero_surf_bin + surf_bin_offset;

if 0
  close all
  %% DEBUG
  figure(1); clf;
  subplot(3,1,1);
  plot(lp(data.surf_vals(surf_bins,:,1)),'.')
  a1 = gca;
  subplot(3,1,2);
  plot(angle(data.surf_vals(surf_bins,:,4) .* conj(data.surf_vals(surf_bins,:,1)) ),'.')
  a2 = gca;
  subplot(3,1,3);
  plot(data.roll);
  a3 = gca;
  linkaxes([a1 a2 a3],'x');
  figure(2); clf;
  imagesc(lp(data.surf_vals(:,:,1)));
  return;
end

%% Retrack surface
if 1
  ref_idx = roll_to_ant_mapping(ref_ant);
  
  threshold = lp(mean(abs(data.surf_vals(1,:,ref_idx)).^2)) + 4;
  surf_bin = NaN*zeros(1,size(data.surf_vals,2));
  ml_data = lp(fir_dec(abs(data.surf_vals(:,:,ref_idx)).^2,ones(1,5)/5,1));
  for rline = 1:size(data.surf_vals,2)
    cur_threshold = max([ml_data(1,rline)+7; ml_data(:,rline)-13]);
    tmp = find(ml_data(:,rline) > cur_threshold,1);
    if ~isempty(tmp)
      [~,max_offset] = max(ml_data(tmp+(0:2),rline));
      tmp = tmp-1 + max_offset;
      surf_bin(rline) = tmp;
    end
  end
  
  figure(1); clf;
  imagesc(ml_data);
  hold on;
  plot(surf_bin);
  
  for rline = 1:size(data.surf_vals,2)
    if ~isnan(surf_bin(rline))
      data.surf_vals(:,rline,:) = circshift(data.surf_vals(:,rline,:),[zero_surf_bin-surf_bin(rline) 0 0]);
    end
  end
  
  ml_data = lp(fir_dec(abs(data.surf_vals(:,:,ref_idx)).^2,ones(1,5)/5,1));
  % Check to make sure surface is flat
  figure(2); clf;
  imagesc(ml_data);
  
  keyboard
end

%% Setup phase centers
wf_adc_list   = data.param_analysis.analysis.imgs{1}(roll_to_ant_mapping,:);

lever_arm_param.season_name = data.param_records.season_name;
lever_arm_param.radar_name  = ct_output_dir(data.param_records.radar_name);
lever_arm_param.gps_source  = data.param_records.gps_source;

phase_centers = [];

for wf_adc_idx = 1:length(wf_adc_list)
  
  wf  = abs(wf_adc_list(wf_adc_idx,1));
  adc = wf_adc_list(wf_adc_idx,2);
  phase_centers(:,end+1) = lever_arm(lever_arm_param,data.param_records.radar.wfs(wf).tx_weights, ...
    data.param_records.radar.wfs(wf).rx_paths(adc));
end

%% Setup for WDOA

% Compute maximum propagation delay between any two phase centers
for pc_idx = 1:size(phase_centers,2)
  for pc_comp_idx = pc_idx+1:size(phase_centers,2)
    pc_dist(pc_idx,pc_comp_idx) ...
      = sqrt(sum(abs(phase_centers(:,pc_idx) - phase_centers(:,pc_comp_idx)).^2));
  end
end
max_array_dim = max(max(pc_dist));
tau_max       = 2*max_array_dim/c;

% Compute range impulse response
Nt            = size(data.surf_vals,1);
Nx            = size(data.surf_vals,2);
Nc            = size(phase_centers,2);
W             = ceil(tau_max/data.wfs(wf).dt);
Mt            = 10;
Hwin          = data.param_analysis.get_heights.ft_wind(length(data.wfs(wf).time));
Hwin          = interpft(real(ifft(ifftshift(Hwin))),Mt*length(Hwin));
Hwin_num_samp = 2 * Mt * W;
imp_resp.vals = fftshift(Hwin([1:1+Hwin_num_samp, end - Hwin_num_samp + 1:end]));
imp_resp.time_vec = data.wfs(wf).dt/Mt * (-Hwin_num_samp:Hwin_num_samp);

% Setup wideband param
wb_param.h              = imp_resp.vals(:);
wb_param.t0             = imp_resp.time_vec(1);
wb_param.dt             = imp_resp.time_vec(2) - wb_param.t0;
wb_param.fs             = data.wfs(wf).fs;
wb_param.fc             = fc;
wb_param.Nsig           = 1;
wb_param.options        = optimset('Display','off');
wb_param.y_pc           = phase_centers(2,:).';
wb_param.z_pc           = phase_centers(3,:).';
wb_param.src_limits{1}  = [-pi/2 pi/2];

%% Create ideal steering vectors
[theta,sv]              = array_proc_sv(Nsv,fc,phase_centers(2,:).',phase_centers(3,:).');
SV                      = fftshift(sv,2);
theta                   = fftshift(theta);

%% Finish setting up control structures for DOA estimators
wb_param.theta          = theta;
nb_param.y_pc          = phase_centers(2,:).';
nb_param.z_pc          = phase_centers(3,:).';
nb_param.theta         = theta;
nb_param.SV            = SV;
nb_param.Nsig          = 1;
nb_param.fs            = data.wfs(wf).fs;
nb_param.fc            = fc;
nb_param.options       = wb_param.options;
nb_param.src_limits{1} = [-pi/2 pi/2];

%% Setup for looping over along-track and defining conditions for edges
% Preallocate memory for estDOAs. Nx by 3 matrix used to store DOA
% estimates at each along-track position from MUSIC, MLE and WDOA
estDOAs     = nan(numel(data.roll),3); 
start_line  = max(1 - rline_rng);
stop_line   = min(Nx - rline_rng);
start_bin   = 1-min(bin_rng - floor((W-1)/2));

%% Loop over slow time and estimate DOA
for lineIdx = start_line:stop_line
  
  line    = lineIdx;
%   bin     = surf_bin(line);
  bin = surf_bins;
  
%   start_bin = min(bin + bin_rng);
%   if start_bin < 1
%     bin = 1 - min(bin_rng);
%   end
    
  if ~mod(lineIdx-1,50)
    fprintf('Line %d of %d (%s)\n',line,Nx,datestr(now));
  end
  
  if ~isnan(bin) && bin >= start_bin
    % Estimate covariance matrix for narrowband methods
    dataSample_nb   = data.surf_vals(bin+bin_rng,line+rline_rng,roll_to_ant_mapping);
    dataSample_nb   = reshape(dataSample_nb,[length(bin_rng)*length(rline_rng) Nc]);
    dataSample_nb   = dataSample_nb.';
    Rxx             = (1/size(dataSample_nb,2)) * (dataSample_nb * dataSample_nb');
    
    nb_param.Rxx   = Rxx; 
    
    % MUSIC
    % =======================================================================
    if 1
      music_param   = nb_param;
      % Evaluate pseudospectrum over a coarse grid and estimate peaks
      theta0_music  = music_initialization(Rxx,nb_param);
      [doa_music, fval,exitflag]  = fminsearch(@(X) music_cost_function(X,nb_param), theta0_music, nb_param.options);
      estDOAs(line,1) = doa_music;
    else
      
      [V,D]                 = eig(Rxx);
      eigenVals             = diag(D);
      [eigenVals noiseIdxs] = sort(eigenVals);
      noiseIdxs             = noiseIdxs(1:end-Nsig);
      Sarray                = 1./(mean(abs(SV(:,:)'*V(:,noiseIdxs)).^2,2));
      [~,doa_idx]           = max(Sarray);
      estDOAs(line,1)       = theta(doa_idx);
    end
    
    % MLE
    % =======================================================================
    theta0_mle                = mle_initialization(Rxx,nb_param);
    [doa_mle, fval,exitflag]  = fminsearch(@(theta_hat) mle_cost_function(theta_hat,nb_param), theta0_mle,nb_param.options);
    estDOAs(line,2)           = doa_mle;
    
    % Wideband DOA
    % =======================================================================
    
    % Estimate space-time covariance matrix
    % ----------------------------------------------------------------
    dataSample        = [];
    for W_offset = -floor((W-1)/2):floor((W-1)/2)
      offset_bin      = bin + W_offset;
%       bins = offset_bin + bin_rng;
%       bins = bins(bins>0);
%       dataSample_tmp  = double(data.surf_vals(bins,line+rline_rng,roll_to_ant_mapping));
      dataSample_tmp  = double(data.surf_vals(offset_bin + bin_rng,line+rline_rng,roll_to_ant_mapping));
      dataSample_tmp  = reshape(dataSample_tmp,[length(bin_rng)*length(rline_rng) Nc]).';
      dataSample      = cat(1,dataSample,dataSample_tmp);
    end
    
    DCM           = (1/size(dataSample,2))*dataSample*dataSample';
    wb_param.DCM  = DCM;
    theta0        = wb_initialization(DCM,wb_param);
    
    [doa_wb,Jval,exitflag,OUTPUT] = ...
      fminsearch(@(theta_hat) wb_cost_function(theta_hat,wb_param), theta0,wb_param.options);
    estDOAs(line,3) = doa_wb;
    
  end
end

%% Plot INS roll and DOA estimated Rolls
% =========================================================================
figure(6);clf
plot(data.roll.*180/pi,'.k','LineWidth',2);
hold on
grid on
plot(estDOAs.*180/pi,'.','LineWidth',2);
legend('INS','MUSIC','MLE','WDOA','Location','Best')


%% Estimate RMS error versus DOA
% =========================================================================
[trueDOAs, sort_idxs]  = sort(data.roll);

sort_idxs = sort_idxs(:);
trueDOAs  = trueDOAs(:).*180/pi;
estDOAs   = estDOAs.*(180/pi);
sortDOAs  = estDOAs(sort_idxs,:);

% Squared error (e2 sorted to correspond with trueDOAs)
% -------------------------------------------------------------------------
e2        = abs(estDOAs(sort_idxs,:) - repmat(trueDOAs,1,size(estDOAs,2))).^2;

[roll_binned,~,roll_idxs] = unique(round(trueDOAs));

for roll_idx = 1:length(roll_binned)
  sample_size(roll_idx) = sum(roll_idxs == roll_idx);
  e2_vals = e2((roll_idxs == roll_idx),:);
  binnedDOAs = sortDOAs((roll_idxs == roll_idx),:);
  theta_hat(roll_idx,:) = nanmean(binnedDOAs,1);
  std_dev(roll_idx,:)   = nanstd(binnedDOAs,1,1);
  RMSE(roll_idx,:)      = (nanmean(e2_vals,1)).^0.5;
end


roll_binned = roll_binned(:);

figure(7);
clf;
plot(roll_binned,RMSE,'x')
grid on
xlabel('Roll ( \circ )')
ylabel('RMSE ( \circ )')
title('Measured Error of DOA Estimators')
legend('MUSIC','MLE','WDOA')


figure(8);
clf;
errorbar(repmat(roll_binned,1,size(theta_hat,2)),theta_hat,std_dev,'x');
grid on
hold on
xlabel('Roll ( \circ )')
ylabel('DOA ( \circ )')
legend('MUSIC','MLE','WDOA')


 %% Expected DCM
 if 0
   test_theta = data.roll(line); % angle you want to test
   t     = imp_resp.time_vec;
   h     = imp_resp.vals;
   uy    = sin(test_theta).'; % make directional sines and cosines into row vecs
   uz    = cos(test_theta).';
   tau   = (2/c)*(wb_param.y_pc*uy - wb_param.z_pc*uz);
%    Nc     = lenghth(wb_param.y_pc);
   
   shift_vector = [0*ones(8,1);-1*ones(8,1);-2*ones(8,1)];
   columns = 1:size(wb_param.DCM,2);
   for c_idx = columns
     
     if c_idx > 8 && c_idx < 17
       shift_vector = [0*ones(8,1);-1*ones(8,1);-2*ones(8,1)];
       shift_vector = shift_vector + 1;
     end
     
     if c_idx > 16 && c_idx <= 24;
       shift_vector = [0*ones(8,1);-1*ones(8,1);-2*ones(8,1)];
       shift_vector = shift_vector + 2;
     end
     ref_index = mod(c_idx - 1,8) + 1;
     tau_sub   = tau - repmat(tau(ref_index,:),8,1);
     tau_vec   = repmat(tau_sub,3,1)- repmat((shift_vector ./ wb_param.fs),1,size(tau,2));
     H_i       = interp1(t,h,tau_vec);
     A         = (H_i).*exp(1i*2*pi*wb_param.fc*repmat(tau_sub,3,1)); %%% original
     DCM_model(:,c_idx) = A;
   end
   figure(25);clf
   imagesc(10*log10(abs(DCM_model)))
   
   figure(26);clf
   imagesc(angle(DCM_model))
 end

