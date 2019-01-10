function [sim_data] = crosstrack_data(param,target_model)
% sim_data = sim.crosstrack_data(param,target_model)
%
% Creates 2-D simulated radar data specifically for testing array
% processing algorithms. This is equivalent to SAR-processed
% data where there is no contamination or interference from other along
% track targets (i.e. resolution in along track is perfect).
% Note that spherical spreading loss and attenuation are not accounted for.
% Noise is assumed to be complex gaussian and uncorrelated across snapshots.
%
% TODO: individual element complex beam patterns, currently ideal and
%   isotropic sensors are assumed
% TODO: correlated noise in time and between channels, currently additive
%   white Gaussian noise is assumed in time and between channels
% TODO: support physical optics facet target model, currently target model
%   is isotropic scatterers
%
% param: structure describing how to create the simulated data
%  .src = structure describing signal source characteristics
%    .fs = sampling frequency (Hz)
%    .t0 = start time of range gate (sec)
%    .t1 = stop time of range gate (sec)
%    .y_pc: Nc length vector of sensor positions in cross track with
%           positive pointing to the left
%    .z_pc: Nc length vector of sensor position in elevation with positive
%           pointing up
%    .ft_func: function handle to function generating fast time signal
%              Duration should fall within -(t0+t1)/2 to (t0+t1)/2 and take
%              time in seconds as the only argument.
%    .ft_wind: function handle to function generating fast time frequency
%              domain window. Should take number of samples for window.
%    .noise_power: normalized (to target RCS) noise power
% target_model: structure describing target positions for each snapshot
%  .z: Ns by Nx array of target z-positions in meters
%  .y: Ns by Nx array of target y-positions in meters
%  .rcs: Ns by Nx array of target radar cross sections (normalized)
%      Ns is number of targets in this snapshot
%      Nx is number of snapshots
%
% Author: Sean Holloway, John Paden, Theresa Stumpf
%
% See also: sim.crosstrack.m, sim.crosstrack_data.m, sim.crosstrack_example.m

%% Setup

% Load standard physical constants like c = speed of light
physical_constants;

% Nx: number of snapshots
Nx = size(target_model.z,2);
% Nc: number of channels
Nc = length(param.src.y_pc);
% dt: Fast time sample spacing (sec)
dt = 1/param.src.fs;
% Nt: number of range bins (fast time samples)
Nt = floor((param.src.t1-param.src.t0)/dt);

% time = fast time axis (sec)
time = param.src.t0 + dt*(0:Nt-1).';

% freq = fast time frequency axis (Hz)
df = 1/(Nt*dt);
freq = df*fftshift([-floor(Nt/2) : floor((Nt-1)/2)].');

% ref = reference function (eventually in frequency domain)
tref = (param.src.t0+param.src.t1)/2;
ref = param.src.ft_func(time - tref);
ref = ref./max(ref);
ref = fft(ref);

% From Sravya: SCALING WHEN USING FFT.
% Mohanad: Scaling ref in f-domain means that norm(ref./sqrt(N))^2=1, which is nice
% because it will do the job without affecting the total SNR of the data.
ref = ref/norm(ref);%sqrt(Nt);
if 0
  % Sravya's way: I don't think it matters if rcs is real or complex
  if isreal(target_model.rcs(1))
    ref = ref/Nt;
  else
    ref = ref/sqrt(Nt);
  end
end

% Transmit beamforming
Q = size(target_model.z,1);
if isfield(param.src,'tx_weights') && ~isempty(param.src.tx_weights)
  comp_tx_weight = param.src.tx_weights;
else
  comp_tx_weight = ones(Q,1);
end

% Array calibration errors
if isfield(param,'error_params') && ~isempty(param.error_params)
  error_params   = param.error_params;
  error_ypc      = error_params.error_ypc;
  error_zpc      = error_params.error_zpc;
  error_phase    = error_params.error_phase;
  error_g_s      = error_params.error_g_s;
  error_g_p      = error_params.error_g_p;
  error_g_offset = error_params.error_g_offset;
else
  error_ypc      = zeros(Nc,1);
  error_zpc      = zeros(Nc,1);
  error_phase    = zeros(Nc,1);
  error_g_s      = zeros(Nc,1);
  error_g_p      = zeros(Nc,1);
  error_g_offset = zeros(Nc,1);
end

if isfield(param,'array_gain_model_fh') && ~isempty(param.array_gain_model_fh)
 array_gain_model_fh =  param.array_gain_model_fh;
else
  array_gain_model_fh = @(x) (10.^(x/20));
end
  
% Mutual coupling matrix
if isfield(param.src,'mutual_coup_mtx') && ~isempty(param.src.mutual_coup_mtx)
  C = param.src.mutual_coup_mtx;
else
  C = eye(Nc);
end

% Multipath weights
if isfield(param,'MP_params')
  MP_params = param.MP_params;
end

% error_ypc = param.error_ypc;
% error_zpc = param.error_zpc;

%% Loop through each snapshot
sim_data = zeros(Nt, Nc, Nx);
td_min = inf;
td_max = -inf;

Num_targets = size(target_model.z,1);
for snapshot = 1:Nx
  %% Loop through each receiver
  for chan = 1:Nc
    %% Loop through each target
    for target = 1:Num_targets
      % Determine the range from the receiver to the target
      Rvec = [target_model.y(target) - (param.src.y_pc(chan)+error_ypc(chan));
        target_model.z(target,snapshot) - (param.src.z_pc(chan)+error_zpc(chan))];
      R = norm(Rvec,2);
      td = 2*R/c;
      if td < td_min
        td_min = td;
      end
      if td > td_max
        td_max = td;
      end
      
      % Incorporate gain and phase array errors. Location errors were
      % incorporated previously in Rvec. 
      % Multiplying theta by the sign of y guarantees that +/-theta goes with +/-y
      target_doa = sign(target_model.y(target)) * atan(abs(target_model.y(target))/abs(target_model.z(target,snapshot)));
      gain_error_exp = error_g_s(chan).*(sin(target_doa)-sin(error_g_p(chan))).^2 + error_g_offset(chan);
%       gain_error = array_gain_model_fh(gain_error_exp);
      gain_error  = 10.^(gain_error_exp./20);
      phase_error = error_phase(chan);
      pg_error = gain_error .* exp(1i*phase_error);
        
      amp = target_model.rcs(target,snapshot); 
      tx_weight = comp_tx_weight(target);
      
      % Add this target's energy to the simulated data matrix 
      if 1
        % Use this for comparisons of DOA estimation methods. The idea here
        % is that all targets have the same SNR for the purpose of
        % comparison.
        target_sim_data = amp *tx_weight*pg_error*exp(-1i*2*pi*param.src.fc*td) ...
          * exp(-1i*2*pi*freq*(td-tref)) .* ref;
        target_sim_data = sqrt((10^(param.snr_db/10))/2)*target_sim_data./(norm(target_sim_data));
        %       sim_data_norm(target) = norm(target_sim_data);
        sim_data(:,chan,snapshot) = sim_data(:,chan,snapshot) + target_sim_data;
      else
        % Use this for normal 2D simulations
        sim_data(:,chan,snapshot) = sim_data(:,chan,snapshot) ...
          + amp *tx_weight*pg_error*exp(-1i*2*pi*param.src.fc*td) ...
          * exp(-1i*2*pi*freq*(td-tref)) .* ref;
      end
      
      if isfield(param,'MP_params')
        % Incorporate MPCs associated with this target, if any
        y_pos_MP = MP_params{target}.y_pos_MP;
        z_pos_MP = MP_params{target}.z_pos_MP;
        w_MP     = MP_params{target}.w_MP;
        w_MP = (10.^(w_MP./20));
        
        for MP_idx = 1:length(w_MP)  
          MPC_signal = [];
          % Distance from the main target location to the reflctor
          R_MP_tr = sqrt((target_model.y(target)-y_pos_MP(MP_idx)).^2 + ...
            (target_model.z(target,snapshot)-z_pos_MP(MP_idx)).^2);
          
          % Distance from the reflector location to the radar
          R_MP_rr = sqrt((y_pos_MP(MP_idx) - (param.src.y_pc(chan)+error_ypc(chan))).^2 + ...
            (z_pos_MP(MP_idx) - (param.src.z_pc(chan)+error_zpc(chan))).^2);
          
          doa_MP = sign(y_pos_MP(MP_idx)) * atan(abs(y_pos_MP(MP_idx)/z_pos_MP(MP_idx)));
          gain_error_exp = error_g_s(chan).*(sin(doa_MP)-sin(error_g_p(chan))).^2 + error_g_offset(chan);
          gain_error = array_gain_model_fh(gain_error_exp);
%           gain_error  = exp(-gain_error_exp./2);
          phase_error = error_phase(chan);
          pg_error = gain_error .* exp(1i*phase_error);
          
          % MPC  signal: the signal arrives at the target's location after
          % traveling a distance R. Then this signal is delayed by the
          % travel time from the target location to the reflector location
          % then to the back to the radar receiver.
          % amp is the target's signal. Using different amp means different
          % target.
          MPC_signal = amp *tx_weight*w_MP(MP_idx)*pg_error* ...                               % Amplitude reduces by w_MP
            exp(-1i*2*pi*param.src.fc*R/c) * exp(-1i*2*pi*freq*(R/c-tref))* ...                % Radar to target location
            exp(-1i*2*pi*param.src.fc*R_MP_tr/c) .* exp(-1i*2*pi*freq*(R_MP_tr/c-tref))* ...   % Target to reflector
            exp(-1i*2*pi*param.src.fc*R_MP_rr/c) .* exp(-1i*2*pi*freq*(R_MP_rr/c-tref)).* ...  % Reflector to receiver
            ref;
          
          % Contaminated signal
          sim_data(:,chan,snapshot) = sim_data(:,chan,snapshot) + MPC_signal;
%         end
        end
      end

    end
  end
  
  % Incorporate mutual coupling matrix
  for binIdx = 1:Nt
    sim_data(binIdx,:,snapshot) = C * sim_data(binIdx,:,snapshot).';
  end
  
end

fprintf('    time gate ranges from %.0f ns to %.0f ns\n', time(1)*1e9, time(end)*1e9);
fprintf('    td ranges from %.0f ns to %.0f ns\n', td_min*1e9, td_max*1e9);

%% Add noise to sim_data
for chan = 1:Nc
  noise_data = 10.^(param.src.noise_power(chan)/20)*(randn(Nt,1,Nx) ...
    +1i*randn(Nt,1,Nx))/sqrt(2);
  sim_data(:,chan,:) = sim_data(:,chan,:) + noise_data;
end

%% Apply fast time frequency domain window to sim_data and convert to time domain
if 1 ||(isfield(param,'optimal_test') || isfield(param,'suboptimal_test'))
  % For MOE: Sravya used flipped Hanning window
  sim_data = ifft(sim_data .* repmat(ifftshift(param.src.ft_wind(Nt)),[1 Nc Nx]));
else
  % Other simulations: John used normal Hanning window
  sim_data = ifft(sim_data .* repmat(param.src.ft_wind(Nt),[1 Nc Nx]));
end

% figure;imagesc(10*log10(abs(squeeze(sim_data(:,8,:)))));
%% Debug:
if 0
  % Debug: Plot surface model. Here we plot the range-bins, but the targets are plotted in crosstrack.m
  h = figure(1000);
  hold on,circle(0,0,R_shell_values);
  
  ylim([-2500 -900])
  xlabel('Cross-track, y (m) ')
  ylabel('Range, z (m)')
  set(h,'position',[194   280   655   586])
  if param.dist_target_factor >1
    title('Distributed targets in constant range cylinders, Uniformly sampled  in wavenumber','interpreter','Latex')
  else
    title('Point targets in constant range cylinders, Uniformly sampled in wavenumber','interpreter','Latex')
  end
  grid on
end

end
