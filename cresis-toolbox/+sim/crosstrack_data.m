function [sim_data, sources_true_all] = crosstrack_data(param,target_model)
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

%%%% Y DIMENSION DESCRIPTION CHECK Ns by 1


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
ref = fft(ref);


%%%%%%%%%%%%%plotting fft (impulse fft uniform over all f)
% L = length(ref);
% P2 = abs(ref/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = fs*(-L/4:(L/4))/L;
% figure(1),clf;
% plot(f,P1) 
% title('Single-Sided Spectrum of ref')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% freq1 = fs*(-L/2:(L/2))/L;
% figure(2),clf;
% plot(freq1(1:end-1),P2) 
% title('two-Sided  Spectrum of ref')
% xlabel('f (Hz)')
% ylabel('|P2(f)|')

%%%%%%%%%%%%%%%%%%%%%%%%




sim_data = zeros(Nt, Nc, Nx);
td_min = inf;
td_max = -inf;
%% Loop through each snapshot
for snapshot = 1:Nx
  %% Loop through each receiver
  for chan = 1:Nc
    %% Loop through each target
    for target = 1:size(target_model.z,1)
      % Determine the range from the receiver to the target
      Rvec = [target_model.y(target) - param.src.y_pc(chan);
        target_model.z(target,snapshot) - param.src.z_pc(chan)];
      R = norm(Rvec,2);
      td = 2*R/c;
      if td < td_min
        td_min = td;
      end
      if td > td_max
        td_max = td;
      end
      amp = target_model.rcs(target,snapshot);
      % Add this target's energy to the simulated data matrix
      sim_data(:,chan,snapshot) = sim_data(:,chan,snapshot) ...
        + amp * exp(-1i*2*pi*param.src.fc*td) ...
        * exp(-1i*2*pi*freq*(td-tref)) .* ref;
% %%%SRAVYA    
% mag(:,target) = abs(ifft(sim_data(:,1,1)));
% figure(1),clf;
% plot(mag(:,target))


    end
     
    
  
  mag = abs(ifft(sim_data(:,1,chan)));
  figure(10);
  subplot(4,2,chan)
  plot(mag)
  
    
    
  end
  
  
end
fprintf('    time gate ranges from %.0f ns to %.0f ns\n', time(1)*1e9, time(end)*1e9);
fprintf('    td ranges from %.0f ns to %.0f ns\n', td_min*1e9, td_max*1e9);

%% Add noise to sim_data
for chan = 1:Nc
    
    
%%% param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
  
 noise_data = 10.^(param.src.noise_power(chan)/20)*(randn(Nt,1,Nx) ...
   +1i*randn(Nt,1,Nx))/sqrt(2);



 % I THINK /20 IS WRONG INSTEAD WE SHOULD DO /10 above
 
%   noise_data = 10.^(param.src.noise_power(chan)/10)*(randn(Nt,1,Nx) ...
%     +1i*randn(Nt,1,Nx))/sqrt(2);


%   noise_data = sqrt(param.src.noise_power(chan))*(randn(Nt,1,Nx) ...
%     +1i*randn(Nt,1,Nx))/sqrt(2);


%   noise_data = sqrt(1/2)*(randn(Nt,1,Nx) ...
%     +1i*randn(Nt,1,Nx));


% 
%   noise_data = sqrt(Nx/2).*(randn(Nt,1,Nx) ...
%     +1i*randn(Nt,1,Nx));

  sim_data(:,chan,:) = sim_data(:,chan,:) + noise_data;
  noise(:,chan,:) =  noise_data;
  
end

%% Apply fast time frequency domain window to sim_data and convert to time domain
sim_data = ifft(sim_data .* repmat(param.src.ft_wind(Nt),[1 Nc Nx]));





%% SRAVYA  (TRUE SOURCES)
if 1
for snapshot =1:Nx;

R_bins_values = [0;time*c/2];

% 
% fig_idx = 1
% figure(fig_idx),clf;


% for chan = 1:8
     chan =1;
  
    
    clear Rvec 
    
    for target = 1:size(target_model.z,1)
        % Determine the range from the receiver to the target
        Rvec = [target_model.y(target) - param.src.y_pc(chan);
            target_model.z(target,snapshot) - param.src.z_pc(chan)];
        R_targets(target) = norm(Rvec,2);
        
    end
    
    clear Rvec
    
    target_model.R_targets = R_targets;
    for dR_idx = 1:Nt
        
        bin_sort_idxs =  find ( target_model.R_targets <= R_bins_values(dR_idx +1) & target_model.R_targets > R_bins_values(dR_idx));
        target_model.R_targets(bin_sort_idxs) = dR_idx-0.5;
        sources_true(dR_idx) = length(bin_sort_idxs);
        
    end
    
    doa_target = atan(target_model.y./target_model.z(:,snapshot))*180/pi;
    
    clear bin_sort_idxs 
    
    
% %     figure(fig_idx)
% %     plot(-1*doa_target,target_model.R_targets,'*');grid on
% %     
% %     %%% CHECK doa - and +
% %     set(gca,'Ydir','reverse')
% %     
% %     set(gca,'ytick',1:100)
% %     ylim([0 100])
% %     
    
    
    %
    % figure(6);plot(target_model.y,target_model.R_targets,'*');grid on
    % set(gca,'Ydir','reverse')
    %
    % set(gca,'ytick',1:100)
    % ylim([0 100])
    
%     fig_idx = fig_idx+1
 
sources_true_all(:,snapshot) = sources_true.';

clear sources_true
end

end








% end

%%
end

