% script sim.doa_example_multipath
%
% Example setup scripts for the sim.doa function. This includes tutorial
% examples that illustrate how to simulate with multipath and mutual
% coupling.
%
% Author: John Paden, Theresa Stumpf
%
% See also: sim.doa.m

%% Multipath and Mutual coupling effects
% =========================================================================
physical_constants
param = [];

% -----------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------
% Source parameters
Nc                = 7;
param.src.fc      = 195e6;%195e9;
param.src.BW      = 10e6;%10e6;
param.src.fs      = param.src.BW;
param.src.f0      = param.src.fc-param.src.BW/2;
param.src.f1      = param.src.fc+param.src.BW/2;
param.src.ft_wind = @(N) hanning(N);
param.src.ft_wind  = @boxcar;

% DOA method parameters
param.method.list                   = [7];
param.method.Nsv                    = 32;
param.method.OneD_Nsv               = 128;
param.method.theta_guard            = 1.5/180*pi;
param.method.nb_nd.init             = 'ap';
param.method.wb_td.init             = 'ap';
param.method.wb_td.widening_factor  = 1;
param.method.wb_fd.init             = 'ap';
param.method.wb_fd.filter_banks     = 1;

% -----------------------------------------------------------------------
% Beam steering and 3dB beamwidth
% -----------------------------------------------------------------------
% theta contains the DoAs from a single range-bin. The goal here is to
% test the effect of multipath on data structure. So, 1 range-bin is
% enough.
lambda = c/param.src.fc;
d_y    = lambda/2;

phase_center      =  zeros(3,Nc);
phase_center(2,:) = d_y/2*(0:Nc-1);

y_pc = phase_center(2,:).';
z_pc = phase_center(3,:).';

param.src.y_pc       = y_pc;
param.src.z_pc       = z_pc;
k = 2*pi/(lambda/2);

% Transmit beamforming
if 0
  % Steering martix function handle
  %   A = @(theta) sim.steering_mtx(theta,param);
  A = @(theta) steering_mtx(theta,param);
  
  % Test if the SVs are really orthogonal if orthogonal DoAs were used
  if 0
    SVs = A(DOA_orthogonal);
    for idx = 1:size(SVs,2)
      orth_check(idx) = SVs(:,1)'*SVs(:,idx);
    end
    orth_check = orth_check
  end
  
  % Beam steering angle
  steer_ang = 0*pi/180;
  
  % Complex transmit weights
  tx_weights     = hanning(Nc);
  %   tx_weights     = [hanning(4).',0 0 0 0].';
  steering_delay = (y_pc*sin(steer_ang) - z_pc*cos(steer_ang))/c;
  w = tx_weights.*exp(-1i*4*pi*param.src.fc*steering_delay);
  
  % Set the limits of the DoAs to within 3dB of the antenna beampattern
  theta_RP = linspace(-90,90,2048)*pi/180;
  RP = abs(w.'*A(theta_RP)).^2;
  RP = RP./max(RP);
  RP_dB = 10*log10(RP);
  
  idxs_3dB  = find(RP_dB>=-3.5);
  RP_3dB    = RP_dB(idxs_3dB);
  theta_3dB = theta_RP(idxs_3dB);
  
  %   beam_doa_lims = [min(theta_3dB) max(theta_3dB)];
end

% -----------------------------------------------------------------------
% Multipath effect
% -----------------------------------------------------------------------
% Set the multipath componens, MPCs, associated with each actual DoA. Number of
% multipath components can be different for each actual DoA. Assume same
% number for now.

if 1
  % Orthogonal DoAs (to test that MP is working fine)
  % Also use it to set angles and delays rather than locations of targets
  L = Nc*d_y;
  k_spacing = [-(Nc-1)/2:1:(Nc-1)/2]*(lambda/L);
  %     k_spacing = [-Nc/2:1:Nc/2-1]*(lambda/L); % John's
  DOA_orthogonal = asin(k_spacing);
  %     theta = DOA_orthogonal(4);
  %     theta_MP = DOA_orthogonal([3 5]);
  theta = [10]*pi/180;
  theta_MP = [-50 40]*pi/180;
  z_pos = -500;    % z-axis of the main target
  Range = abs(z_pos)/cos(theta); % Range to the main target
  y_pos = Range*sin(theta); % y-axis of the main target
  
  z_pos_MP = [z_pos+400, z_pos+480];
  y_pos_MP = z_pos_MP.*tan(theta_MP);
  
  % MPCs delay relative to the BW (or relative to the main path delay)
  %     tau_MP = 1/param.src.BW*ones(1,2);
  target_to_reflector_dist   = sqrt((z_pos-z_pos_MP).^2 + (y_pos_MP-y_pos).^2);
  reflector_to_receiver_dist = sqrt(z_pos_MP.^2 + y_pos_MP.^2);
  total_distance_MP = target_to_reflector_dist + reflector_to_receiver_dist;
  tau_MP = (total_distance_MP - Range)./c; % Relative MPCs delay
else
  theta = [+5]*pi/180;
  
  % Main path location and DoA
  z_pos = -500;    % z-axis of the main target
  Range = abs(z_pos)/cos(theta); % Range to the main target
  y_pos = Range*sin(theta); % y-axis of the main target
  
  % MPCs location, delay, and DoA
  z_pos_MP = [z_pos+400, z_pos+480]; % z-axis of the reflectors
  y_pos_MP = [y_pos+200, y_pos+200]; % y-axis of the reflectors
  Range_MP = sqrt(z_pos_MP.^2 + y_pos_MP.^2); % Range to the reflectors
  
  target_to_reflector_dist   = sqrt(z_pos_MP.^2 + abs((y_pos_MP-y_pos)).^2);
  reflector_to_receiver_dist = sqrt((z_pos-z_pos_MP).^2 + y_pos_MP.^2);
  total_distance_MP = target_to_reflector_dist + reflector_to_receiver_dist;
  tau_MP = (total_distance_MP - Range)./c; % Relative MPCs delay
  
  theta_MP = sign(y_pos_MP).*atan(y_pos_MP./(abs(z_pos_MP))); % MPCs DoA relative to the nadir
end

TBP = tau_MP/(1/param.src.BW);
fract_BW = param.src.BW/param.src.fc;

if 1
  % Debug: plot problem geometry
  figure(100);clf
  hold on
  h1 = plot(y_pos,z_pos,'*b','LineWidth',10);  % Main target
  h2 = plot(y_pos_MP(1),z_pos_MP(1),'*r','LineWidth',10); % Reflector 1
  h3 = plot(y_pos_MP(2),z_pos_MP(2),'*m','LineWidth',10); % Reflector 2
  h4 = plot(0,0,'sk','LineWidth',10); % Radar
  
  plot([y_pos,y_pos_MP(1)],[z_pos,z_pos_MP(1)],'r') % Target to reflector 1
  plot([y_pos,y_pos_MP(2)],[z_pos,z_pos_MP(2)],'m') % Target to reflector 2
  plot([y_pos 0],[z_pos 0],'b')                     % Radar to target
  plot([0 y_pos_MP(1)],[0 z_pos_MP(1)],'r')         % Radar to reflector 1
  plot([0 y_pos_MP(2)],[0 z_pos_MP(2)],'m')         % Radar to reflector 2
  plot([0 0],[0 z_pos],'-.k')                       % nadir
  
  grid on
  xlabel('y-axis')
  ylabel('z-axis')
  title('Multipath problem geometry')
  legend([h1 h2 h3 h4],'Target','Reflector 1','Reflector 2','Radar','Location','southwest')
end
% -----------------------------------------------------------------------
% Mutual coupling effect
% -----------------------------------------------------------------------

% Define S-parameters matrix (symmetric and Toepletz) assuming matched
% dipoles
if 0
  S = [0.000   0.200  -0.40  -0.2j   0.10   0.05j   0.010;...
    0.200   0.000   0.20  -0.40  -0.2j   0.100   0.05j;...
    -0.400   0.200   0.00   0.20  -0.40  -0.20j   0.100;...
    -0.20j  -0.400   0.20   0.00   0.20   0.400  -0.20j;...
    0.100  -0.20j  -0.40   0.20   0.00   0.200  -0.400;...
    0.05j   0.100   0.2j  -0.40   0.20   0.000   0.200;...
    0.010   0.05j   0.10  -0.2j  -0.40   0.200   0.000];
  % Define Mutual coupling matrix
  C = eye(Nc) - S;
  %   C = eye(Nc);
  param.src.mutual_coup_mtx = C;
else
  param.src.mutual_coup_mtx = [];
end

% -----------------------------------------------------------------------
% Run the simulator
% -----------------------------------------------------------------------

% DOA monte carlo setup
SNR   = 70;
Nsnap = 1*101;
param.monte.runs  = 1;
param.monte.random_seed_offset = 0;

param.monte.DOA = theta*180/pi;
param.monte.SNR = SNR*ones(length(theta),1);
param.monte.Nsnap = Nsnap*ones(length(theta),1);
comp_tx_weight = ones(length(theta)+length(theta_MP),1);

param.monte.DOA_MP = theta_MP*180/pi;
delay_MP = [0 tau_MP];

% Setting the MPC weight to -inf ignores the MPC (controls number of
% MPCs). In dB.
w_MP     = [0, -10 -10];

% Transmit beamforming weights
param.src.tx_weights = comp_tx_weight;
% MPCs relative delays
param.src.tau_MP = delay_MP;
% MPCs weights (in dB)
param.monte.w_MP = w_MP;
% Set MCPs DoA limits
param.method.src_limits = [];
for MPC_idx = 1:length(param.monte.DOA)+length(param.monte.DOA_MP)
  param.method.src_limits{MPC_idx} = [-90 +90]*pi/180;
  %       param.method.src_limits{idx} = [beam_doa_lims(1) beam_doa_lims(2)];
end

if 1
  clear doa
  results = sim.doa(param);
  est_doa = squeeze(results.theta_est{param.method.list}(1,1,:));
end

%     param.src.Nsnap = param.monte.Nsnap(1);
param.src.SNR   = param.monte.SNR;
param.src.Nsnap = param.monte.Nsnap;
param.src.DOAs  = param.monte.DOA;
param.src.DOAs_MP  = param.monte.DOA_MP;

[Data,DCM] = sim.doa_wideband_data(param);

eigval = real(sort(eig(DCM),'descend')).';

%     sigma_n = mean(eigval(4:end));
%     sigma_s = mean(abs(Data(:)).^2);
%     SNR = 10*log10(sigma_s/sigma_n - 1);

fprintf('\n\n TBP = [%2.4f %2.4f]\n \nFractional BW = %2.4f\n \nDoA = %2.2f\n \nMPCs DOA = [%2.2f %2.2f]\n\n',...
  TBP(1),TBP(2),fract_BW,theta*180/pi,theta_MP(1)*180/pi,theta_MP(2)*180/pi)

if 0
  % Plot true vs estimated DoA
  figure(101);clf;
  doa_true = param.src.DOAs.';
  doa_hat = est_doa*180/pi;
  scatter(doa_true,doa_hat,20,'fill');%,'MarkerFaceColor','b')
  
  xlim([-60 +60])
  ylim([-60 +60])
  xlabel('True DoA (deg.)')
  ylabel('Estimated DoA (deg.)')
  title('1D simulator: multipath')
  grid on
end

% Plot DCM eigen values in dB
figure(1);clf
stem(10*log10(abs(eigval./max(eigval))),'filled','b','LineWidth',1.5)
xlabel('Eigenvalue index')
ylabel('Eigenvalue (dB)')
title('Normalized eigenvalues of R')
grid on

% Plot the magnitude and phase of the DCM
figure(2);clf
subplot(211)
imagesc(10*log10(abs(DCM)./max(abs(DCM(:)))))
xlabel('Sensor index')
ylabel('Sensor index')
h = colorbar;
ylabel(h,'Normalized |R| (dB)')
colormap parula

subplot(212)
imagesc(angle(DCM)*180/pi)
xlabel('Sensor index')
ylabel('Sensor index')
h = colorbar;
ylabel(h,'Angle of R (deg.)')
colormap parula

% -----------------------------------------------------------------------
% Plot true vs estimated DoA
% -----------------------------------------------------------------------
%   figure(101);clf;
%   hold on
%   for doa_idx = 1:length(actual_doa)
%     doa_true = theta_multipath{doa_idx}.'*180/pi;
%     doa_hat = est_doa{doa_idx}(:)*180/pi;
%     scatter(doa_true,doa_hat,20,'fill');%,'MarkerFaceColor','b')
%   end
%
%   xlim([-60 +60])
%   ylim([-60 +60])
%   xlabel('True DoA (deg.)')
%   ylabel('Estimated DoA (deg.)')
%   title('1D simulator: multipath')
%    grid on
