% script sim.crosstrack_example_multipath
%
% Example multipath and mutual coupling simulation. Calls crosstrack.m and
% crosstrack_rmse.m.
%
% Author: Mohanad Al-Ibadi, Sravya Athinarapu, Sean Holloway, John Paden,
% Theresa Stumpf
%
% See also: sim.crosstrack.m, sim.crosstrack_data.m,
% sim.crosstrack_example_*.m, sim.crosstrack_rmse.m

physical_constants;

%% Multipath and Mutual coupling effects
% ----------------------------------------------------------------------
% Here is how we simulate the effect of the MPCs created by a target
% at some range-bin on a target at another range-bin:
% 1) Set the location of the main target that is going to create the MPCs.
% 2) Create the reflectors locations by defining DOA_MP and z_MP then
% calculating y_MP.
% The MPC may not affect the same range-bin that created this MPC, unless
% the MPC arrives at the same time (within the range-resolution) as the
% main path.

%% Setup simulation parameters
param = [];

% Debug level of 3 causes this function to stop early and output
% simulated data, simulation parameters, and array processing results
% in a "results" structure.
param.debug_level = 3;

%% Source parameters
fc = 195e6;
BW = 30e6;
% Nc: number of sensors
Nc = 7;
flight_height  = 0;
% Set the DoAs and weights of the main paths and MPCs
% Orthogonal DoAs (to test if it is actually working)
lambda = c/fc;
d_y    = lambda/2;
L = Nc*d_y;
k_spacing = [-(Nc-1)/2:1:(Nc-1)/2]*(lambda/L);
DOA_orthogonal = asin(k_spacing);
%     DOA = DOA_orthogonal(5).';
%     DOA_MP = DOA_orthogonal([1 7]).';
DOA = [30]*pi/180;
DOA_MP =[-5 8].'*pi/180;

% Setting the MPC weight to -inf ignores the MPC (controls number of
% MPCs).
MP_weights = [-inf -inf]; % in dB

% Main target location
z_pos = -1500;
Range = abs(z_pos)/cos(DOA); % Range to the main target
y_pos = Range*sin(DOA);      % y-axis of the main target

% MPCs location
z_pos_MP = [z_pos+400, z_pos+480];
y_pos_MP = z_pos_MP.*tan(DOA_MP.');
Range_MP = abs(z_pos_MP)./cos(DOA_MP.');

%     Range_MP = Range*ones(1,length(DOA_MP));
%     z_pos_MP = -Range_MP.*cos(DOA_MP.');
%     y_pos_MP = Range_MP.*sin(DOA_MP.');
%     R_main_target = Range*ones(1,length(DOA_MP));

target_to_reflector_dist   = sqrt((z_pos-z_pos_MP).^2 + (y_pos_MP-y_pos).^2);
reflector_to_receiver_dist = sqrt(z_pos_MP.^2 + y_pos_MP.^2);
total_distance_MP = Range + target_to_reflector_dist + reflector_to_receiver_dist;
tau_MP = (total_distance_MP - 2*Range)./c; % Relative MPCs delay

% Place a target at the same range as the MPCs (this target doesn't
% creat a MPC). Make sure that this target and the one created the MPCs
% are not within one range-resolution distance from each other.
%     Range_2 = total_distance_MP(1)/2;
%     y_pos_2 = sqrt(Range_2^2 - z_pos^2);
%     DOA_2 = sign(y_pos_2)*atan(abs(y_pos_2/z_pos));

%     DOA = [DOA DOA_2].';
DOA = DOA;

if 0
  % Debug: plot problem geometry
  figure(100);clf
  hold on
  h1 = plot(y_pos,z_pos,'*b','LineWidth',10);             % Main target
  h2 = plot(y_pos_MP(1),z_pos_MP(1),'*r','LineWidth',10); % Reflector 1
  h3 = plot(y_pos_MP(2),z_pos_MP(2),'*m','LineWidth',10); % Reflector 2
  h4 = plot(0,0,'sk','LineWidth',10);                     % Radar
  
  plot([y_pos,y_pos_MP(1)],[z_pos,z_pos_MP(1)],'r','LineWidth',1.5) % Target to reflector 1
  plot([y_pos,y_pos_MP(2)],[z_pos,z_pos_MP(2)],'m','LineWidth',1.5) % Target to reflector 2
  plot([y_pos 0],[z_pos 0],'b','LineWidth',1.5)                     % Radar to target
  plot([0 y_pos_MP(1)],[0 z_pos_MP(1)],'r','LineWidth',1.5)         % Radar to reflector 1
  plot([0 y_pos_MP(2)],[0 z_pos_MP(2)],'m','LineWidth',1.5)         % Radar to reflector 2
  plot([0 0],[0 z_pos],'-.k')                                       % nadir
  
  grid on
  xlabel('y-axis')
  ylabel('z-axis')
  title('Multipath problem geometry')
  legend([h1 h2 h3 h4],'Main target','Reflector 1','Reflector 2','Radar','Location','southwest')
end

% MP parameters associated with each target
MP_params{1}.y_pos_MP = y_pos_MP;
MP_params{1}.z_pos_MP = z_pos_MP;
MP_params{1}.w_MP     = MP_weights;

MP_params{2}.y_pos_MP = y_pos_MP;
MP_params{2}.z_pos_MP = z_pos_MP;
MP_params{2}.w_MP     = [-inf -inf]; % Second target doesn't create MPCs

param.MP_params = MP_params;

% Set the range gate based on the flight height and targets DoA
max_range      = max(Range + target_to_reflector_dist + reflector_to_receiver_dist);
min_range      = min(Range);

param.src.f0                      = fc-BW/2;
param.src.f1                      = fc+BW/2;
param.src.t0                      = 2*(min_range-500)/c;
param.src.t1                      = 2*(max_range+1000)/c;
if param.src.t0 < 0
  param.src.t0 = 0;
end
param.src.ft_func                 = @(t) tukeywin_cont(t * BW);
param.src.ft_wind                 = @(N) hanning(N);
param.src.lever_arm.fh            = @sim.lever_arm_example;

% Nsep: number of lambda/4 steps between sensors
Nsep = 1;
% Arguments for a linear array in y-dimension
% param.src.lever_arm.fh_args: Arguments for a linear array in y-dimension
param.src.lever_arm.fh_args       = {[], 1, 1:Nc, [0; Nsep*c/fc/2; 0]};
phase_center = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});

% Change reference sensor if you want to
ref_chan = 1;
phase_center(2,:) = phase_center(2,:) - phase_center(2,ref_chan);
phase_center(3,:) = phase_center(3,:) - phase_center(3,ref_chan);

% Account for roll angle. Default is 0
theta_roll = 0 * pi/180;
phase_center(2,:) = phase_center(2,:) * cos(theta_roll);
for chan_idx = 1:Nc
  d_chan_ref      = sqrt(abs(diff(phase_center(2,[ref_chan  chan_idx])))^2 + abs(diff(phase_center(3,[ref_chan  chan_idx])))^2);
  phase_center(3,chan_idx) = phase_center(3,chan_idx) + d_chan_ref * sin(theta_roll);
end
param.src.phase_center = phase_center;
ypc = -phase_center(2,:).';
zpc = -phase_center(3,:).';

dt = 1/BW;
% Nt: number of range bins (fast time samples)
Nt = floor((param.src.t1-param.src.t0)/dt);

% Time = fast time axis (sec)
Time = param.src.t0 + dt*(0:Nt-1).';
%   param.src.noise_power             = 10*log10(BoltzmannConst*290*abs(param.src.f0-param.src.f1)) + 2*ones(1,Nc);
param.src.noise_power             = zeros(1,Nc);

% DOA method parameters
param.method.list                   = [7];

%% Simulation Runs Setup

% Cross track monte carlo setup
param.monte.target_func = @sim.surface_gen;
param.monte.runs = 1;
param.monte.random_seed_offset = 0;

% Target surface parameters
surf_param = [];
surf_param.z.mean = z_pos;
surf_param.z.rms_height = 0;
surf_param.z.corr_length_x = 400;
surf_param.z.corr_length_y = 400;
surf_param.rcs_in.mean = 0;
surf_param.rcs_in.var = 1e2;
surf_param.dy = 10;
%   surf_param.y_range = [-2500 2500];
surf_param.dx = 10;
surf_param.x_range = [-2500 2500];
surf_param.x = [-1500:-1500+100]; % 199...100 range-lines
%   surf_param.y = [-1500:20:1500].';
%   surf_param.y = [y_pos y_pos_2].';
surf_param.y = y_pos;

% y_range should >= maximum y (i.e. targets shoulb be inside the imaged
% swath)
y_range = max(target_to_reflector_dist + reflector_to_receiver_dist+Range);
if y_range == 0
  surf_param.y_range = [-1000 1000];
else
  surf_param.y_range = [-1.5*y_range 1.5*y_range];
end
param.monte.target_param{1} = surf_param;


% surf_param = [];
% surf_param.z.mean = 0.7;
% surf_param.z.rms_height = 0.5;
% surf_param.z.corr_length_x = 1000;
% surf_param.z.corr_length_y = 40;
% surf_param.rcs_in.mean = 0;
% surf_param.rcs_in.var = 100;
% surf_param.dy = 1;
% surf_param.y_range = [-1500 1500];
% surf_param.dx = 1;
% surf_param.x_range = [-1500 1500];
% surf_param.x = [-1:0.1:1];
% surf_param.y = [-50:1:50].';
% param.monte.target_param{2} = surf_param;

%% Array Processing parameters
array_param = [];

fc = (param.src.f0+param.src.f1)/2;
fs = param.src.f1-param.src.f0;

%% NN: total length of the sensor array (in lambda/4 units)
NN = Nc*Nsep;
% lambda: wavelength
lambda = c/fc;
% k: wavenumber
k = 2*pi/(lambda/2);
% My: over sampling factor
My = 4;
% dy: phase center spacing
dy = Nsep*lambda/4;
% dky and ky: y-component of wavenumber (spacing and axis)
dky = 2*pi / (Nc*dy) / My;
ky = dky * ifftshift(-floor(My*Nc/2) : floor((My*Nc-1)/2));
% theta: theta values associated with ky axis
theta = fftshift(asin(ky/k));
theta = [theta,pi/2];
array_param.Nsv = {'theta', asin(ky/k)};

array_param.sv_fh = @array_proc_sv;

array_param.dbin = 1;
array_param.dline = 1;

array_param.bin_rng = 0;
array_param.line_rng = -50:50;

array_param.Nsrc = 2;
Nsrc = array_param.Nsrc;

array_param.init = 'ap';
array_param.doa_theta_guard = 2*pi/180; %(max(theta)-min(theta))/(4*Nc);

array_param.Nsubband = 1;
dt = 1/fs;
array_param.imp_resp.time_vec = -3*array_param.Nsubband*dt : dt/8 : 3*array_param.Nsubband*dt;
BW = abs(param.src.f1 - param.src.f0);
array_param.imp_resp.vals = tukeywin_cont(array_param.imp_resp.time_vec / BW);

for idx = 1:array_param.Nsrc
  array_param.doa_constraints(idx).method = 'fixed';
  array_param.doa_constraints(idx).init_src_limits = [min(theta) max(theta)]*180/pi;
  array_param.doa_constraints(idx).src_limits = [min(theta) max(theta)]*180/pi;
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

param.array_param = array_param;
clear array_param;

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

sim_data = squeeze(results.sim_data{1}); % Nt-by-Nx-by-Nc
%   est_doa         = results.tomo.doa;
actual_doa_cell = results.actual_doa;
array_param     = results.array_param;

% Convert actual_doa from cell array into a matrix (maximum Nsrc targets
% per range-bin).
actual_doa = NaN(length(array_param.bins),Nsrc,length(array_param.lines));
for lineIdx = 1:length(array_param.lines)
  for binIdx = 1:length(array_param.bins)
    if ~isempty(actual_doa_cell{binIdx,lineIdx})
      doa_tmp = actual_doa_cell{binIdx,lineIdx};
      if length(doa_tmp)<Nsrc
        doa = NaN(Nsrc,1);
        doa(1:length(doa_tmp)) = doa_tmp;
      else
        doa = doa_tmp(1:Nsrc);
      end
      actual_doa(binIdx,1:Nsrc,lineIdx) = actual_doa_cell{binIdx,lineIdx}; % Nt*Nsrc*Nx
    end
  end
end

if 1
  % Calculate the range-bins to which the MPCs blong to. The range is
  % calculated as half of the total distance from the send to receive.
  range = Time*c/2;
  dr = c/(2*BW);
  MPC_rbins = [];
  for lineIdx = 1:length(array_param.lines)
    for MP_idx = 1:length(DOA_MP)
      %       rbins = find(abs(Range_MP(MP_idx)-range) < dr/2*2);
      [~,min_idx] = min(abs(total_distance_MP(MP_idx)/2 - range));
      MPC_rbins(MP_idx,lineIdx) = min_idx;
    end
  end
end

% Calculate the eigenvalues of each test range-bin
eigenvalue_all  = [];
sample_data_all = [];
sample_data     = [];
for lineIdx_idx = 1:length(array_param.lines)
  line_idx = array_param.lines(lineIdx_idx);
  good_rbins = sort(find(~isnan(actual_doa(:,2,lineIdx_idx))),'ascend');
  %     first_rbin = MPC_rbins(1);
  rqd_rbin = good_rbins(1);%min(good_rbins);
  range_bin_idxs = rqd_rbin + [0];
  for binIdx_idx = 1:length(range_bin_idxs)
    bin_idx = range_bin_idxs(binIdx_idx);
    sample_data = squeeze(sim_data(bin_idx,line_idx+array_param.line_rng,:)).';
    
    if 0
      % FB averaging and spatial smoothing
      DCM = 1/size(sample_data,2) * sample_data*sample_data';
      Rxx1 = 1/(size(sample_data(1:end-1,:),2))*sample_data(1:end-1,:) * sample_data(1:end-1,:)';
      Rxx2 = 1/(size(sample_data(2:end,:),2))*sample_data(2:end,:) * sample_data(2:end,:)';
      
      % Apply FB averaging for each subarray
      reflect_mtx = flipud(eye(size(Rxx1)));
      Rxx1_fb = (1/2)*(Rxx1 + reflect_mtx * conj(Rxx1) * reflect_mtx);
      Rxx2_fb = (1/2)*(Rxx2 + reflect_mtx * conj(Rxx2) * reflect_mtx);
      
      % Average the two DCMs
      Rxx_ss = 1/2*(Rxx1_fb+Rxx2_fb); % (Nc-1)-by-(Nc-1)
      
      % Handle the lost sensor (or DOF) such that the final DCM is Nc-by-Nc
      Rxx_tmp = zeros(Nc);
      Rxx_tmp(1:end-1,1:end-1) = Rxx_ss;
      Rxx_tmp(end,:) = DCM(end,:);
      Rxx_tmp(:,end) = DCM(:,end);
      
      DCM = Rxx_tmp; % Nc-by-Nc matrix
    elseif 0
      % FB averaging only
      DCM = 1/size(sample_data,2) * sample_data*sample_data';
      reflect_mtx = flipud(eye(size(Rxx_calib)));
      DCM = (1/2)*(DCM + reflect_mtx * conj(DCM) * reflect_mtx);
    else
      DCM = 1/size(sample_data,2) * sample_data*sample_data';
    end
    if bin_idx == rqd_rbin
      DCM_tmp = DCM;
    end
    
    eigenvalue = real(sort(eig(DCM),'descend'));
    eigenvalue_all(:,end+1) = eigenvalue;
    
    sample_data_all  = [sample_data_all;sample_data(:)];
  end
end
mean_eigval = mean(real(eigenvalue_all),2);
mean_eigval = mean_eigval./max(mean_eigval);

if 1
  % Debug
  % Estimate the SNR
  SNR = 10*log10(mean(abs(sample_data_all).^2));
  
  TBP = tau_MP/(1/BW);
  fract_BW = BW/fc;
  
  %     good_rbins_1 = find(~isnan(results.tomo.doa(:,1,lineIdx_idx)));
  %     good_rbins_2 = find(~isnan(results.tomo.doa(:,2,lineIdx_idx)));
  %     est_doa_1 = results.tomo.doa(good_rbins_1,1,lineIdx_idx)*180/pi;
  %     est_doa_2 = results.tomo.doa(good_rbins_2,2,lineIdx_idx)*180/pi;
  %     est_doa = sort(unique([est_doa_1' est_doa_2']),'ascend');
  
  if length(DOA_MP) == 2
    sprintf('\nTBP = [%2.4f  %2.4f]\n \nFractional BW = %2.4e\n \nSNR = %2.1f\n \n DoA = %2.2f deg.\n \n MPCs DoA = [%2.2f %2.2f] deg.\n',...
      TBP(1), TBP(2), fract_BW, SNR, DOA(1)*180/pi,DOA_MP(1)*180/pi,DOA_MP(2)*180/pi)
    
  elseif length(DOA_MP) == 1
    sprintf('\nTBP = %2.4f\n \nFractional BW = %2.4e\n \nSNR = %2.1f\n \n DoA = %2.2f\n \n MPCs DoA = %2.2f\n',...
      TBP(1), fract_BW, SNR, DOA*180/pi,DOA_MP*180/pi)
  end
  
  range_bins = good_rbins
  MPC_rbins = MPC_rbins(:,1)
end

if 0
  % Phase difference between sensors that we expect to see in the phase
  % of the DCM. It is calculated from geometry.
  %     ypc = -phase_center(2,:).';
  %     zpc = -phase_center(3,:).';
  clear max_expected_phase
  for chan_idx = 1:Nc
    % Distance from the a given sensor to the reference sensor
    L_chan      = sqrt(abs(diff(ypc([ref_chan  chan_idx])))^2 + abs(diff(zpc([ref_chan  chan_idx])))^2);
    % Array length in the direction of DOA
    L_doa   = L_chan*cos(DOA);
    % Extra delay with respect to the reference sensor
    L_extra = L_chan*sin(DOA);
    max_expected_phase(chan_idx,1) = 4*pi/lambda * L_extra * 180/pi;
  end
  max_expected_phase
  
  DCM_phase_UTM = zeros(Nc);
  for chan_idx = 1:Nc
    DCM_phase_UTM(chan_idx,chan_idx:end) = max_expected_phase(1:Nc-chan_idx+1);
  end
  DCM_phase_LTM = -triu(DCM_phase_UTM).';
  expected_DCM_phase = DCM_phase_UTM + DCM_phase_LTM;
end
%% Plots
% Plot the magnitude and phase of the DCM
% ---------------------------------------
figure(3000);clf
suptitle(sprintf('2D sim: DOA = %2.2f deg.\n',DOA(1)*180/pi))
subplot(211)
imagesc(10*log10(abs(DCM_tmp)./max(abs(DCM_tmp(:)))))
%   xlabel('Sensor index')
ylabel('Sensor index')
h = colorbar;
ylabel(h,'Normalized |R| (dB)')
colormap parula

subplot(212)
DCM_phase = angle(DCM_tmp);
%   DCM_phase = unwrap(DCM_phase,[],2);
%   DCM_phase = DCM_phase - diag(diag(DCM_phase));
imagesc(DCM_phase*180/pi);
%   imagesc(DCM_phase*180/pi)
%   xlabel('Sensor index')
ylabel('Sensor index')
h = colorbar;
ylabel(h,'\angle R (\circ)')
colormap parula
h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
if ~isempty(Ticks)
  h.Ticks = Ticks;
end

% Plot the phase of the DCM derived from array geometry
% -----------------------------------------------------
if 0
  figure(3001);clf
  suptitle(sprintf('From geometry: DOA = %2.2f deg.\n',DOA(1)*180/pi))
  subplot(211)
  imagesc(expected_DCM_phase)
  %   xlabel('Sensor index')
  ylabel('Sensor index')
  h = colorbar;
  ylabel(h,'\angle R (\circ) - Before unwrapping')
  colormap parula
  h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  if ~isempty(Ticks)
    h.Ticks = Ticks;
  end
  
  subplot(212)
  imagesc(wrapToPi(expected_DCM_phase*pi/180)*180/pi)
  xlabel('Sensor index')
  ylabel('Sensor index')
  h = colorbar;
  ylabel(h,'\angle R (\circ) - After unwrapping')
  colormap parula
  h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
  if ~isempty(Ticks)
    h.Ticks = Ticks;
  end
end

% Plot the eigenvalues
% --------------------
figure(3002);clf
stem(10*log10(mean_eigval),'filled','LineWidth',2)
xlabel('Eigenvalue index')
ylabel('Normalized eigenvalue (dB)')
title(sprintf('2D sim: DOA = %2.2f deg.',DOA(1)*180/pi))
%   title('Mean eigenvalues of R')
grid on
grid minor

% Plot eigenvaectors spectrum (or MUSIC pseudo-spectrum)
% ----------------------------------------------------
theta_sv = linspace(-90,90,2048)*pi/180;

% Create steering vectors
sv_params = [];
%   phase_center = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
sv_params.src.y_pc = -phase_center(2,:).';
sv_params.src.z_pc = -phase_center(3,:).';
sv_params.src.fc   = fc;

A = @(theta,param)steering_mtx(theta,sv_params);
SVs = A(theta_sv,sv_params)./sqrt(Nc);

% MUSIC pesudo-spectrum
[V, D] = eig(DCM);
[D, D_idxs] = sort(real(diag(D)),'descend');
V = V(:,D_idxs);

N = 1;
f_music = 1./sum(abs(V(:,N+1:end)' * SVs).^2,1);
f_music = 10*log10(f_music./max(abs(f_music)));

figure(3003);clf
plot(theta_sv*180/pi,f_music,'b')
xlim(DOA(1)*180/pi+[-60  60])
xlabel('\theta ^\circ')
ylabel('Pseudo-power (dB)')
title(sprintf('2D sim: DOA = %2.2f deg.',DOA(1)*180/pi))
grid on
grid minor

% Plot array geometry
% -------------------
figure(3004);clf
%   ypc = -phase_center(2,:).';
%   zpc = -phase_center(3,:).';
plot(ypc,zpc ,'b*','LineWidth',2)
if max(ypc) > min(ypc)
  xlim([min(ypc)  max(ypc)])
end
if max(zpc) > min(zpc)
  ylim([min(zpc)  max(zpc)])
end

grid on
xlabel('y-axis')
ylabel('z-axis')
title(sprintf('Phase centers of the sensors for %2.2f deg. roll angle',theta_roll*180/pi))

if 0
  % Plot the eigenvalues similar to 1D simulator (i.e. different
  % range-bins are the snapshots, even though there is correlation
  % between range snapshots here, which is ignored in the 1D case)
  sample_data = squeeze(sim_data(:,1,:)).';
  DCM = 1/size(sample_data,2) * sample_data*sample_data';
  eigenvalue = real(sort(eig(DCM),'descend'));
  
  figure(3);clf
  stem(10*log10(eigenvalue),'filled','LineWidth',2)
  xlabel('Eigenvalue index')
  ylabel('Normalized eigenvalue (dB)')
  title('1D-equivalent eigenvalues of R')
  grid on
  grid minor
end

if param.debug_level >= 3
  return
end
