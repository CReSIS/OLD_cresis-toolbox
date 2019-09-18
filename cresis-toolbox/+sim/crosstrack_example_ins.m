% script sim.crosstrack_example_ins
%
% Example inertial navigation system (INS) simulation; specifically roll.
% Calls crosstrack.m and crosstrack_rmse.m.
%
% Author: Mohanad Al-Ibadi, Sravya Athinarapu, Sean Holloway, John Paden,
% Theresa Stumpf
%
% See also: sim.crosstrack.m, sim.crosstrack_data.m,
% sim.crosstrack_example_*.m, sim.crosstrack_rmse.m

%% Studying the effect roll and the actual phase of the DCM: for array calibration prurpos
% ----------------------------------------------------------------------------------------

%% Setup simulation parameters
physical_constants;
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

lambda = c/fc;

DOA = [20]*pi/180;

% Main target location
z_pos = -1500;
Range = abs(z_pos)/cos(DOA); % Range to the main target
y_pos = Range*sin(DOA);      % y-axis of the main target

param.src.f0                      = fc-BW/2;
param.src.f1                      = fc+BW/2;
param.src.t0                      = 2*(abs(z_pos)-500)/c;
param.src.t1                      = 2*(abs(z_pos)+1000)/c;
if param.src.t0 < 0
  param.src.t0 = 0;
end
param.src.ft_func                 = @(t) tukeywin_cont(t * BW);
param.src.ft_wind                 = @(N) hanning(N);
param.src.lever_arm.fh            = @sim.lever_arm_example;

% Nsep: number of lambda/4 steps between sensors
Nsep = 1;
param.src.lever_arm.fh_args       = {[], 1, 1:Nc, [0; Nsep*c/fc/2; 0]};
phase_center = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});

% Change reference sensor if you want to
ref_chan = 1;
phase_center(2,:) = phase_center(2,:) - phase_center(2,ref_chan);
phase_center(3,:) = phase_center(3,:) - phase_center(3,ref_chan);

% Account for roll angle. Default is 0
%   theta_roll = DOA;
theta_roll = 0 * pi/180;
for chan_idx = 1:Nc
  d_chan_ref = sqrt(abs(diff(phase_center(2,[ref_chan  chan_idx])))^2 + abs(diff(phase_center(3,[ref_chan  chan_idx])))^2);
  phase_center(2,chan_idx) = d_chan_ref .* cos(theta_roll);
  phase_center(3,chan_idx) = d_chan_ref .* sin(theta_roll);
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
surf_param.y = y_pos;

% y_range should >= maximum y (i.e. targets shoulb be inside the imaged
% swath)
y_range = 2*Range;
if y_range == 0
  surf_param.y_range = [-1000 1000];
else
  surf_param.y_range = [-1.5*y_range 1.5*y_range];
end
param.monte.target_param{1} = surf_param;

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
  % Phase difference between sensors that we expect to see in the phase
  % of the DCM. It is calculated from geometry.
  %     ypc = -phase_center(2,:).';
  %     zpc = -phase_center(3,:).';
  
  % DOA_roll is the DOA wrt the range vector connecting the array and target
  % (i.e. the new nadir).
  if theta_roll >= DOA
    DOA_roll = theta_roll - DOA;
  elseif theta_roll < DOA
    DOA_roll = DOA - theta_roll;
  end
  
  clear max_expected_phase
  for chan_idx = 1:Nc
    % Distance from the a given sensor to the reference sensor
    L_chan      = sqrt(abs(diff(ypc([ref_chan  chan_idx])))^2 + abs(diff(zpc([ref_chan  chan_idx])))^2);
    % Array length in the direction of DOA
    L_doa   = L_chan*cos(DOA_roll);
    % Extra delay with respect to the reference sensor
    L_extra = L_chan*sin(DOA_roll);
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

% Calculate TBP across the array in the direction of DOA
TBP = 2*L_doa/c * BW;

sprintf('TBP across the array in the direction of %2.1f deg. is %2.2f',DOA*180/pi,TBP)

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
% Unwrap along both dimsnsions
DCM_phase = unwrap(DCM_phase,[],2);
DCM_phase = DCM_phase - diag(diag(DCM_phase));

DCM_phase = unwrap(DCM_phase,[],1);
DCM_phase = DCM_phase - diag(diag(DCM_phase));

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
sv_params.src.y_pc = ypc;
sv_params.src.z_pc = zpc;
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
