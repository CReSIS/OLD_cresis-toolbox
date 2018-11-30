function fk_data = fk_migration(data,time,freq,kx,eps_r,param)
% fk_data = fk_migration(data,time,freq,kx,eps_r,param)
%
% Author: William Blake, John Paden

physical_constants;

% Calculate variables for later use
dt = time(2)-time(1);
w = 2*pi*freq;

%% Shift to time zero
data = data.*repmat(exp(-1i*w*time(1)),[1 size(data,2)]);

%% Compensate for nonzero time start
v_p     = c/sqrt(eps_r(1));
dz      = v_p*time(1)/2;
k       = 2*w/v_p; % rad/m
z_shift = exp(1i*dz*sqrt((repmat(k,1,size(data,2))).^2-(repmat(kx,size(data,1),1)).^2));

data = data.*z_shift;

%% Move zero slow-time frequency to center
data = fftshift(data,2);
kx = fftshift(kx);

%% Create slow time window
num_subapertures = length(param.sar.sub_aperture_steering);
proc_oversample = (1+num_subapertures)/2;
Nx_out = size(data,2)/proc_oversample;
Hwindow = param.sar.st_wind(Nx_out).';

%% Preallocate output matrix
fk_data = zeros(size(data,1),size(data,2)/proc_oversample,num_subapertures,class(data));

%% Perform FK Migration
for idx = 1:length(eps_r)
  % Calculate kz
  if (idx == 1) || (eps_r(idx) ~= eps_r(idx-1))
    v_p     = c/sqrt(eps_r(idx)); % m/s - Propogation Velocity
    dz      = v_p*dt/2;
    k       = 2*w/v_p; % rad/m
    z_shift = exp(1i*dz*sqrt((repmat(k,1,size(data,2))).^2-(repmat(kx,size(data,1),1)).^2));
  end
  
  for subap = 1:num_subapertures
    good_idxs = 1 + round((subap-1)/2*Nx_out) + (0:Nx_out-1);
    if 0
      % Apply frequency dependent slow-time window
      error('Not supported');
      
      % Calculate fields at time zero
      data_time_zero = mean(data);
      
      % Apply subaperture slow-time window and take slow-time IFFT
      fk_data(idx,:,subap) = ifft(ifftshift(data_time_zero(good_idxs)));
    else
      % Calculate fields at time zero
      data_time_zero = mean(data);
      
      % Apply subaperture slow-time window and take slow-time IFFT
      fk_data(idx,:,subap) = ifft(ifftshift(data_time_zero(good_idxs) .* Hwindow));
    end
  end
  
  % Apply phase correction so that SAR processing only corrects phase to
  % closest approach rather than all the way back to the sensor.
  fk_data(idx,:,:) = fk_data(idx,:,:) * exp(-1i*2*pi*freq(1)*time(idx));
  
  % Spatially shift plane waves
  data    = data.*z_shift;
end

% figure; imagesc(lp(ifft2(g_data)))
% figure; imagesc(lp(ifft2(tmp_fk_data)))
% figure; imagesc(lp(fk_data))

return