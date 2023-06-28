function A = steering_mtx(theta,sv_params)

physical_constants;

% -------------------------------------------------------------------------
% Uncalibrated array parameters
% -------------------------------------------------------------------------
doa  = theta;%sv_params.src.actual_doa;
y_pc = sv_params.src.y_pc;
z_pc = sv_params.src.z_pc;

Nc = length(y_pc);
Q  = length(doa);

lambda = c/sv_params.src.fc;
k      = 2*pi/(lambda/2);

% if isfield(sv_params,'array_gain_model_fh') && ~isempty(sv_params.array_gain_model_fh)
%  array_gain_model_fh =  sv_params.array_gain_model_fh;
% else
%   array_gain_model_fh = @(x) (exp(-x/20));
% end

% phase_equaliz = [-72.8 157.9 98.2 -33.3 -107.1 -5.9 97.1].'*pi/180;
% gain_equaliz = [4 2.7 0.1 2.4 2.4 3 -1.1].';
% gain_equaliz = 10.^(gain_equaliz./10);
% -------------------------------------------------------------------------
% Calibration error parameters
% -------------------------------------------------------------------------
% Array errors
if isfield(sv_params,'extra_error_params') && ~isempty(sv_params.extra_error_params)
  extra_error_params = sv_params.extra_error_params;
  
  error_ypc      = extra_error_params.error_ypc; % in meters (i.e. already multiplied by lambda)
  error_zpc      = extra_error_params.error_zpc; % in meters (i.e. already multiplied by lambda)
  error_phase    = extra_error_params.error_phase ;
  error_g_s      = extra_error_params.error_g_s;
  error_g_p      = extra_error_params.error_g_p;
  error_g_offset = extra_error_params.error_g_offset;
else
  error_ypc      = zeros(Nc,1);
  error_zpc      = zeros(Nc,1);
  error_phase    = zeros(Nc,1);
  error_g_s      = zeros(Nc,1);
  error_g_p      = zeros(Nc,1);
  error_g_offset = zeros(Nc,1);
end

% Calibration parameters (for simulation only): subtract the actual array  
% errors from these errors. Ideally, you should get all zeros.
% In real data case, we don't know the actual array errors, but we know the
% estimate of them (after calibration).
if isfield(sv_params,'calib_params') && ~isempty(sv_params.calib_params)
  calib_params = sv_params.calib_params;
  
  calib_ypc      = calib_params.calib_ypc;
  calib_zpc      = calib_params.calib_zpc;
  calib_phase    = calib_params.calib_phase ;
  calib_g_s      = calib_params.calib_g_s;
  calib_g_p      = calib_params.calib_g_p;
  calib_g_offset = calib_params.calib_g_offset;
else
  calib_ypc      = zeros(Nc,1);
  calib_zpc      = zeros(Nc,1);
  calib_phase    = zeros(Nc,1);
  calib_g_s      = zeros(Nc,1);
  calib_g_p      = zeros(Nc,1);
  calib_g_offset = zeros(Nc,1);
end

if 1
  % Relative phase centers is what matters
  y_pc = y_pc - repmat(y_pc(1,:),[Nc  1]);
  z_pc = z_pc - repmat(z_pc(1,:),[Nc  1]);
%   y_pc = cumsum([0;diff(y_pc,[],1)],1,'forward');
%   z_pc = cumsum([0;diff(z_pc,[],1)],1,'forward');
  
  error_ypc = error_ypc - repmat(error_ypc(1,:),[Nc  1]);
  error_zpc = error_zpc - repmat(error_zpc(1,:),[Nc  1]);
%   error_ypc = cumsum([0;diff(error_ypc,[],1)],1,'forward');
%   error_zpc = cumsum([0;diff(error_zpc,[],1)],1,'forward');
  
  calib_ypc = calib_ypc - repmat(calib_ypc(1,:),[Nc  1]);
  calib_zpc = calib_zpc - repmat(calib_zpc(1,:),[Nc  1]);
%   calib_ypc = cumsum([0;diff(calib_ypc,[],1)],1,'forward');
%   calib_zpc = cumsum([0;diff(calib_zpc,[],1)],1,'forward');
end
% Determine the array steering matrix
% -------------------------------------------------------------------------
A = zeros(Nc,Q);
for doa_idx = 1:Q
  tmp_doa = doa(doa_idx);
  ky = k*sin(tmp_doa);
  kz = k*cos(tmp_doa);
  
  gain_error_offset = error_g_offset; % In dB
  gain_calib_offset = calib_g_offset; % In dB
  
%   gain_error_other = array_gain_model_fh(gain_error_exp); 
%   gain_calib_other = array_gain_model_fh(gain_calib_exp); 

  gain_error_other = error_g_s.*(sin(tmp_doa)-sin(error_g_p)).^2; % In dB
  gain_calib_other = calib_g_s.*(sin(tmp_doa)-sin(calib_g_p)).^2; % In dB

%   gain_error_other = (sin(error_g_s * tmp_doa)-sin(error_g_s .* error_g_p)).^2; % In dB
%   gain_calib_other = (sin(calib_g_s * tmp_doa)-sin(calib_g_s .* calib_g_p)).^2; % In dB
  
  gain_error_dB = gain_error_other + gain_error_offset;  
  gain_calib_dB = gain_calib_other + gain_calib_offset; 

  gain_error = 10.^(gain_error_dB./20);
  gain_calib = 10.^(gain_calib_dB./20);
  gain = gain_error./gain_calib;
  
  if 0
    % If you would like to make the center sensor the referenc sensor.
    error_ypc = error_ypc - repmat(error_ypc(ceil(Nc/2),:),[Nc 1]);
    error_zpc = error_zpc - repmat(error_zpc(ceil(Nc/2),:),[Nc 1]);
  end  
  
  % Adding the error_zpc*kz will not affect the final results, but we keep
  % it here because mathematically it should be added.
  phase_error = exp(1i*(error_ypc*ky - error_zpc*kz + error_phase));% + error_zpc*kz));
  phase_calib = exp(1i*(calib_ypc*ky - calib_zpc*kz + calib_phase));% + error_zpc*kz));
  phase = phase_error .* conj(phase_calib);
  
  A(:,doa_idx) = gain .* phase .* exp(1i*(y_pc*ky - z_pc*kz));
%   A(:,doa_idx) = diag(gain .* phase) * exp(1i*(y_pc*ky - z_pc*kz));
end
% A = (1/sqrt(Nc))*A; % Don't normalize when calibrating the array.

return