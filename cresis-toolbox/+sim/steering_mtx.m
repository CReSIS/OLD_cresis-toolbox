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

% -------------------------------------------------------------------------
% Calibration error parameters
% -------------------------------------------------------------------------
if isfield(sv_params,'extra_error_params') && ~isempty(sv_params.extra_error_params)
  extra_error_params = sv_params.extra_error_params;
  
  error_ypc      = extra_error_params.error_ypc;
  error_zpc      = extra_error_params.error_zpc;
  error_phase    = extra_error_params.error_phase;
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

% Determine the array steering matrix
% -------------------------------------------------------------------------
A = zeros(Nc,Q);
for doa_idx = 1:Q
  tmp_doa = doa(doa_idx);
  gain_error_exp = error_g_s.*(sin(tmp_doa)-sin(error_g_p)).^2 + error_g_offset;
%   gain_error  = 10.^(gain_error_exp./20);
%   gain_error  = 1 + gain_error_exp;
  gain_error  = exp(-gain_error_exp./2);
%   gain_error  = 10.^(-gain_error_exp./20);
  ky = (y_pc+error_ypc)*k*sin(tmp_doa);
  kz = (z_pc+error_zpc)*k*cos(tmp_doa);
  phase_error = error_phase;
  
  A(:,doa_idx) = gain_error .* exp(1i*(ky-kz+phase_error));
end
A = (1/sqrt(Nc))*A;

return