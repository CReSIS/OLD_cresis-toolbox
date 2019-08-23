function ac_cost = array_calibration_cost_1D(est_errors,ac_cost_params)
% This function calculates the cost associated with a combination of array
% errors. The cost is then used to calibrate the array for these errors.
% This function is called from doa_example.m
%
% Author: Sravya
% Modified/updated by Mohanad Al-Ibadi
%%
physical_constants;
  % -----------------------------------------------------------------------
  % Call steering matrix function
  % -----------------------------------------------------------------------
if 1
  % Define the relevant parameters
  actual_doa         = ac_cost_params.actual_doa;
  sv_params.src.y_pc = ac_cost_params.y_pc;
  sv_params.src.z_pc = ac_cost_params.z_pc;
  sv_params.src.fc   = ac_cost_params.fc;
  DCM                = ac_cost_params.DCM;
  
  Nc     = length(sv_params.src.y_pc);
  lambda = c/sv_params.src.fc;
  
  % est_errors is an Nc*6 matrix, where each row represents the errors
  % associated with one sensor AND each column represents one error type for
  % all the sensors.
  est_errors = reshape(est_errors,[Nc, 6])
  
  extra_error_params.error_ypc      = est_errors(:,1)*lambda;
  extra_error_params.error_zpc      = est_errors(:,2)*lambda;
  extra_error_params.error_phase    = est_errors(:,3);
  extra_error_params.error_g_s      = est_errors(:,4);
  extra_error_params.error_g_p      = est_errors(:,5);
  extra_error_params.error_g_offset = est_errors(:,6);
  
  sv_params.extra_error_params      = extra_error_params;
  sv = @(theta) sim.steering_mtx(theta,sv_params);
  ac_cost = 0;
  for  t_idx = 0:ceil(length(actual_doa)/2)-1    
    doa = actual_doa([t_idx+1,end-t_idx]).';
    if doa(1) == doa(2)
      doa = doa(1);
    end
    
    R = DCM{t_idx+1};
    
    A = sv(doa);
    proj_mtx = A*((A'*A)\A');
    I = eye(size(proj_mtx));
    
    ac_cost = ac_cost + log10(abs((trace((I-proj_mtx)*R))));
  end
  
  % ac_cost = Nc*M*ac_cost;
  
  return
  
  % -----------------------------------------------------------------------
  % Generate steering matrix locally
  % -----------------------------------------------------------------------
elseif 0
  
  
  % Define the relevant parameters
  actual_doa = ac_cost_params.actual_doa;
  y_pc       = ac_cost_params.y_pc;
  z_pc       = ac_cost_params.z_pc;
  fc         = ac_cost_params.fc;
  DCM        = ac_cost_params.DCM;
  
  Nc     = length(y_pc);
  lambda = c/fc;
  k      = 2*pi/(lambda/2);
  
  % est_errors is an Nc*6 matrix, where each row represents the errors
  % associated with one sensor. OR, each column represents one error type for
  % all the sensors.
  est_errors = reshape(est_errors,[Nc, 6]);
  
  est_error_ypc      = est_errors(:,1)*lambda;
  est_error_zpc      = est_errors(:,2)*lambda;
  est_error_phase    = est_errors(:,3);
  est_error_g_s      = est_errors(:,4);
  est_error_g_p      = est_errors(:,5);
  error_g_offset     = est_errors(:,6);
  
  ac_cost = 0;
  for  t_idx = 0:ceil(length(actual_doa)/2)-1
    
    doa = actual_doa([t_idx+1,end-t_idx]).';
    if doa(1) == doa(2)
      doa = doa(1);
    end
    
    R = DCM{t_idx+1};
    
    for doa_idx = 1:length(doa)
      tmp_doa = doa(doa_idx);
      gain_error_exp = est_error_g_s.*(sin(tmp_doa)-sin(est_error_g_p)).^2 + error_g_offset;
      gain_error  = 10.^(gain_error_exp./20);
      ky = (y_pc+est_error_ypc)*k*sin(tmp_doa);
      kz = (z_pc+est_error_zpc)*k*cos(tmp_doa);
      phase_error = est_error_phase;
      
      A(:,doa_idx) = gain_error .* exp(1i*(ky-kz+phase_error));
    end
    A = (1/sqrt(Nc))*A;
    proj_mtx = A*((A'*A)\A');
    I = eye(size(proj_mtx));
    
    ac_cost = ac_cost + log10(abs((trace((I-proj_mtx)*R))));
  end
  
  % ac_cost = Nc*M*ac_cost;
  return
  
end


