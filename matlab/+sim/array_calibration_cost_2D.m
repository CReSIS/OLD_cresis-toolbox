function total_ac_cost = array_calibration_cost_2D(est_errors,ac_cost_params)
% This function calculates the cost associated with a combination of array
% errors. The cost is then used to calibrate the array for these errors.
% This function is called from .crosstrack_example.m.
%
% Note that in 2D simulator, we would like to minimize the cost over all
% all range-lines (i.e. over all available data set).

% Author: Mohanad Al-Ibadi and Sravya Athinarapu

physical_constants;

% -----------------------------------------------------------------------
%               Call steering matrix function
% -----------------------------------------------------------------------

% Define the relevant parameters
roll_angles         = ac_cost_params.roll_angles;  
sv_params.src.y_pc = ac_cost_params.y_pc;
sv_params.src.z_pc = ac_cost_params.z_pc;
% sv_params.src.fc   = ac_cost_params.fc;
if ~isfield(ac_cost_params,'ref_chan') || isempty(ac_cost_params.ref_chan)
  ref_chan = 1;
else
  ref_chan = ac_cost_params.ref_chan;
end

if ~isfield(ac_cost_params,'Nb') || isempty(ac_cost_params.Nb)
  Nb = 1;
else
  Nb = ac_cost_params.Nb;
end
BW = ac_cost_params.BW;
fc = ac_cost_params.fc;
fc = fc + BW*[0:floor((Nb-1)/2), -floor(Nb/2):-1]/Nb;

% array_param        = ac_cost_params.array_param;
if isfield(ac_cost_params,'sim_data') && ~isempty(ac_cost_params.sim_data)
  sim_data        = ac_cost_params.sim_data; % Nt*Nx*Nc matrix
  Nt = size(sim_data,1);
  Nx = size(sim_data,2);
elseif isfield(ac_cost_params,'Rxx_all') && ~isempty(ac_cost_params.Rxx_all)
  Rxx_all = ac_cost_params.Rxx_all; % Nt*Nc*Nc*Nx matrix
%   Nt = size(Rxx_all,1);
%   Nx = size(Rxx_all,4);
else
  sprintf('No data or DCM is provided')
end
Nc     = length(sv_params.src.y_pc(:,1));

% est_errors is an Nc*6 matrix, where each row represents the errors
% associated with one sensor AND each column represents one error type for
% all the sensors.
est_errors = reshape(est_errors,[Nc, length(est_errors)/Nc]);
% est_errors_for_display = est_errors;
% est_errors_for_display(:,[3 5]) = est_errors_for_display(:,[3 5])*180/pi

% extra_error_params.error_ypc      = est_errors(:,1);% *lambda;
% extra_error_params.error_zpc      = est_errors(:,2);% *lambda;
extra_error_params.error_phase    = est_errors(:,3);
extra_error_params.error_g_s      = est_errors(:,4);
extra_error_params.error_g_p      = est_errors(:,5);
extra_error_params.error_g_offset = est_errors(:,6);

% sv_params.extra_error_params      = extra_error_params;

% if isfield(ac_cost_params,'array_gain_model_fh') && ~isempty(ac_cost_params.array_gain_model_fh)
%  sv_params.array_gain_model_fh =  ac_cost_params.array_gain_model_fh;
% else
%   sv_params.array_gain_model_fh = @(x) (exp(-x/20));
% end

sv = @(theta) steering_mtx(theta,sv_params);

% -----------------------------------------------------------------------
%                    Calculate the cost
% -----------------------------------------------------------------------
if 1
  % Use when calibrating using real data.
  % For a given DOA, find all the rnage-lines that have range-bins containing   
  % this DOA. Then use these range-lines to form a DCM for that DOA. 
  m = length(roll_angles);
  ac_cost = 0;
  for doa_idx = 1:m    
      theta_roll = roll_angles(doa_idx);
      theta_roll = theta_roll(~isnan(theta_roll));
      if ~isempty(theta_roll)
        % Phase centers for this range-line. First, rotate the array to the
        % nominal orientation (0 degree roll). To rotate, you should first
        % make the center of the array as the point of rotation, and then
        % change back to your nominal reference sensor.
        rotation_mtx = [cos(-theta_roll), -sin(-theta_roll);sin(-theta_roll), cos(-theta_roll)];
        phase_center(1,:) = ac_cost_params.y_pc(:,doa_idx);
        phase_center(2,:) = ac_cost_params.z_pc(:,doa_idx);
%         phase_center(1,:) = ac_cost_params.y_pc(:,doa_idx) - ac_cost_params.y_pc(ceil(Nc/2),doa_idx);
%         phase_center(2,:) = ac_cost_params.z_pc(:,doa_idx) - ac_cost_params.z_pc(ceil(Nc/2),doa_idx);
       
        rotated_phase_center = rotation_mtx*phase_center;
%         rotated_phase_center(1,:) = rotated_phase_center(1,:) - rotated_phase_center(1,ref_chan);
%         rotated_phase_center(2,:) = rotated_phase_center(2,:) - rotated_phase_center(2,ref_chan);
        
        sv_params.src.y_pc = rotated_phase_center(1,:).'; 
        sv_params.src.z_pc = rotated_phase_center(2,:).';
        
%         sv_params.src.y_pc = ac_cost_params.y_pc(:,doa_idx); 
%         sv_params.src.z_pc = ac_cost_params.z_pc(:,doa_idx); 
               
        extra_error_params.error_ypc = est_errors(:,1);% * cos(theta_roll); 
        extra_error_params.error_zpc = est_errors(:,2);% * sin(theta_roll); 
        
        ac_cost_tmp = 0;
        for nb_idx = 1:Nb
          sv_params.src.fc   = fc(nb_idx);
          sv_params.extra_error_params = extra_error_params;
          sv = @(theta) steering_mtx(theta,sv_params);
%           theta_roll = unique(theta_roll);
%           A = sv(0);  
          A = sv(theta_roll); 
          
          R = squeeze(Rxx_all{doa_idx}((nb_idx-1)*Nc+(1:Nc),:));
          if 1
            % MLE cost function: Gives the best results
            proj_mtx = A*((A'*A)\A');
            I = eye(size(proj_mtx));
%           ac_cost_tmp = ac_cost_tmp + abs((trace((I-proj_mtx)*R))); % Mohanad's ==> Remove the - from total_ac_cost if you use it
%             ac_cost_tmp = ac_cost_tmp + abs((trace(proj_mtx*R)));% Wax's
            ac_cost_tmp = ac_cost_tmp + abs(sum(sum(proj_mtx .* R.'))); % John's (another form of Wax's)
          elseif 0
            % MAP cost function (Pr(x,error)=Pr(x|error) * Pr(error))
            % STILL NOT COMPLETE
            proj_mtx = A*((A'*A)\A');
            I = eye(size(proj_mtx));
%             L_mle = trace(proj_mtx*R);
            L_mle = abs(sum(sum(proj_mtx*R.')));
            
            % Aprior density of y and z errors
            dy = extra_error_params.error_ypc;
            dz = extra_error_params.error_zpc;
            sigma_y = abs(10/100)^2; % in meters
            sigma_z = abs(20/100)^2; % in meters
            L_y = log(1/(sqrt(sigma_y)) * exp(-(norm(dy,2)^2)./(2*sigma_y)));
            L_z = log(1/(sqrt(sigma_z)) * exp(-(norm(dz,2)^2)./(2*sigma_z)));
            L_aprior = L_y + L_z;
            
            % Posterior density 
            L_posterior = L_mle + L_aprior;
%             L_posterior = 10*log10(abs(L_mle)) + 10*log10(abs(L_aprior));
            ac_cost_tmp = ac_cost_tmp - 10*log10(abs(L_posterior));
          elseif 0
            % MUSIC-like cost function
            [V,D] = eig(R);
            [D D_idx] = sort(real(diag(D)),'descend');
            Ds = D(1);
            Dn = D(2:end);
            Vs = V(:,D_idx(1));
            Vn = V(:,D_idx(2:end));
            P_v = Vs*((Vs'*Vs)\Vs');
            P_n = Vn*((Vn'*Vn)\Vn');
            ac_cost_tmp = ac_cost_tmp + norm(P_n*A)^2;
            %           ac_cost_tmp = -norm(P_v*A)^2;
          end
        end
        ac_cost = ac_cost + ac_cost_tmp; % in dB  
      end
  end
  total_ac_cost = -10*log10(ac_cost);
  
elseif 0
  % Use for 2D simulations
  % The total cost is the mean of the costs over range-lines (eventually in
  % dB) -- This method seems to give a better results
  ac_cost = zeros(Nx,1);
  for lineIdx = 1:Nx %length(array_param.lines)
    lineIdx_idx = lineIdx; % array_param.lines(lineIdx);   
    ac_cost_line = 0;
    for binIdx = 1:Nt %length(array_param.bins)
      binIdx_idx = binIdx; % array_param.bins(binIdx);
      % ADD A CONSTRAINT TO NOT PROCESS DOAS WITH 2 DEGREES FROM EACH OTHER
      doa = actual_doa(binIdx,:,lineIdx);
      doa = doa(~isnan(doa));
      if ~isempty(doa)
        % Estimate the DCM
        if exist('sim_data','var')
          dataSample  = sim_data(binIdx_idx+array_param.bin_rng,lineIdx_idx+array_param.rline_rng,:);
          dataSample  = reshape(dataSample,[length(array_param.bin_rng)*length(array_param.rline_rng) Nc]);
          array_data  = dataSample.';
          R = (1/size(array_data,2)) * (array_data * array_data');
        else
          R = squeeze(Rxx_all(binIdx,:,:,lineIdx));
        end
        
        doa = unique(doa);
        A = sv(doa);
        proj_mtx = A*((A'*A)\A');
        I = eye(size(proj_mtx));
        
        ac_cost_line = ac_cost_line + 10*log10(abs((trace((I-proj_mtx)*R)))); % in dB to avoide numerical issues
      end
    end
    ac_cost(lineIdx) = 10^(ac_cost_line/10); % In linear units
  end
  
  total_ac_cost = 10*log10(mean(ac_cost));
  
elseif 0
  % Use for 2D simulations
  % The total cost is the sum of the costs from each range-line (in dB)
  ac_cost = zeros(length(array_param.lines),1);
  for lineIdx = 1:length(array_param.lines)
    lineIdx_idx = array_param.lines(lineIdx);
    
    ac_cost_line = 0;
    for binIdx = 1:length(array_param.bins)
      binIdx_idx = array_param.bins(binIdx);
      
      doa = actual_doa(binIdx,:,lineIdx_idx);
      doa = doa(~isnan(doa));
      if ~isempty(doa)
        % Estimate the DCM
        dataSample  = sim_data(binIdx_idx+array_param.bin_rng,lineIdx_idx+array_param.rline_rng,:);
        dataSample  = reshape(dataSample,[length(array_param.bin_rng)*length(array_param.rline_rng) Nc]);
        array_data  = dataSample.';
        R = (1/size(array_data,2)) * (array_data * array_data');
        
        doa = unique(doa);
        A = sv(doa);
        proj_mtx = A*((A'*A)\A');
        I = eye(size(proj_mtx));
        
        ac_cost_line = ac_cost_line + log10(abs((trace((I-proj_mtx)*R))));
      end
    end
    ac_cost(lineIdx) = ac_cost_line;
  end
  
  total_ac_cost = sum(ac_cost);
end

return