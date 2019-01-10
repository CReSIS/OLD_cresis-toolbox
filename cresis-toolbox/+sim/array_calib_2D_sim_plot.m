%     actual_doa = 0;    
    roll_angle = 10* pi/180;
    actual_doa = roll_angle;
    N = 1;
    [~,roll_angle_idx] = min(abs(roll_angle-roll_angles));
    
    phase_unwrapping = 1;
    k = 4*pi/(c/fc);
    % Rotate phase centers and errors back to nominal orientation
    rotation_mtx = [cos(-roll_angle), -sin(-roll_angle);sin(-roll_angle), cos(-roll_angle)];
    rotated_est_error = rotation_mtx*[ac_est_error(:,1) ac_est_error(:,2)].';
    
    calib_params = [];
    calib_params.calib_ypc      = rotated_est_error(1,:).';
    calib_params.calib_zpc      = rotated_est_error(2,:).';
    calib_params.calib_phase    = ac_est_error(:,3);% + k*ac_est_error(:,2);
    calib_params.calib_g_s      = ac_est_error(:,4);
    calib_params.calib_g_p      = ac_est_error(:,5);
    calib_params.calib_g_offset = ac_est_error(:,6);
    
    rotated_error = rotation_mtx*[error_ypc error_zpc].';
    ypc_error = rotated_error(1,:).';
    zpc_error = rotated_error(2,:).';
    
    extra_error_params.error_ypc = ypc_error;
    extra_error_params.error_zpc = zpc_error;
    
%     ideal_ypc = -ac_cost_params.y_pc(:,roll_angle_idx);
%     ideal_zpc = -ac_cost_params.z_pc(:,roll_angle_idx);
%     ideal_ypc = ac_cost_params.y_pc(:,roll_angle_idx);
%     ideal_zpc = ac_cost_params.z_pc(:,roll_angle_idx);
        
    rotated_pc = rotation_mtx*[ ac_cost_params.y_pc(:,roll_angle_idx) ac_cost_params.z_pc(:,roll_angle_idx)].';
    ideal_ypc = rotated_pc(1,:).';
    ideal_zpc = rotated_pc(2,:).';
    
    uncalib_ypc = ideal_ypc + Err(:,1);
    uncalib_zpc = ideal_zpc + Err(:,2);

    calib_ypc = uncalib_ypc - calib_params.calib_ypc ;
    calib_zpc = uncalib_zpc - calib_params.calib_zpc;% - ac_est_error(:,2);
    
    if 0
      % Debug: calculate the angle between y and z coordinates from geometry.
      % This angle should match the roll angle (or DOA defined above).
      DOA_hat_uncalib = nanmean(atan(uncalib_zpc./uncalib_ypc));
      DOA_hat_calib = nanmean(atan(calib_zpc./calib_ypc));
      
      h1 = sprintf('\nActual roll angle is %2.2f\n',roll_angles(doa_idx)*180/pi);
      h2 = sprintf('Uncalib. array (average) roll angle from geometry is %2.2f deg. \n',DOA_hat_uncalib*180/pi);
      h3 = sprintf('Calib. array (average) roll angle from geometry is %2.2f deg. \n',DOA_hat_calib*180/pi);
      disp(h1)
      disp(h2)
      disp(h3)
    end
    
    %% Plot phase centers of the sensors
    if 1
        ideal_ypc_plot = ac_cost_params.y_pc(:,roll_angle_idx);
        deal_zpc_plot = ac_cost_params.z_pc(:,roll_angle_idx);
        
        rotated_Err = [cos(roll_angle), -sin(roll_angle);sin(roll_angle), cos(roll_angle)] * [Err(:,1) Err(:,2)].';
        uncalib_ypc_plot = ideal_ypc_plot + rotated_Err(1,:).';
        uncalib_zpc_plot = ideal_zpc_plot + rotated_Err(2,:).';
        
        calib_ypc_plot = uncalib_ypc_plot - 
        
      figure(9998);clf
      hold on
      % Ideal phace centers
      plot(ideal_ypc_plot,deal_zpc_plot,'k*','LineWidth',1.5)
      % Actual phase centers
      plot(uncalib_ypc_plot,uncalib_zpc_plot,'b+','LineWidth',3)
      % Calibrated phase centers
      plot(calib_ypc,calib_zpc ,'rx','LineWidth',1.5)
%       xlim([min(min(uncalib_ypc),min(calib_ypc))  max(max(uncalib_ypc),max(calib_ypc))])
%       ylim([min(min(uncalib_zpc),min(calib_zpc))  max(max(uncalib_zpc),max(calib_zpc))])
      
      grid on
      xlabel('y-axis')
      ylabel('z-axis')
      title('Phase centers, PCs, of the sensors')
      legend('Ideal PCs','Uncalib. PCs','Calib. PCs','Location','southeast')
      %   legend('Actual PCs','Calib. PCs','calb. PCs projected onto y-axis','Location','best')
    end
    %% Plot DCM magnitude and angle
    if 1
      % 1- Before calibration
      % --------------------------
      if 1
        % L0:array length, L_doa:array length projected in the direction of
        % DOA, L_extra:extra distance the signal travels from sensor 1 to
        % sensor Nc (from which we determine the maximum phase difference
        % accross the array for a target at angle DOA). The DCM should show the
        % same phase between the first and last sensors (Sanity check).
        % Remember that the phase centers for a roll=DOA are projected onto the
        % y-axis (the nominal array axis) and now the array is no longer orthogomal
        % to the range vector in the direction of DOA. So, the new array length
        % is shorter and the maximum phase accross the array from a target at DOA
        % is the two way phase (i.e. 4*pi/lambda * L_extra).
        
        % 1) Maximum phase of the array before calibration
        %   L0 = lambda/4*(Nc-1);
        L0 = sqrt(abs(diff(ideal_ypc([1 end])))^2 + abs(diff(ideal_zpc([1 end])))^2);
%         L0 = sqrt(abs(diff(phase_center(2,[1 end])))^2 + abs(diff(phase_center(3,[1 end])))^2);
       
        % The DOA here is measured wrt vector orthogonal to the array axis.
%         doa_roll = actual_doa-roll_angle;
        doa_roll = roll_angle;
        L_doa  = L0*cos(doa_roll);
        L_extra = L0*sin(doa_roll);
        max_actual_expected_phase = 4*pi/lambda * L_extra * 180/pi;
        h1 = sprintf('\nMaximum phase (ideal) is %4.2f deg. (%4.2f after wraping to pi) \n',max_actual_expected_phase,wrapToPi(max_actual_expected_phase*pi/180)*180/pi);
        
        % 2) Maximum phase of the array after calibration
        L0 = sqrt(abs(diff(calib_ypc([1 end])))^2 + abs(diff(calib_zpc([1 end])))^2);
        L_doa  = L0*cos(doa_roll);
        L_extra = L0*sin(doa_roll);
        max_new_expected_phase = 4*pi/lambda * L_extra * 180/pi;
        h2 = sprintf('\nMaximum phase after calib. is %4.2f deg. (%4.2f after wraping to pi) \n',max_new_expected_phase,wrapToPi(max_new_expected_phase*pi/180)*180/pi);
        
        disp(h1)
        disp(h2)
      end
      
      DataSample = squeeze(Data{roll_angle_idx});      
      Rxx_uncalib = Rxx_all{roll_angle_idx};
      
      figure(9990);clf
      Title = sprintf('Roll = %2.2f deg.',roll_angles(roll_angle_idx)*180/pi);      
      t = suptitle(strcat('Before calib.: ',Title));
      
      subplot(211)
      abs_Rxx_uncalib = 10*log10(abs(Rxx_uncalib)./max(abs(Rxx_uncalib(:))));
      imagesc(abs_Rxx_uncalib)
      h = colorbar;
      colormap parula
      ylabel(h,'Normalized |R| (dB)')
      %   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      
      subplot(212)
      phase_Rxx_uncalib = angle(Rxx_uncalib);
      if phase_unwrapping
        phase_Rxx_uncalib = unwrap(phase_Rxx_uncalib,[],1);
        phase_Rxx_uncalib = phase_Rxx_uncalib - diag(diag(phase_Rxx_uncalib));
        
        phase_Rxx_uncalib = unwrap(phase_Rxx_uncalib,[],2);
        phase_Rxx_uncalib = phase_Rxx_uncalib - diag(diag(phase_Rxx_uncalib));
      end
      imagesc(phase_Rxx_uncalib*180/pi)
      h = colorbar;
       colormap parula
      ylabel(h,'\angle{R} (\circ)')
%       h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
%       Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
%       if ~isempty(Ticks)
%         h.Ticks = Ticks;
%       end
      
      % 2- After calibration
      % -------------------------
      gain_dB = calib_params.calib_g_s.*(sin(roll_angles(roll_angle_idx)) - sin(calib_params.calib_g_p)).^2 + calib_params.calib_g_offset;
      phase_deg = (calib_params.calib_ypc*k*sin(roll_angles(roll_angle_idx)) - calib_params.calib_zpc*k*cos(roll_angles(roll_angle_idx)) + ...
        calib_params.calib_phase)*180/pi;
      
      % Correct the DCM using the estimated calibration parameters.
      % The following two options are mathematically identical
      if 1
        calib_vec = 10.^(-gain_dB./20) .*exp(-1i*phase_deg*pi/180);
        DataSample_calib = repmat(calib_vec,[1 size(DataSample,2)]) .* DataSample;
        Rxx_calib =  1/size(DataSample_calib,2) * DataSample_calib*DataSample_calib';
      else
        calib_vec = 10.^(gain_dB./20) .*exp(1i*phase_deg*pi/180);
        H = diag(calib_vec);
        DataSample_calib = inv(H) * DataSample; % flipud(DataSample);
        Rxx_calib =  1/size(DataSample_calib,2) * DataSample_calib*DataSample_calib';
        % OR, do this
        %     Rxx_calib = inv(H) * Rxx_uncalib *inv(H');
      end
      
      figure(9991);clf
      suptitle(strcat('After calib.: ',Title))
      
      subplot(211)
      abs_Rxx_calib = 10*log10(abs(Rxx_calib)./max(abs(Rxx_calib(:))));
      imagesc(abs_Rxx_calib)
      h = colorbar;
       colormap parula
      ylabel(h,'Normalized |R| (dB)')
      %   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      
      subplot(212)
      phase_Rxx_calib = angle(Rxx_calib);
      if phase_unwrapping
        phase_Rxx_calib = unwrap(phase_Rxx_calib,[],1);
        phase_Rxx_calib = phase_Rxx_calib - diag(diag(phase_Rxx_calib));
        
        phase_Rxx_calib = unwrap(phase_Rxx_calib,[],2);
        phase_Rxx_calib = phase_Rxx_calib - diag(diag(phase_Rxx_calib));
      end
      imagesc(phase_Rxx_calib*180/pi)
      h = colorbar;
       colormap parula
      ylabel(h,'\angle{R} (\circ)')
%       Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
%       if ~isempty(Ticks)
%         h.Ticks = Ticks;
%       end
      
      % 3) Before calibration with forward-backward averaging
      % ------------------------------------------------------
      if 1
        % FB averaging
        reflect_mtx = flipud(eye(size(Rxx_uncalib)));
        Rxx_uncalib_fb = (1/2)*(Rxx_uncalib + reflect_mtx * conj(Rxx_uncalib) * reflect_mtx);
      elseif 0
        % FB averaging AND spatial smoothing
        Rxx1 = 1/(size(DataSample(1:end-1,:),2))*DataSample(1:end-1,:) * DataSample(1:end-1,:)';
        Rxx2 = 1/(size(DataSample(2:end,:),2))*DataSample(2:end,:) * DataSample(2:end,:)';
        
        % Apply FB averaging for each subarray
        reflect_mtx = flipud(eye(size(Rxx1)));
        Rxx1_fb = (1/2)*(Rxx1 + reflect_mtx * conj(Rxx1) * reflect_mtx);
        Rxx2_fb = (1/2)*(Rxx2 + reflect_mtx * conj(Rxx2) * reflect_mtx);
        
        % Average the two DCMs
        Rxx_ss = 1/2*(Rxx1_fb+Rxx2_fb); % (Nc-1)-by-(Nc-1)
        
        % Handle the lost sensor (or DOF) such that the final DCM is Nc-by-Nc
        Rxx_tmp = zeros(Nc);
        Rxx_tmp(1:end-1,1:end-1) = Rxx_ss;
        Rxx_tmp(end,:) = Rxx_uncalib(end,:);
        Rxx_tmp(:,end) = Rxx_uncalib(:,end);
        
        Rxx_uncalib_fb = Rxx_tmp; % Nc-by-Nc matrix
      end
      
      figure(9992);clf
      suptitle(strcat('Before calib. and FB avging:',',', Title))
      %      suptitle(strcat('After array calib. and FB avging: ',frame_name(1:8),'\_',frame_name(10:11),'\_',frame_name(13:end),',',Title))
      
      subplot(211)
      imagesc(10*log10(abs(Rxx_uncalib_fb)./max(abs(Rxx_uncalib_fb(:)))))
      h = colorbar;
       colormap parula
      ylabel(h,'Normalized |R| (dB)')
      %   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      
      subplot(212)
      phase_Rxx_uncalib_fb = angle(Rxx_uncalib_fb);
      if phase_unwrapping
        phase_Rxx_uncalib_fb = unwrap(phase_Rxx_uncalib_fb,[],1);
        phase_Rxx_uncalib_fb = phase_Rxx_uncalib_fb - diag(diag(phase_Rxx_uncalib_fb));
        
        phase_Rxx_uncalib_fb = unwrap(phase_Rxx_uncalib_fb,[],2);
        phase_Rxx_uncalib_fb = phase_Rxx_uncalib_fb - diag(diag(phase_Rxx_uncalib_fb));
      end
      imagesc(phase_Rxx_uncalib_fb*180/pi)
      h = colorbar;
       colormap parula
      ylabel(h,'\angle{R} (\circ)')
%       Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
%       if ~isempty(Ticks)
%         h.Ticks = Ticks;
%       end
      
      % 4) After calibration with forward-backward averaging
      % ----------------------------------------------------
      if 1
        % FB averaging
        reflect_mtx = flipud(eye(size(Rxx_calib)));
        Rxx_calib_fb = (1/2)*(Rxx_calib + reflect_mtx * conj(Rxx_calib) * reflect_mtx);
      elseif 0
        % FB averaging AND spatial smoothing
        Rxx1 = 1/(size(DataSample_calib(1:end-1,:),2))*DataSample_calib(1:end-1,:) * DataSample_calib(1:end-1,:)';
        Rxx2 = 1/(size(DataSample_calib(2:end,:),2))*DataSample_calib(2:end,:) * DataSample_calib(2:end,:)';
        
        % Apply FB averaging for each subarray
        reflect_mtx = flipud(eye(size(Rxx1)));
        Rxx1_fb = (1/2)*(Rxx1 + reflect_mtx * conj(Rxx1) * reflect_mtx);
        Rxx2_fb = (1/2)*(Rxx2 + reflect_mtx * conj(Rxx2) * reflect_mtx);
        
        % Average the two DCMs
        Rxx_ss = 1/2*(Rxx1_fb+Rxx2_fb); % (Nc-1)-by-(Nc-1)
        
        % Handle the lost sensor (or DOF) such that the final DCM is Nc-by-Nc
        Rxx_tmp = zeros(Nc);
        Rxx_tmp(1:end-1,1:end-1) = Rxx_ss;
        Rxx_tmp(end,:) = Rxx_calib(end,:);
        Rxx_tmp(:,end) = Rxx_calib(:,end);
        
        Rxx_calib_fb = Rxx_tmp; % Nc-by-Nc matrix
      end
      
      figure(9993);clf
      suptitle(strcat('After calib. and FB avging:',',', Title))
      %      suptitle(strcat('After array calib. and FB avging: ',frame_name(1:8),'\_',frame_name(10:11),'\_',frame_name(13:end),',',Title))
      
      subplot(211)
      imagesc(10*log10(abs(Rxx_calib_fb)./max(abs(Rxx_calib_fb(:)))))
      h = colorbar;
       colormap parula
      ylabel(h,'Normalized |R| (dB)')
      %   h.Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
      
      subplot(212)
      phase_Rxx_calib_fb = angle(Rxx_calib_fb);
      if phase_unwrapping
        phase_Rxx_calib_fb = unwrap(phase_Rxx_calib_fb,[],1);
        phase_Rxx_calib_fb = phase_Rxx_calib_fb - diag(diag(phase_Rxx_calib_fb));
        
        phase_Rxx_calib_fb = unwrap(phase_Rxx_calib_fb,[],2);
        phase_Rxx_calib_fb = phase_Rxx_calib_fb - diag(diag(phase_Rxx_calib_fb));
      end
      imagesc(phase_Rxx_calib_fb*180/pi)
      h = colorbar;
       colormap parula
      ylabel(h,'\angle{R} (\circ)')
%       Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
%       if ~isempty(Ticks)
%         h.Ticks = Ticks;
%       end      
      
      % 5) Ideal case
      % --------------
      figure(9999);clf
      suptitle(strcat('Ideal case (no errors):', Title))
      subplot(211)
      imagesc(10*log10(abs(Rxx_ideal_all{roll_angle_idx})./max(abs(Rxx_ideal_all{roll_angle_idx}(:)))))
      h = colorbar;
       colormap parula
      ylabel(h,'Normalized |R| (dB)')
      
      subplot(212)
      phase_Rxx_ideal = angle(Rxx_ideal_all{roll_angle_idx});
      if phase_unwrapping
        phase_Rxx_ideal = unwrap(phase_Rxx_ideal,[],1);
        phase_Rxx_ideal = phase_Rxx_ideal - diag(diag(phase_Rxx_ideal));
        
        phase_Rxx_ideal = unwrap(phase_Rxx_ideal,[],2);
        phase_Rxx_ideal = phase_Rxx_ideal - diag(diag(phase_Rxx_ideal));
      end
      imagesc(phase_Rxx_ideal*180/pi)
      h = colorbar;
       colormap parula
      ylabel(h,'\angle{R} (\circ)')
%       Ticks = [ceil(h.Limits(1)):(-ceil(h.Limits(1))+floor(h.Limits(2)))/4:floor(h.Limits(2))];
%       if ~isempty(Ticks)
%         h.Ticks = Ticks;
%       end            
    end
    
%     return    
    
    %% Eigenvalues plots (in dB)
    if 1
      % 1) Before calibrarion
      % ----------------------
      [V_uncalib, D_uncalib] = eig(Rxx_uncalib);
      [D_uncalib, D_uncalib_idxs] = sort(real(diag(D_uncalib)),'descend');
      V_uncalib = V_uncalib(:,D_uncalib_idxs);
      
      figure(9994);clf
      subplot(121)
      stem(10*log10(D_uncalib./max(D_uncalib)),'fill','b','LineWidth',2)
      xlabel('Eigenvalue index')
      ylabel('Eigenvalue (dB)')
      title('Before calib.')
      xlim([1 Nc])
      grid on
      grid minor
      
      % 2) After calibrarion
      % ----------------------
      [V_calib, D_calib] = eig(Rxx_calib);
      [D_calib, D_calib_idxs] = sort(real(diag(D_calib)),'descend');
      V_calib = V_calib(:,D_calib_idxs);
      
      figure(9994);
      subplot(122)
      stem(10*log10(D_calib./max(D_calib)),'fill','b','LineWidth',2)
      xlabel('Eigenvalue index')
      %    ylabel('Eigenvalue (dB)')
      title('After calib.')
      xlim([1 Nc])
      grid on
      grid minor
      
      % 3) Before calibration with forward-backward averaging
      % ------------------------------------------------------
      [V_uncalib_fb, D_uncalib_fb] = eig(Rxx_uncalib_fb);
      [D_uncalib_fb, D_uncalib_fb_idxs] = sort(real(diag(D_uncalib_fb)),'descend');
      V_uncalib_fb = V_uncalib_fb(:,D_uncalib_fb_idxs);
      
      figure(9995);clf
      subplot(121)
      stem(10*log10(D_uncalib_fb./max(D_uncalib_fb)),'fill','b','LineWidth',2)
      xlabel('Eigenvalue index')
      ylabel('Eigenvalue (dB)')
      title('Before calib. & FB')
      xlim([1 Nc])
      grid on
      grid minor
      
      % 4) After calibration with forward-backward averaging
      % ------------------------------------------------------
      [V_calib_fb, D_calib_fb] = eig(Rxx_calib_fb);
      [D_calib_fb, D_calib_fb_idxs] = sort(real(diag(D_calib_fb)),'descend');
      V_calib_fb = V_calib_fb(:,D_calib_fb_idxs);
      
      figure(9995);
      subplot(122)
      stem(10*log10(D_calib_fb./max(D_calib_fb)),'fill','b','LineWidth',2)
      xlabel('Eigenvalue index')
      %    ylabel('Eigenvalue (dB)')
      title('After calib. & FB')
      xlim([1 Nc])
      grid on
      grid minor
      
      % 5) Ideal case
      % -------------
      [V_ideal, D_ideal] = eig(Rxx_ideal_all{roll_angle_idx});
      [D_ideal, D_ideal_idxs] = sort(real(diag(D_ideal)),'descend');
      V_ideal = V_ideal(:,D_ideal_idxs);
    end
%     return
    %% Plot eigenvectors spectrum (or MUSIC pseudospectrum)
    if 1
      theta = linspace(-90,90,2048)*pi/180;
      
      % Steering vectors before calibration (i.e. uncalibrated array)
      sv_params = [];
      sv_params.src.y_pc = ideal_ypc;%uncalib_ypc;
      sv_params.src.z_pc = ideal_zpc;%uncalib_zpc;
      sv_params.src.fc   = fc;
      
      sv_params.extra_error_params = extra_error_params;
      sv_params.calib_params       = [];
      
      A = @(theta,param)steering_mtx(theta,sv_params);
      SVs_uncalib = A(theta,sv_params)./sqrt(Nc);
      
      % Steering vectors after calibration
      sv_params = [];
      sv_params.src.y_pc = uncalib_ypc; %calib_ypc; 
      sv_params.src.z_pc = uncalib_zpc; %calib_zpc; 
      sv_params.src.fc   = fc;
      
      sv_params.extra_error_params = [];
      sv_params.calib_params       = calib_params;
      
      A = @(theta,param)steering_mtx(theta,sv_params);
      SVs_calib = A(theta,sv_params)./sqrt(Nc);
      
      % Steering vectors of the ideal array (no errors case)
      sv_params = [];
%       rotation_mtx = [cos(-roll_angle), -sin(-roll_angle);sin(-roll_angle), cos(-roll_angle)];
%       rotated_pc_ideal = rotation_mtx*[ideal_ypc ideal_zpc].';
%       sv_params.src.y_pc = rotated_pc_ideal(1,:).';
%       sv_params.src.z_pc = rotated_pc_ideal(2,:).';
      sv_params.src.y_pc = ideal_ypc;
      sv_params.src.z_pc = ideal_zpc;
      sv_params.src.fc   = fc;
      
      sv_params.extra_error_params = [];
      sv_params.calib_params       = [];
      
      A = @(theta,param)steering_mtx(theta,sv_params);
      SVs_ideal = A(theta,sv_params)./sqrt(Nc);
      
      % Calculate the spectrum. FB averaging doesn't change the eigenvectors
      % of the DCM, but it does change the eigenvalues.
      N            = N;
      N_calib      = N;%1;
      N_uncalib_fb = N;%1;
      N_calib_fb   = N;%1;
      
      % 1) Before calibration
      % ----------------------
      f_uncalib = 1./sum(abs(V_uncalib(:,N+1:end)' * SVs_uncalib).^2,1);
      f_uncalib = 10*log10(f_uncalib./max(abs(f_uncalib)));
      
      % 2) After calibration
      % ----------------------
      % Either calibrate SVs or measurements, BUT NOT BOTH. Calibrating the
      % measurements correct the eigenvectors, and thus the DOA estimation, and
      % this is what you need to do here.
      f_calib = 1./sum(abs(V_calib(:,N_calib+1:end)' * SVs_uncalib).^2,1);
      %   f_calib = 1./sum(abs(V_uncalib(:,N_calib+1:end)' * SVs_calib).^2,1);
      f_calib = 10*log10(f_calib./max(abs(f_calib)));
      
      % 3) Before calibration with forward-backward averaging
      % ------------------------------------------------------
      f_uncalib_fb = 1./sum(abs(V_uncalib_fb(:,N_uncalib_fb+1:end)' * SVs_uncalib).^2,1);
      f_uncalib_fb = 10*log10(f_uncalib_fb./max(abs(f_uncalib_fb)));
      
      % 4) After calibration with forward-backward averaging
      % ------------------------------------------------------
      f_calib_fb = 1./sum(abs(V_calib_fb(:,N_calib_fb+1:end)' * SVs_uncalib).^2,1);
      f_calib_fb = 10*log10(f_calib_fb./max(abs(f_calib_fb)));
      
      % 5) Ideal case
      % --------------
      f_ideal = 1./sum(abs(V_ideal(:,N+1:end)' * SVs_ideal).^2,1);
      f_ideal = 10*log10(f_ideal./max(abs(f_ideal)));
      
      % DOA at maximum point in the spectrum
      [~ ,theta_uncalib_idx] = max(f_uncalib);
      theta_max_uncalib = theta(theta_uncalib_idx)*180/pi;
      
      [~, theta_max_f_calib_idx] = max(f_calib);
      theta_max_calib = theta(theta_max_f_calib_idx)*180/pi;
      
      [~, theta_uncalib_fb_idx] = max(f_uncalib_fb);
      theta_max_uncalib_fb = theta(theta_uncalib_fb_idx)*180/pi;
      
      [~, theta_calib_fb_idx] = max(f_calib_fb);
      theta_max_calib_fb = theta(theta_calib_fb_idx)*180/pi;
      
      sprintf('\nUncalib.: %2.2f deg. \n Calib.: %2.2f deg. \n Uncalib. with FB: %2.2f deg. \n Calib. with FB: %2.2f deg. \n',...
        theta_max_uncalib,theta_max_calib,theta_max_uncalib_fb,theta_max_calib_fb)
      
      figure(9996);clf
      subplot(131)
      plot(theta*180/pi,f_uncalib,'b')
      xlim(actual_doa*180/pi+[-60  60])
      xlabel('\theta ^\circ')
      ylabel('Pseudo-power (dB)')
      title('Before calib.')
      grid on
      grid minor
      
      subplot(132)
      plot(theta*180/pi,f_calib,'b')
      xlim(actual_doa*180/pi+[-60  60])
      xlabel('\theta ^\circ')
      %     ylabel('Pseudo-power (dB)')
      title('After calib.')
      grid on
      grid minor
      
      subplot(133)
      plot(theta*180/pi,f_ideal,'b')
      xlim(actual_doa*180/pi+[-60  60])
      xlabel('\theta ^\circ')
%       ylabel('Pseudo-power (dB)')
      title('Ideal')
      grid on
      grid minor
      suptitle(sprintf('Eigenvectors pattern: Nsig=%2d',N))

      figure(9997);clf
      subplot(121)
      plot(theta*180/pi,f_uncalib_fb,'b')
      xlim(actual_doa*180/pi+[-60  60])
      xlabel('\theta ^\circ')
      ylabel('Pseudo-power (dB)')
      title('Before calib.&FB')
      grid on
      grid minor
      
      subplot(122)
      plot(theta*180/pi,f_calib_fb,'b')
      xlim(actual_doa*180/pi+[-60  60])
      xlabel('\theta ^\circ')
      %   ylabel('Pseudo-power (dB)')
      title('After calib.&FB')
      grid on
      grid minor
      suptitle(sprintf('Eigenvectors pattern: Nsig=%2d',N))
    end
    
    return
    
    %%
    Color  = {'b','r','k','c','m','g',[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560]};
    Marker = {'-*','-o','-+','-.','-s','-^','-<','->'};
    
    if 1     
      theta_RP = linspace(-90,90,2048)*pi/180;
      [~, DOA_idx] = min(abs(DOA-theta_RP)); % For main beam direction
      
      w = hanning(Nc);
      %     w = ones(Nc,1);
      
      % Plot 3dB radiation pattern for 2 cases
      % ---------------------------------------------------------------------
      % Case 1: Before array calibration
      sv_params.src.y_pc = uncalib_ypc;%y_pos{fn_idx}{doa_idx};
      sv_params.src.z_pc = uncalib_zpc;%z_pos{fn_idx}{doa_idx};   % z_pc(:,doa_idx);
      sv_params.src.fc   = fc;
      
      sv_params.extra_error_params = extra_error_params;
      sv_params.calib_params       = [];
      
      A = @(theta,param)steering_mtx(theta,sv_params);
      SVs_uncalib = A(theta_RP,sv_params)./sqrt(Nc);
      A0 = SVs_uncalib(:,DOA_idx); % Main-beam direction
      
      RP = abs((w.*A0)'*SVs_uncalib).^2;
      RP = RP./max(RP);
      RP_dB = 10*log10(RP);
      
      idxs_3dB  = find(RP_dB>=-3);
      RP_3dB_uncalib    = RP_dB(idxs_3dB);
      theta_3dB_uncalib = theta_RP(idxs_3dB);
      
      beam_doa_lims = [min(theta_3dB_uncalib) max(theta_3dB_uncalib)];
      
      RP_dB_uncalib = RP_dB;
      beam_doa_lims_vec(:,1) = beam_doa_lims;
      
      figure(100);clf
      subplot(121)
      plot(theta_3dB_uncalib*180/pi,RP_3dB_uncalib,'b')
      xlabel('\theta^\circ')
      ylabel('Power (dB)')
      title('Before calibration')
      xlim([min(theta_3dB_uncalib*180/pi)  max(theta_3dB_uncalib*180/pi)])
      %     xlim([beam_doa_lims(1) beam_doa_lims(2)]*180/pi)
      grid on
      
      % Case 2: After array calibration
      sv_params.src.y_pc = calib_ypc;
      sv_params.src.z_pc = calib_zpc;
      sv_params.src.fc   = fc;
      
      sv_params.extra_error_params = [];
      sv_params.calib_params       = [];
      
      A = @(theta,param)steering_mtx(theta,sv_params);
      SVs_calib = A(theta_RP,sv_params)./sqrt(Nc);
      A0 = SVs_calib(:,DOA_idx);
      
      RP = abs((w.*A0)'*SVs_calib).^2;
      RP = RP./max(RP);
      RP_dB = 10*log10(RP);
      
      idxs_3dB  = find(RP_dB>=-3);
      RP_3dB_calib    = RP_dB(idxs_3dB);
      theta_3dB_calib = theta_RP(idxs_3dB);
      
      beam_doa_lims = [min(theta_3dB_calib) max(theta_3dB_calib)];
      
      RP_dB_calib = RP_dB;
      beam_doa_lims_vec(:,2) = beam_doa_lims;
      
      figure(100);
      subplot(122)
      plot(theta_3dB_calib*180/pi,RP_3dB_calib,'b')
      xlabel('\theta^\circ')
      ylabel('Power (dB)')
      title('After calibration')
      xlim([min(theta_3dB_calib*180/pi)  max(theta_3dB_calib*180/pi)])
      %     xlim([min(theta_3dB*180/pi)  max(theta_3dB*180/pi)])
      %     xlim([beam_doa_lims(1) beam_doa_lims(2)]*180/pi)
      grid on
      suptitle(sprintf('Array 3dB radiation pattern of %2.2f deg. target',DOA*180/pi))
      
      close(figure(100))
      
      % Plot the three radiation patterns above in one plot
      figure(101);clf
      hold on
      plot(theta_RP*180/pi,RP_dB_calib,'b')
      plot(theta_RP*180/pi,RP_dB_uncalib,'r')
      %     xlim([beam_doa_lims_vec(1,1) beam_doa_lims_vec(2,1)]*180/pi)
      %     xlim([min(beam_doa_lims_vec(1,:)) max(beam_doa_lims_vec(2,:))]*180/pi)
      ylim([-3 0])
      xlabel('\theta^\circ')
      ylabel('Power (dB)')
      title(sprintf('Array 3dB radiation pattern of %2.2f deg. target',DOA*180/pi))
      grid on
      legend('After calibration','Before calibration','Location','best')
      
      % Plot sensors gain pattern
      % ---------------------------------------------------------------------
      % After array calibration (before calibration it was 1 or uniform)
      theta_mtx = repmat(theta_RP,[Nc 1]);
      error_g_p_mtx = repmat(calib_params.calib_g_p,[1 length(theta_RP)]);
      error_g_s_mtx = repmat(calib_params.calib_g_s,[1 length(theta_RP)]);
      error_g_offset_mtx = repmat(calib_params.calib_g_offset,[1 length(theta_RP)]);
      gain_dB = error_g_s_mtx.* (sin(theta_mtx) - sin(error_g_p_mtx)).^2 + error_g_offset_mtx; % In dB
      
      figure(102);clf
      %     subplot(121)
      hold on;
      for chan_idx = 1:Nc
        chan_resp = 10.^(gain_dB(chan_idx,:));
        %     chan_resp = SVs_calib(chan_idx,:);
        chan_gain = abs(chan_resp).^2;
        %     chan_gain = chan_gain./max(chan_gain);
        chan_gain_dB = 10*log10(chan_gain);
        if all(abs(chan_gain_dB)<1e-4)
          chan_gain_dB = zeros(size(chan_gain_dB));
        end
        plot(theta_RP*180/pi,chan_gain_dB,'Color',Color{chan_idx})
      end
      xlabel('\theta^\circ')
      ylabel('Power (dB)')
      title(sprintf('Estimated gain relative to sensor 1 (after calib.)'))
      xlim([-20  20])
      grid on
      legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
      
      % Plot sensors phase deviation
      % ---------------------------------------------------------------------
      % After array calibration (before calibration it was 0)
      phase_deg = (calib_params.calib_ypc*k*sin(theta_RP) - calib_params.calib_zpc*k*cos(theta_RP) + ...
        repmat(calib_params.calib_phase + ac_est_error(:,2)*k,[1 length(theta_RP)]))*180/pi;
      
      figure(103);clf
      hold on
      for chan_idx = 1:Nc
        % chan_phase = SV_phase_dev(chan_idx,:);
        %     chan_phase = angle(chan_resp)*180/pi;
        chan_phase = phase_deg(chan_idx,:);
        if all(abs(chan_phase)<1e-4)
          chan_phase = zeros(size(chan_phase));
        end
        plot(theta_RP*180/pi,chan_phase,'Color',Color{chan_idx})
      end
      xlabel('\theta^\circ')
      ylabel('Phase (deg.)')
      title(sprintf('Estimated phase deviation relative to sensor 1 (after calib.)'))
      xlim([-20  20])
      %   xlim([ceil(min(theta_3dB_calib)*180/pi)   floor(max(theta_3dB_calib)*180/pi)])
      grid on
      legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
      
      % Plot sensors phase pattern
      % ---------------------------------------------------------------------
      % Case 1: Before  calibration
      figure(104);clf
      subplot(211)
      hold on
      for chan_idx = 1:Nc
        chan_resp = SVs_uncalib(chan_idx,:);
        chan_phase = angle(chan_resp)*180/pi;
        plot(theta_RP*180/pi,chan_phase,'Color',Color{chan_idx})
      end
      %     xlabel('\theta^\circ')
      ylabel('Phase (deg.)')
      title('Before calibration')
%       xlim([ceil(min(theta_3dB_uncalib)*180/pi)   floor(max(theta_3dB_uncalib)*180/pi)])
       xlim([-10 +10])
      grid on
      legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
      
      % Case 2: After calibration
      figure(104);
      subplot(212)
      hold on
      for chan_idx = 1:Nc
        chan_resp = SVs_calib(chan_idx,:) ;
        chan_phase = angle(chan_resp)*180/pi;
        plot(theta_RP*180/pi,chan_phase,'Color',Color{chan_idx})
      end
      xlabel('\theta^\circ')
      ylabel('Phase (deg.)')
      title('After calibration')
%       xlim([ceil(min(theta_3dB_calib)*180/pi)   floor(max(theta_3dB_calib)*180/pi)])
       xlim(DOA*180/pi+[-10 +10])
      grid on
      %     legend('Ant 1','Ant 2', 'Ant 3', 'Ant 4', 'Ant 5', 'Ant 6','Ant 7','Location','best');
      suptitle('Phase pattern relative to sensor 1')
      
    end
    