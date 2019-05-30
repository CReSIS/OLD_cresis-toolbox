function [actual_doa, actual_doa_len] = actual_number_of_targets(param)
% Function for determining the actual number of targets. Used in the 2D
% simulator.
%
% This scripts is mainly from the crosstrack_rmse.m script.
%
% Mohanad

physical_constants;

if isfield(param,'optimal_test') ||  isfield(param,'suboptimal_test')
  % This is Sravya's idea. She assumes that targets lying between two
  % range-bins belong to the first/preceding range-bin. It works with MOE
  % training that shee created.
  rlines    = param.rlines;
  rbins     = param.rbins;
  rline_rng = param.rline_rng;
  dline     = param.dline;
  
%   param.src.y_pc = param.phase_center(2,:);
%   param.src.z_pc = param.phase_center(3,:);
  
  surf_model   = param.surf_model;
  
  Nx = length(rlines);
  Nt = length(rbins);
  Nc = length(param.src.y_pc);
  dr = c/(2*param.BW);
  dt = 1/param.fs;
  
  time = param.t0 + dt*(0:Nt-1).';
  R_bins_values = time*c/2;
  R_shell_values = [R_bins_values;R_bins_values(end)+dr];
  
  for snapshot =1:Nx;
    chan =(Nc+1)/2;
        
    for target = 1:size(surf_model.z,1)
      % Determine the range from the receiver to the target
      %CHECK ERROR..............
      Rvec = [surf_model.y(target) - param.src.y_pc(chan);
        surf_model.z(target,snapshot) - param.src.z_pc(chan)];
      R_targets(target) = norm(Rvec,2);
    end
    
    clear Rvec
    
%     surf_model.R_targets = R_targets;
    doa_target = atan(surf_model.y./surf_model.z(:,snapshot))*180/pi;
    
    for dR_idx = 1:Nt
      bin_sort_idxs =  find ( R_targets < R_shell_values(dR_idx +1) & ...
        R_targets >= R_shell_values(dR_idx));
      R_targets(bin_sort_idxs) = dR_idx; % target sorted into this bin number
      sources_true(dR_idx) = length(bin_sort_idxs);
      doa_target_bins{snapshot}(dR_idx,1:length(bin_sort_idxs)) =  doa_target(bin_sort_idxs);
    end
    
    clear bin_sort_idxs
    
    sources_true_temp(:,snapshot) = sources_true.';
    clear sources_true
    
  end
  
  %SOURCES_TRUE for all valid/active range-lines
  records_idx = [max(rline_rng)+1:dline: Nx-max(rline_rng)];
  actual_doa_len = sources_true_temp(:,records_idx);
  
  if isfield(param,'SS') && ~isempty(param.SS) && param.SS==0
    actual_doa_len = actual_doa_len./param.dist_target_factor;
  end
  
  actual_doa = doa_target_bins{records_idx};
  
else
  % This is John's idea. He assumes that targets lying between two
  % consequtive range-bins belong to both range-bins. It makes more sense
  % than the previous one.
  rlines = param.rlines;
  rbins = param.rbins;
  Nx = length(rlines);
  Nt = length(rbins);
  dr = c/(2*param.BW);
  dt = 1/param.fs;
%   phase_center = -param.phase_center;
  
  %% Determine the number of sources and their DOA for each image pixel
  time = param.t0 + dt*(0:Nt-1).';
  range = time*c/2;
  actual_doa = cell(Nt,Nx);
  
  for x_idx = 1:Nx
    line = rlines(x_idx);
    for y_idx = 1:length(param.surf_model.y)
      % Determine range to the target
      R = sqrt(param.surf_model.y(y_idx).^2 + param.surf_model.z(y_idx,line).^2);
      %     src_doa = atan2(surf_model.y(y_idx),surf_model.z(y_idx,line));
      src_doa = atan(param.surf_model.y(y_idx)/abs(param.surf_model.z(y_idx,line)));
%       rbins = find(abs(R-range) < dr/2*2);
      [~,min_idx] = min(abs(R-range));
      rbins = min_idx;
%       if length(rbins) == 2
%         [~,min_idx] = min(abs(R-range));
%         rbins = min_idx;
%       end
      
      for rbin = rbins(:).'
        actual_doa{rbin,x_idx} = [actual_doa{rbin,x_idx} src_doa];
      end
    end
  end
    
  %% Group and filter DOAs that are close together
  
  % Convert doa from theta to ky
  k = 4*pi*param.fc/c;
  actual_doa = cellfun(@(theta) k*sin(theta),actual_doa,'UniformOutput',false);
  
  % Determine dky resolution to use
  My = 4;
  
  % phase_center = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
  Y = max(param.src.y_pc) - min(param.src.y_pc);
  dky = 2*pi / Y / My;
  
  for x_idx = 1:Nx
    for rbin = 1:Nt
      
      actual_doa{rbin,x_idx} = sort(actual_doa{rbin,x_idx});
      doa_idx = 1;
      while doa_idx <= length(actual_doa{rbin,x_idx})
        % Collect all doa's within the resolution cell and average into one
        % doa bin
        doa_bin_edge = actual_doa{rbin,x_idx}(doa_idx) + dky;
        doa_end_idx = doa_idx;
        while doa_end_idx+1 <= length(actual_doa{rbin,x_idx}) ...
            && actual_doa{rbin,x_idx}(doa_end_idx+1) < doa_bin_edge
          doa_end_idx = doa_end_idx + 1;
        end
        actual_doa{rbin,x_idx} ...
          = [actual_doa{rbin,x_idx}(1:doa_idx-1), ...
          mean(actual_doa{rbin,x_idx}(doa_idx:doa_end_idx)), ...
          actual_doa{rbin,x_idx}(doa_end_idx+1:end)];
        
        % Skip to the next doa index that was not binned
        doa_idx = doa_end_idx + 1;
      end
      
      % Ignore DoAs outside the 3dB antenna beampattern
      %     actual_doa{rbin,x_idx}(actual_doa{rbin,x_idx} < param.plot_doa_lims(1)) = [];
      %     actual_doa{rbin,x_idx}(actual_doa{rbin,x_idx} > param.plot_doa_lims(2)) = [];
    end
  end
  
  
  % Convert doa from ky to theta
  actual_doa = cellfun(@(ky) asin(ky/k),actual_doa,'UniformOutput',false);
  
  % actual_doa_len: Number of targets in each pixel
  actual_doa_len = cellfun(@length,actual_doa);
  
end
end