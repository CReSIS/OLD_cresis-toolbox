function [theta,sv] = array_proc_sv(Nsv, fc, yAnt, zAnt, roll, LUT, rx_paths)
% [theta,sv] = array_proc_sv(Nsv, fc, yAnt, zAnt, roll, LUT, rx_paths)
%
% Function specified in sv_fh column of array worksheet of param
% spreadsheet for use with array_proc.
%
% Nsv: Several Options:
%   Nsv is a cell array:
%     Nsv{1}: 'theta'
%     Nsv{2}: theta vector (radians) of angles at which steering vectors
%       will be made.
%   Nsv is numeric scalar:
%     Contains the number of uniformily sampled in wavenumber space
%     steering vectors will be created.
% fc: center frequency (Hz)
%   To adjust the dielectric of the steering vectors, pass in the effective
%   center frequency after considering the dielectric. For example:
%     fc_effective = fc_actual*sqrt(relative_dielectric)
% yAnt,zAnt: phase centers from lever arm (note this includes tx and
%   rx positions, so nyquist sampling is quarter wavelength for the
%   phase centers because of two way propagation... i.e. k = 4*pi*fc/c)
%   These are Ny by 1 vectors
% roll: roll information (aircraft attitude)
% LUT: measured steering vectors lookup table
% rx_paths: active receive adcs
%
% sv = steering vector of size Ny by Nsv
% theta = incidence angle vector of size 1 by Nsv, defined by atan2(ky,kz)
%
% y is increasing to the right
% z is increasing downwards
% sv = sqrt(1/length(yAnt)) * exp(1i*(-zAnt*kz + yAnt*ky));
% A positive ky with a positive y (y points to the left) implies a positive
% phase for the steering vector. This means the measurement is closer to
% the target the bigger y gets (i.e. the more to the left you go)
% and therefore positive ky implies a target from the left. Positive ky
% corresponds to positive theta.
% kz is always positive. A positive z is always moving away from the target.
%
% Author: John Paden

% Decide ideal or measured steering vectors generation
ni = nargin;
% Creation of linear steering vector for 2D arbitrary array
c = 2.997924580003452e+08; % physical_constants too slow
% Wavenumber for two way propagation
k = 4*pi*fc/c;

if iscell(Nsv)
  if strcmpi(Nsv{1},'theta')
    theta = Nsv{2};
    
    % shape doas into row vector
    theta = theta(:).';
    
    if ni == 4        % Ideal generation
      ky = k*sin(theta);
      kz = k*cos(theta);
      
    elseif ni > 4     % Measured generation
      
      if 1
        theta_lut = theta - roll;
        sv_corr = (interp1(LUT.bins,LUT.sv_real,theta_lut,'linear','extrap') + 1i*interp1(LUT.bins,LUT.sv_imag,theta_lut,'linear','extrap')).';
        
        ky = k*sin(theta);
        kz = k*cos(theta);
        sv = sqrt(1/length(yAnt)) * exp(1i*(-zAnt*kz + yAnt*ky));
        sv = sv .* sv_corr;
        
      else
        theta = theta + roll;
        ky = k*sin(theta);
        kz = k*cos(theta);
        range = -45:1:45;   % Default sv correction table range
        theta_deg = round(theta*180/pi);
        sv = zeros(size(yAnt,1),length(theta_deg));
        if length(theta_deg) > length(yAnt)
          angle_all = zeros(size(sv));
          %         angle_all(:,1) = -zAnt*kz(1) + yAnt*ky(1);
          %         angle_all(:,end) = -zAnt*kz(end) + yAnt*ky(end);
          sv(:,1) = sqrt(1/length(yAnt)) * exp(1i*(-zAnt*kz(1) + yAnt*ky(1)));
          sv(:,end) = sqrt(1/length(yAnt)) * exp(1i*(-zAnt*kz(end) + yAnt*ky(end)));
          
          idx1 = find(theta_deg >= range(1),1,'first');
          idx2 = find(theta_deg <= range(end),1,'last');
          
          [~,~,roll_idx] = intersect(theta_deg(idx1:idx2),range);
          %         tmp = zeros(length(yAnt),length(roll_idx));
          tmp = conj(LUT.sv_table(rx_paths,roll_idx));
          %         angle_all(:,idx1:idx2) = angle(conj(LUT.sv_table(rx_paths,roll_idx)));
          sv(:,idx1:idx2) = sqrt(1/length(yAnt)) * abs(tmp) .* exp(1i * angle(tmp));
          
          %         for ant_idx = 1:length(yAnt)
          angle_all = interp1([theta_deg(1) range theta_deg(end)], [angle(sv(:,1)) angle(sv(:,idx1:idx2)) angle(sv(:,end))].', theta_deg(1):theta_deg(end)).';
          %           angle_all = interp1([theta_deg(1) range theta_deg(end)], [angle_all(:,1) angle_all(:,idx1:idx2) angle_all(:,end)].', theta_deg(1):theta_deg(end)).';
          %           mag_tmp = interp1([theta_deg(1) range theta_deg(end)], [abs(sv(ant_idx,1)) abs(sv(ant_idx,idx1:idx2)) abs(sv(ant_idx,end))], theta_deg(1):theta_deg(end));
          sv = sqrt(1/length(yAnt)) .* exp(1i * angle_all);
          %         end
          
        else
          % Individual sv generation, no use in array_proc, just for
          % integrity
          for theta_idx = 1:length(theta_deg)
            tmp = zeros(length(yAnt),1);
            if (theta_deg(theta_idx) >= range(1) & theta_deg(theta_idx) <= range(end))
              roll_idx = find(theta_deg(theta_idx) == range);
              tmp = conj(LUT.sv_table(:,roll_idx));
              sv(:,theta_idx) = sqrt(1/length(yAnt)) * abs(tmp) .* exp(1i * angle(tmp));
            else
              sv(:,theta_idx) = sqrt(1/length(yAnt)) * exp(1i*(-zAnt*kz(theta_idx) + yAnt*ky(theta_idx)));
            end
          end
          
        end
      end
      
      return;
    end
  end
  
else
  
  % Choose equally spaced y-dimension (cross-track) wavenumbers
  dNsv = 2/Nsv;
  ky = dNsv *k* [0 : floor((Nsv-1)/2), -floor(Nsv/2) : -1];
  
  % Determine z-dimension (elevation) dimension wavenumbers for each ky
  kz = sqrt(k^2 - ky.^2);
  
  % Calculate the angle of arrival for each ky
  theta = atan2(ky,kz);
  %
  %   if ni > 4     % Measured generation
  %       [theta,sort_idx] = sort(theta,'ascend');
  %       theta = theta + roll;
  %       ky = k*sin(theta);
  %       kz = k*cos(theta);
  %       range = -45:1:45;   % Default sv correction table range
  %       theta_deg = round(theta*180/pi);
  %       sv = zeros(size(yAnt,1),length(theta_deg));
  %
  %       angle_all = zeros(size(sv));
  %       angle_all(:,1) = -zAnt*kz(1) + yAnt*ky(1);
  %       angle_all(:,end) = -zAnt*kz(end) + yAnt*ky(end);
  %       idx1 = find(theta_deg >= range(1),1,'first');
  %       idx2 = find(theta_deg <= range(end),1,'last');
  %
  %       [~,~,roll_idx] = intersect(theta_deg(idx1:idx2),range);
  % %         tmp = zeros(length(yAnt),length(roll_idx));
  % %         tmp = conj(LUT.sv_table(rx_paths,roll_idx));
  %       angle_all(:,idx1:idx2) = interp1(range, angle(conj(LUT.sv_table(rx_paths,roll_idx))).',theta_deg(idx1:idx2)).'
  % %       angle_all(:,idx1:idx2) = angle(conj(LUT.sv_table(rx_paths,roll_idx)));
  % %         sv(:,idx1:idx2) = sqrt(1/length(yAnt)) * abs(tmp) .* exp(1i * angle(tmp));
  %
  % %         for ant_idx = 1:length(yAnt)
  %       angle_all = interp1([theta_deg(1) theta_deg(idx1:idx2) theta_deg(end)], [angle_all(:,1) angle_all(:,idx1:idx2) angle_all(:,end)].', theta_deg(1):theta_deg(end)).';
  % %           mag_tmp = interp1([theta_deg(1) range theta_deg(end)], [abs(sv(ant_idx,1)) abs(sv(ant_idx,idx1:idx2)) abs(sv(ant_idx,end))], theta_deg(1):theta_deg(end));
  %       sv = sqrt(1/length(yAnt)) .* exp(1i * angle_all);
  %         return;
  %   end
  
end
if nargout < 2
  return;
end

% Take the outer product of the antenna positions with the trig(theta)
% to create 2D matrix. Normalize the steering vector lengths.
sv = sqrt(1/length(yAnt)) * exp(1i*(-zAnt*kz + yAnt*ky));
% Equivalent: sv = sqrt(1/length(yAnt)) * exp(1i*k*(-zAnt*cos(theta) + yAnt*sin(theta)));
