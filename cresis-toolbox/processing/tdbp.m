function tdbp_data =  tdbp(param,data)
% tdbp_data = tdbp(param,data)
% time domain back projection SAR processor
% param = structure controlling the time domain processor
%   .c = speed of light in free space
%   .eta_ice = refraction index of air-ice interface
%   .refraction_method
%      = 1 linear approximation for refraction, faster and
%          valid for linear surface and small surface slope
%      = 2 resolve 4th order polynomial equation for refraction
%          valid for linear surface with surface slope
%      = 3 Newton iterative method, valid for nonlinear surface
%   .surfBins = range bin indexes of ice surface
%   .fc = transmit center frequency
%   .time = two-way propagation time at each range bin
%   .proc = structure setting the following processing parameters
%      .Nfft = length of fft of upsampling for interpolation
%      .skip_surf = 1, skip processing pixels above ice surface
%       indcluding surface pixels)
%                 = 0 or empty (default),not to skip surface pixels
%      .start_range_bin_above_surf = the range bin number of the first pixel
%       to be processed by time domain porcessor; it is counted from
%       the ice surface to the pixel
%      .start_range_bin = the range bin number of the first pixel
%       to be processed by time domain porcessor when skip processing
%       the surface and the pixels above; it is counted from the first
%       range bin. This is used for just processing the pixels around the
%       ice bottom.
%      .end_range_bin = the range bin number of the last pixel
%       to be processed by time domain porcessor; it is counted from
%       the first range bin
% data = the pulse compressed data of the chunk being processed
%      
% tdbp_data = output of the time domain processor. Unprocessed pixels have
%             the same values with the pulse compressed data.

% Author: Jilu Li
% See also: csarp_task.m


% constant parameters
c_inv = 1/param.c;
freq_shift = 1j*2*pi*param.fc;
sub_apt_shift_num = param.sub_aperture_steering(end) - param.sub_aperture_steering(1) + 1;
[n,m] = size(data);

% initilizatoin for output
tdbp_data = data(:,param.output_along_track_idxs);
surfBins_at_output = param.surfBins(param.output_along_track_idxs);
out_idxs = param.output_along_track_idxs;
out_idxs_step = out_idxs(2)-out_idxs(1);
m1 = length(param.output_along_track_idxs);

if isempty(param.refraction_method)
  param.refraction_method = 1;               % default value
end
if isempty(param.proc.skip_surf)
  param.proc.skip_surf = 0;                   % default value
end
if isempty(param.proc.start_range_bin_above_surf)
  param.proc.start_range_bin_above_surf = 10;  % default value
end
if param.proc.skip_surf
  if isempty(param.proc.start_range_bin)
    param.proc.start_range_bin = max(surfBins_at_output) + 5;   % default value
  end
else
  param.proc.start_range_bin = [];
end
if isempty(param.proc.end_range_bin)
  param.proc.end_range_bin = n;  % default value
end

% number of aperture indexes varies with depth
d_along_track_inv = 1/mean(diff(param.along_track));
apt_length = zeros(n,1);
for iz = 1: n
  apt_length(iz) = ceil(-1.1*tan(param.HbeamWidth)*param.pixel(3,iz,1)*d_along_track_inv);
end
apt_length(apt_length<0) = 0;

% the aperture index for the first output
aperture_idxs0 = out_idxs(1) - apt_length(n) : out_idxs(1) + apt_length(n);
if aperture_idxs0(1) < 1
  aperture_idxs0(aperture_idxs0<1) = [];
end
if aperture_idxs0(end) > m
  aperture_idxs0(aperture_idxs0>m) = [];
end

% upsampling data in range dimention
time_fine = linspace(param.time(1),param.time(end),param.proc.Nfft);
dt_fine_inv = 1/(time_fine(2)-time_fine(1));
signal = interpft(data(:,aperture_idxs0),param.proc.Nfft);   % data buffer for interpolation
along_track = param.along_track(aperture_idxs0);

for ix = 1:m1    % along track for loop
  out_idx = param.output_along_track_idxs(ix);
  ps = polyfit(along_track,param.surf(3,aperture_idxs0),1);
  surfAngle = atan(ps(1));
  
  % depth for loop
  if ~param.proc.skip_surf  %pixel above the ice surface
    for iz = surfBins_at_output(ix)-param.proc.start_range_bin_above_surf:surfBins_at_output(ix)
      % moving aperture with along-track location of output
      aperture_idxs = out_idx - apt_length(iz) : out_idx + apt_length(iz);
      if aperture_idxs(1) < 1
        aperture_idxs(aperture_idxs<1) = [];
      end
      if aperture_idxs(end) > m
        aperture_idxs(aperture_idxs>m) =[];
      end
      m2 = length(aperture_idxs);
      m3 = find(aperture_idxs0==aperture_idxs(1));
      if isempty(m3)
        m3 = 1;
      end
      
      signal_ip = [];
      t_ixiz = [];
      for iu = 1:m2  % aperture for loop
        if m3>size(signal,2)
          break
        end
        Ra = sqrt((param.pixel(1,iz,ix)- param.phase_center(1,aperture_idxs(iu)))^2 + ...
          (param.pixel(2,iz,ix)- param.phase_center(2,aperture_idxs(iu)))^2 +...
          (param.pixel(3,iz,ix)- param.phase_center(3,aperture_idxs(iu)))^2);
        theta_ixiz = asin((param.along_track(out_idx)- along_track(m3))/Ra);
        if abs(theta_ixiz)<=param.HbeamWidth
          t_ixiz_tmp = 2*Ra*c_inv;
          z_idx = (t_ixiz_tmp-time_fine(1))*dt_fine_inv+1;
          if z_idx > size(signal,1)
            m3 = m3 + 1;
            continue
          end
          z_idx1 = floor(z_idx);
          diz = z_idx-z_idx1;
          z_idx2 = z_idx1 + 1;
          signal_ip_tmp = signal(z_idx1,m3) + diz*(signal(z_idx2,m3)-signal(z_idx1,m3));
          t_ixiz = [t_ixiz,t_ixiz_tmp];
          signal_ip = [signal_ip,signal_ip_tmp];
        end
        m3 = m3 + 1;
      end
      su = signal_ip.*exp(freq_shift*(t_ixiz-param.time(iz)));
      apt_len = length(su);
      if apt_len==0
        continue
      else
        % calculate subaperture idxs
        if sub_apt_shift_num>1 & apt_len>2*sub_apt_shift_num 
          sub_apt_len = 2*ceil(floor(apt_len/sub_apt_shift_num)/2)-1;
          sub_apt_len_half = floor(0.5*sub_apt_len);
          sub_apt_center = floor(param.sub_aperture_steering(param.proc.sub_apt_idx)*sub_apt_len);
          sub_apt_idxs = [sub_apt_center - sub_apt_len_half : sub_apt_center + sub_apt_len_half] + ceil(0.5*apt_len);        
          tdbp_data(iz,ix)=sum(param.st_wind(sub_apt_len)'.*su(sub_apt_idxs))/sqrt(sum(param.st_wind(sub_apt_len)));
        else
          tdbp_data(iz,ix)=sum(param.st_wind(apt_len)'.*su)/sqrt(sum(param.st_wind(apt_len)));
        end
      end
    end
  end
  % pixel in ice
  if isempty(param.proc.start_range_bin)
    start_iz = surfBins_at_output(ix)+1;
  else
    start_iz = param.proc.start_range_bin;
  end
  for iz = start_iz:param.proc.end_range_bin
    % moving aperture with along-track location of output
    aperture_idxs = out_idx - apt_length(iz) : out_idx + apt_length(iz);
    if aperture_idxs(1) < 1
      aperture_idxs(aperture_idxs<1) = [];
    end
    if aperture_idxs(end) > m
      aperture_idxs(aperture_idxs>m) =[];
    end
    m2 = length(aperture_idxs);
    m3 = find(aperture_idxs0==aperture_idxs(1));
    if isempty(m3)
      m3 = 1;
    end
    m3_start = m3;
    iu_step = floor(m2/30);
    theta_ixiz = [];
    t_ixiz = [];
    iu_good_idxs = [];
    for iu = 1:iu_step:m2  % aperture for loop
      if m3>size(signal,2)
        break
      end
      if param.refraction_method == 1      % linear surface with zero or very small surface slope
        if param.along_track(out_idx)-along_track(m3) == 0
          along_track_i = along_track(m3);
          surf_i = param.surf(3,aperture_idxs(iu));
        else
          k_p2pix = (param.pixel(3,iz,ix)-param.phase_center(3,aperture_idxs(iu)))/...
            (param.along_track(out_idx)-along_track(m3));
          A = [ps(1) -1;k_p2pix,-1];
          B = [-ps(2);k_p2pix*along_track(m3)-param.phase_center(3,aperture_idxs(iu))];
          [xx] = linsolve(A,B);
          along_track_i = param.along_track(out_idx)+ (xx(1)-param.along_track(out_idx))/param.eta_ice;
          surf_i = xx(2);
        end
        Ra = sqrt((param.phase_center(3,aperture_idxs(iu))-surf_i)^2 + (along_track(m3)-along_track_i)^2);
        Ri = sqrt((param.pixel(3,iz,ix)-surf_i)^2 + (param.along_track(out_idx)-along_track_i)^2);
        theta_tmp = asin((along_track_i-along_track(m3))/Ra);
      elseif param.refraction_method == 2  % linear surface with surface slope
                [theta_tmp,Ra,Ri,err]=refraction(param.h(aperture_idxs(iu)),...
                  param.surf(3,aperture_idxs(iu)) - param.pixel(3,iz,ix),...
                  param.along_track(out_idx)- along_track(m3),surfAngle,param.eta_ice);       
      elseif param.refraction_method == 3  % nolinear surface, Newton iteration
      end
      if ~isempty(theta_tmp)
        iu_good_idxs = [iu_good_idxs,iu];
        theta_ixiz = [theta_ixiz,theta_tmp];
        t_ixiz = [t_ixiz,2*(Ra+param.eta_ice*Ri)*c_inv];
      end
      m3 = m3 + iu_step;
    end
    theta_ixiz = interp1(iu_good_idxs,theta_ixiz,[1:m2],'spline','extrap');
    t_ixiz = interp1(iu_good_idxs,t_ixiz,[1:m2],'spline','extrap');
    
    signal_ip = [];
    t_delay = [];
    m3 = m3_start;
    for iu = 1:m2
      if m3>size(signal,2)
        break
      end
      if abs(theta_ixiz(iu))+sign(theta_ixiz(iu))*surfAngle<=param.HbeamWidth
        z_idx = (t_ixiz(iu)-time_fine(1))*dt_fine_inv+1;
        if z_idx >= size(signal,1)
          m3 = m3 + 1;
          continue 
        end
        z_idx1 = floor(z_idx);
        diz = z_idx-z_idx1;
        z_idx2 = z_idx1 + 1;
        signal_ip_tmp = signal(z_idx1,m3) + diz*(signal(z_idx2,m3)-signal(z_idx1,m3));
        t_delay = [t_delay t_ixiz(iu)];
        signal_ip = [signal_ip,signal_ip_tmp];
      end
      m3 = m3 + 1;
    end
    if isempty(signal_ip)
      continue
    else
      su = signal_ip.*exp(freq_shift*(t_delay-param.time(iz)));
      apt_len = length(su);
      % calculate subaperture
      if sub_apt_shift_num>1 & apt_len>2*sub_apt_shift_num 
          sub_apt_len = 2*ceil(floor(apt_len/sub_apt_shift_num)/2)-1;
          sub_apt_len_half = floor(0.5*sub_apt_len);
          sub_apt_center = floor(param.sub_aperture_steering(param.proc.sub_apt_idx)*sub_apt_len);
          sub_apt_idxs = [sub_apt_center - sub_apt_len_half : sub_apt_center + sub_apt_len_half] + ceil(0.5*apt_len);        
          tdbp_data(iz,ix)=sum(param.st_wind(sub_apt_len)'.*su(sub_apt_idxs))/sqrt(sum(param.st_wind(sub_apt_len)));
      else
        tdbp_data(iz,ix)=sum(param.st_wind(apt_len)'.*su)/sqrt(sum(param.st_wind(apt_len)));
      end
    end
  end
  
  % update the interpolation buffer as the aperture moves
  update_idxs = (1:out_idxs_step) + aperture_idxs0(end);
  if update_idxs(end)>m
    extra_idxs = out_idxs_step-(update_idxs(end)-m)+1:out_idxs_step;
    update_idxs(extra_idxs) = [];
  end
  signal_update = interpft(data(:,update_idxs),param.proc.Nfft);
  signal(:,1:out_idxs_step) = [];
  signal = [signal,signal_update];
  aperture_idxs0 = aperture_idxs0 + out_idxs_step;
  if aperture_idxs0(1) < 1
    aperture_idxs0(aperture_idxs0<1) = [];
  end
  if aperture_idxs0(end) > m
    aperture_idxs0(aperture_idxs0>m) = [];
  end
  along_track = param.along_track(aperture_idxs0);
end
return;
