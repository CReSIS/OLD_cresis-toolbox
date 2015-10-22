function csarp_simulator(param)
% CReSIS CSARP compatible radar simulator
%
% INPUTS:
% param.
%  radar = paramaters for the transmitted radar signal
%   radar.fc:   carrier frequency (Hz)
%   radar.fs:   sampling frequency (Hz)
%   radar.prf:  pulse repetition frequency (Hz)
%   radar.Tpri:    pulse period (s)
%   radar.beamwidth:  antenna beamwidth (rads)
%   radar.wfs{wf}.presums:    number of coherent averages
%   radar.wfs.f0:   linear chirp starting frequency (Hz)
%   radar.wfs.f1:   linear chirp ending frequency (Hz)
%   radar.wfs.BW:   linear chirp bandwidth (Hz)
%   radar.wfs.Tpd:  linear chirp pulse duration (s)
%   radar.wfs.t0:   sampling delay (s)
%   radar.chans:    number of receive channels
%   radar.rx_path:  adc-to-antenna ordering
%   radar.adc_bits: number of bits in the adc simulator
%   radar.adc_dynamic_gain: dynamic range of the adc simulator
%   radar.tx_weights: transmit antenna weighting vector
%   radar.tx_amp:   amplitude of transmit signal
%
%  target = parameters for the target(s) of interest
%   target.x:       along-track position of the target(s) (m)
%   target.z:       depth of target(s) (m)
%   target.Z0:      depth of imaged area (m)
%
%  platform = parameters for the radar platform
%   platform.h:     height of the platform above the ice surface
%   platform.v:     ground speed velocity of the platform
%   platform.lever: lever arm vector (or function)
%
%  errors = parameters for various error introduced to the data
%   errors.td:      time-delay errors vector
%   errors.phase:   phase errors vector
%   errors.amp:     amplitude error vector
%   errors.chan_equal: combined phase/amplitude phasor error
%   errors.noise_floor: additional noise power per channel
%
%  scene = parameters for the properties of the imaged area
%   scene.surface_angle: angle above horizontal for linear surface slope
%
% OUTPUTS:
%  gps = simulated gps data
%   gps.lat:        latitude
%   gps.lon:        longitude
%   gps.elev:       elevation
%   gps.roll:       roll
%   gps.pitch:      pitch
%   gps.heading:    heading
%   gps.gps_time:   gps time
%   gps.gps_source: gps source (novatel-simulator)
%   gps.surface:    surface pick (zeros)
%   gps.bottom:     bottom pick (zeros)
%
%  radar data is output in the MCoRDS binary format
%
%
% NOTE: Not yet set up to properly handle multiple waveforms
%
% Author: Logan Smith, Jilu Li
%
%

% ENVIRONMENTAL PARAMETERS
physical_constants
eta_ice=sqrt(er_ice);     %Ice refraction index
c_ice=c/eta_ice;          %Propagation speed in ice,m/sec

% antenna full-beamwidth (rads)
param.radar.beamwidth = param.radar.beamwidth_deg*pi/180;
% antenna half-beamwidth in ice (rads)
param.radar.beamwidth_ice = ...
  2*asin(sin(0.5*param.radar.beamwidth)/eta_ice);

% =========================================================================
% Scene Setup
% =========================================================================

% Calculate min and max footprint
Bmin=param.platform.h*tan(0.5*param.radar.beamwidth);     % minimum half-beamwidth
Bmax=Bmin+2*param.target.Z0*tan(0.5*param.radar.beamwidth_ice);  % maximum half-beamwidth
L=2*Bmax;                                     % synthic aperture length

wfs = [1]; % vector of waveforms to process


for wf_idx = 1:length(wfs)
  wf = wfs(wf_idx);
  % =======================================================================
  % Physical Domain Parameters
  % =======================================================================
  % 2.2)slow time domain parameters and arrays
  % sample spacing in aperture domain
  dx=param.platform.v*param.radar.Tpri*param.radar.wfs{wf}.presums;
  m=2*ceil(L/dx); % number of samples on aperture
  L=(m-1)*dx;     % make L exactly a mutilple integer of du
  
  x=dx*(-m/2:m/2-1); % synthetic aperture array
  
  % 2.3)Fast-time domain parameters and arrays
  Ts = param.radar.wfs{wf}.t0;                         % start time of sampling
  Tf=2*param.platform.h/cos(0.5*param.radar.beamwidth)/c+... % end time of sampling
    2*2*param.target.Z0/cos(0.5*param.radar.beamwidth_ice)/c_ice+...
    param.radar.wfs{wf}.Tpd;
  Tm=Tf-Ts;                         % fast-time interval of measurement
  dt=1/param.radar.fs;                    % fast-time sampling interval
  n=2*ceil((.5*Tm)/dt);             % number of time samples
  t=Ts+(0:n-1)*dt;                  % time vector for data acquisition
  depth=t*c_ice/2-param.platform.h/sqrt(er_ice); % depth vector
  % dw=2*pi/((n-1)*dt);               % frequency domain sampling
  % w=2*pi*param.radar.fc+dw*(-n/2:n/2-1);  % frequency array (centered at carrier)
  freq = linspace(param.radar.fs,2*param.radar.fs,n);
  w = 2*pi*freq;                    % frequency domain sampling
  
  % HACK to keep gps file the same for wider beamwidth
  %   param.radar.beamwidth_deg = 40; % antenna full-beamwidth (deg)
  %   param.radar.beamwidth = param.radar.beamwidth_deg*pi/180; % antenna full-beamwidth (rads)
  %   % antenna half-beamwidth in ice (rads)
  %   param.radar.beamwidth_ice = ...
  %     2*asin(sin(0.5*param.radar.beamwidth)/eta_ice);
  
  % =========================================================================
  % Create records, gps, and frames files
  % =========================================================================
  fprintf('Creating simulated records, gps, and frames files...\n')
  
  records.prf_computed = param.radar.prf;
  records.eprf_computed = ...
    floor(records.prf_computed/length(param.radar.wfs)/param.radar.wfs{wf}.presums);
  
  tot_secs = floor(m/records.eprf_computed); % total number of recorded seconds in profile
  utc_dt = 1/records.eprf_computed; % time between each utc sample
  extra = m-records.eprf_computed*tot_secs; % time after last full second
  utc_sod = [reshape(repmat(1:tot_secs,records.eprf_computed,1),1,m-extra) (tot_secs+1)*ones(1,extra)]-1;
  utc_frac = mod((0:m-1)*1/records.eprf_computed,1).*param.radar.fs/2;
  utc_leap_secs = utc_leap_seconds((datenum(2011,7,8,0,0,0)-datenum(1970,1,1,0,0,0))*86400);
  
  % Create fake lat/lon data and platform motion and navigation errors
  % [dist,az] = distance(69,-51,69.014,-51,WGS84.ellipsoid)
  sd_dx=0;
  sd_mo=0;
  dh_mo=sd_mo*randn(1,m);
  %   dh_mo = 50*sin(2*pi*(x-x(1))/(x(end)-x(1)));
  dx_mo=sd_dx*randn(1,m);  % motion error in aperture locations=zeros(n,m);

  if param.synth_gps
      dist = 0;
      end_pt = 69;
      while dist < L%+fk_overlap
          end_pt = end_pt+.0001;
          [dist,az] = distance(69,-51,end_pt,-51,WGS84.ellipsoid);
      end
      
      synth_gps.lat = linspace(69,end_pt,m);
      synth_gps.lon = linspace(-51,-51,m);
%       d_elev = linspace(0,(x(end)-x(1))*tand(5),m);
%       synth_gps.elev = zeros(1,m) + (-1)*(d_elev - mean(d_elev));
      d_elev = 0*50*sin(5*2*pi.*(1:m)/m);
      synth_gps.elev = zeros(1,m) + d_elev;
      synth_gps.roll = zeros(1,m);
      synth_gps.pitch = zeros(1,m);
      synth_gps.heading = zeros(1,m);
      ref = synth_gps;
      gps_param.lever_arm_fh = @lever_arm_mcords_simulator;
      gps_param.rx_path = param.radar.rx_path;
      gps_param.tx_weights = param.radar.tx_weights;
      gps = trajectory_with_leverarm(synth_gps,gps_param);     
      gps_start_time = (datenum(2011,7,8,0,0,0)-datenum(1970,1,1,0,0,0))*86400+utc_leap_secs;
      gps.gps_time = gps_start_time+utc_sod+(utc_frac/param.radar.fs*2);
      [x_ecef, y_ecef, z_ecef] = geodetic2ecef(gps.lat*pi/180, gps.lon*pi/180, gps.elev, WGS84.ellipsoid);
      lon_ref = mean(gps.lon);
      lat_ref = mean(gps.lat);
      elev_ref = mean(gps.elev);
      [x_enu,y_enu,z_enu] = ecef2lv(x_ecef, y_ecef, z_ecef, lat_ref*pi/180, lon_ref*pi/180, elev_ref, WGS84.ellipsoid);
      along_track =  [0 cumsum(sqrt(diff(x_enu).^2 + diff(y_enu).^2))];
%       along_track = geodetic_to_along_track(gps.lat,gps.lon,gps.elev,L);
      % Surface topography
      surf = along_track.*tan(param.scene.surface_angle);
      surf = surf - mean(surf);
      x = along_track - along_track(m/2+1);
      param.platform.h = z_enu + param.platform.h - surf;
      gps.surface = param.platform.h*2/c;
      gps.bottom = zeros(size(gps.surface));
      gps.gps_source = 'novatel-simulator';     
  else
      % load in real gps 
      gps_start_time = (datenum(2011,7,8,0,0,0)-datenum(1970,1,1,0,0,0))*86400+utc_leap_secs;
      gps.gps_time = gps_start_time+utc_sod+(utc_frac/param.radar.fs*2);
      synth_gps = load(param.real_gps_fn);
      gps_idx1 = 1293000;
      gps_idx2 = find(synth_gps.gps_time(gps_idx1:end)-synth_gps.gps_time(gps_idx1)>gps.gps_time(end)-gps.gps_time(1),1,'first');
      gps_idx2 = gps_idx2 + gps_idx1 -1;
      
      % find lat,lon,elev,roll,pitch and heading at gps.gps_time, this way
      % the eprf is the same as the simulated idea case
      synth_gps.gps_time = synth_gps.gps_time(gps_idx1:gps_idx2) - synth_gps.gps_time(gps_idx1) + gps.gps_time(1);
      synth_gps.lat = interp1(synth_gps.gps_time,synth_gps.lat(gps_idx1:gps_idx2),gps.gps_time);
      synth_gps.lon = interp1(synth_gps.gps_time,synth_gps.lon(gps_idx1:gps_idx2),gps.gps_time);
      synth_gps.elev = interp1(synth_gps.gps_time,synth_gps.elev(gps_idx1:gps_idx2),gps.gps_time);
      if 1 % modulate the elev variation to make it bigger
          synth_gps.elev = 0*mean(synth_gps.elev)*ones(size(synth_gps.elev)) + (max(synth_gps.elev)-min(synth_gps.elev)).*sin(5*2*pi.*(1:m)/m);
      end
      synth_gps.roll = interp1(synth_gps.gps_time,synth_gps.roll(gps_idx1:gps_idx2),gps.gps_time);
      synth_gps.pitch = interp1(synth_gps.gps_time,synth_gps.pitch(gps_idx1:gps_idx2),gps.gps_time);
      synth_gps.heading = interp1(synth_gps.gps_time,synth_gps.heading(gps_idx1:gps_idx2),gps.gps_time);
      synth_gps.gps_time = gps.gps_time;
      if param.debug_level > 1    % plot to determin gps_idx1, the start point where you want use the real gps and ins data
          figure(3); clf;plot(synth_gps.elev);
          figure(4); subplot(3,1,1);plot(synth_gps.pitch*180/pi);
          subplot(3,1,2);plot(synth_gps.roll*180/pi);subplot(3,1,3);plot(synth_gps.heading*180/pi);
      end
      ref = synth_gps;
      gps_param.lever_arm_fh = @lever_arm_mcords_simulator;
      gps_param.rx_path = param.radar.rx_path;
      gps_param.tx_weights = param.radar.tx_weights;
      gps = trajectory_with_leverarm(synth_gps,gps_param);     
      [x_ecef, y_ecef, z_ecef] = geodetic2ecef(gps.lat*pi/180, gps.lon*pi/180, gps.elev, WGS84.ellipsoid);
      lon_ref = mean(gps.lon);
      lat_ref = mean(gps.lat);
      elev_ref = mean(gps.elev);
      [x_enu,y_enu,z_enu] = ecef2lv(x_ecef, y_ecef, z_ecef, lat_ref*pi/180, lon_ref*pi/180, elev_ref, WGS84.ellipsoid);
      along_track =  [0 cumsum(sqrt(diff(x_enu).^2 + diff(y_enu).^2))];
      % Surface topography
      surf = along_track.*tan(param.scene.surface_angle);
      surf = surf - mean(surf);
      x = along_track - along_track(m/2+1);
      param.platform.h = z_enu + param.platform.h-surf;
      gps.surface = param.platform.h*2/c;
      gps.bottom = zeros(size(gps.surface));
      gps.gps_source = 'novatel-simulator';
  end
  save(fullfile(param.gps_path,param.gps_fn),'-struct','gps')
  output_along_track = along_track;
  
  % create flight coordinate system structure for motion compensation
  fcs_param.squint = [0; 0; -1];
  fcs_param.type = 0;
  fcs_param.Lsar = L;
  fcs = SAR_coord_system(fcs_param,gps,ref,along_track,output_along_track);
  
  %2.4)Echoed signals
  theta=zeros(1,m);
  tdtemp=zeros(1,m);
  theta0=15*pi/180;
  sigma=zeros(param.ntarget,m);
  sploss=zeros(param.ntarget,m);
  iceloss=zeros(param.ntarget,m);
  fn=zeros(param.ntarget,m);    % fn: reflectivity including effects from sigma, sploss and iceloss
  
  % transmit power in Watts (5 mW/-26 dB total or 1.25 mw/-38 dB per
  % channel)
  Pt = param.radar.tx_amp./param.radar.tx_weights./sum(param.radar.tx_weights);
  
  alpha = param.radar.wfs{wf}.BW/param.radar.wfs{wf}.Tpd; % chirp rate
  ref_length = floor(param.radar.wfs{wf}.Tpd*param.radar.fs);
  t_win = tukeywin(ref_length,param.radar.tuk_ratio);
  
  % =========================================================================
  % Create simulated data
  % =========================================================================
  % Motion compesation parameters
  mocomp_param.squint = [0,0,-1];
  mocomp_param.type = 0;
  mocomp_param.tx_weights = param.radar.tx_weights;
  fcs.type = fcs_param.type;
  fcs.squint = fcs_param.squint;
  
  % Force noise to be the same each time
  rs = RandStream('mt19937ar','seed',1);
  RandStream.setDefaultStream(rs);
  
  fprintf('Creating simulated data...\n')
  for rx=1:param.radar.chans
    mocomp_param.rx = rx;
    error_delays = motion_comp(fcs,gps,ref,along_track,along_track); % m
    s=zeros(n,m);
    for targ=1:param.ntarget
      % Nt is fast time dim
      % Nx is slow time dim
      % Ts = scalar first time sample of time gate/recording window (sec)
      % dt = sample spacing in fast-time (sec)
      
      % tuk_mask = Nt by Nx matrix of the tukey fast-time time-domain weights
      tuk_mask = zeros(size(s));
      for rline=1:m
        % theta = 1 by Nx, incidence angle in air (radians)
        % Ra = range in air (m)
        % Ri = range in ice (m)
        h(rline) = param.platform.h(rline) + 0*error_delays(rline); % platform height above surface
        d = param.target.z(targ)+surf(rline); % target depth below surface
        yn = param.target.x(targ)-(x(rline)+dx_mo(rline)); % along track offset between radar and target
        surf_angle = param.scene.surface_angle;
        [theta(rline),Ra,Ri] = refraction(h(rline),d,yn,surf_angle,eta_ice);
%         [theta(rline),Ra,Ri] = refraction(h(rline),d,yn,eta_ice);
        
        % tdtemp = 1 by Nx time delay to the target
        tdtemp(1,rline)=2*(Ra+sqrt(er_ice)*Ri)/c;           % two-way time delay to target
        
        % start_bin = scalar, first bin of target return
        start_bin = round((tdtemp(1,rline)-Ts)/dt);
        
        % update tuk_mask for this range line,target
        tuk_mask(start_bin:(start_bin+ref_length-1),rline) = t_win;
        
        if param.visible_surf
          t_surf(rline) = h(rline)*2/c;
          surf_idx(rline) = find(t >= t_surf(rline),1);
        end
        
        % Return power
        if param.validation_mode_en
          Pr(rx) = Pt(rx);
        else
          % sigma(targ,rline)=exp(-2.73*(theta(rline)/theta0));    % point target response
          % sploss(targ,rline)=1/(4*pi*Ra^2)^2;                    % spreading loss
          % iceloss(targ,rline)=10^(-20*2*Ri*0.001/10);            % ice propagation loss
          % fn(targ,rline)=sqrt(sigma(targ,rline)*sploss(targ,rline)*iceloss(targ,rline));
          lambda = (c/param.radar.fs);
          Pr(rx) = Pt(rx).*G.^2*lambda.^2.*sigma(targ,rline)./((4*pi).^3.*(Ra+Ri).^4);
        end
        am(rline)=sqrt(Pr(rx)*2*param.radar.Z0);  % return signal amplitude in peak voltage
      end
      % td = Nt by Nx, a new time axis for each rline that has time zero
      %   set to be the start of target return
      td=repmat(t(:),[1,m])-repmat(tdtemp,[n,1]);
      
      % param.radar.wfs{wf}.Tpd = scalar pulse duration (sec)
      % param.radar.wfs{wf}.f0 = start frequency (Hz)
      % alpha = chirp rate * pi (Hz/sec * pi)
      % sn = Nt by Nx matrix of current target energy
      sn = repmat(am,n,1) .* exp(1i*2*pi*param.radar.wfs{wf}.f0*td + 1i*pi*alpha*td.^2) ...
        .* tuk_mask .* (td >= 0 & td <= param.radar.wfs{wf}.Tpd &...
        repmat(abs(theta)<=0.5*param.radar.beamwidth,[n,1]));
      
      
      %       sn=sn.*exp(-1i*2*pi*param.radar.fc*repmat(t(:),[1,m]));     % Fast-time baseband conversion
      
      if param.visible_surf
%         % start_bin = scalar, first bin of target return
%         start_bin = round((t_surf-Ts)/dt);
%         
%         % update tuk_mask for this range line,target
%         tuk_mask(start_bin:(start_bin+ref_length-1),rline) = t_win;
        t_mat = repmat(t(:),[1,m]);
        t_surf_mat = t_mat - repmat(t_surf(:).',[n,1]);
%         surf = repmat(am,n,1) .* exp(1i*2*pi*param.radar.wfs{wf}.f0*t_mat + 1i*pi*alpha*t_mat.^2) ...
%           .* tuk_mask .* (t_surf_mat >= 0 & t_surf_mat <= param.radar.wfs{wf}.Tpd);
        surf_mat = repmat(1*am,n,1) .* exp(1i*2*pi*param.radar.wfs{wf}.f0*t_surf_mat + 1i*pi*alpha*t_surf_mat.^2) ...
          .* (t_surf_mat >= 0 & t_surf_mat <= param.radar.wfs{wf}.Tpd);
        sn = sn+surf_mat;
        %           for sline=1:m
        %             sn(surf_idx(sline),sline) = am(sline)*exp(1i*2*pi*param.radar.fc*t_surf(sline));
        %           end
      end
      
      % s = Nt by Nx data matrix, apply principle of superposition and just
      %   add in each target's response
      s = s + sn;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADD NOISE
    F_dB = 1.91; k = BoltzmannConst; T0 = 290;
    noise_BW = param.radar.fs;
    %   thermal_noise_pwr = 10*log10(k*T0*param.radar.wfs{wf}.BW)+F_dB; % in dBw
    thermal_noise_pwr = 10*log10(k*T0*noise_BW)+F_dB; % in dBw
    % The quantization noise floor for the ADC is 10 dBm above the resolution
    % of the ADC.  In this case, the ADC resolution is
    % 20*log10(2/2^15) = -84 dBm.  Add 10 dB to reach the quantization noise
    % floor in practice and another 10 dB to avoid any serious quantization
    % effects and the effective noise floor for the simulated ADC is -54 dBm.
    
    % Adjust noise power to ADC noise floor and the noise removed by
    % windowing in the processing
    ADC_noise_floor = 20*log10(2/2^(param.radar.adc_bits+1));
    windowed_noise_factor = param.radar.wfs{wf}.BW/param.radar.fs;
    %     windowed_noise_factor = noise_BW/param.radar.fs;
    
    rx_gain_dB = thermal_noise_pwr - ADC_noise_floor - 10 - 10*log10(2);% + 10*log10(windowed_noise_factor);
    noise_pwr_dB = thermal_noise_pwr-rx_gain_dB+param.errors.noise_floor(rx);
    
    % Add noise data
    noise_mat = sqrt(10^(noise_pwr_dB/10)*param.radar.Z0)*(randn(size(s)) + 1j*randn(size(s)));
    
    s = (s + noise_mat)/sqrt(2) + 0.*max(am(:))*randn(size(s));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Apply time delay
    s = ifft(fft(s).*exp(-1i.*repmat(w.',1,m).*param.errors.td(rx)));
    phase_adjust = 360-mod(2*pi*param.radar.fc*param.errors.td(rx)*180/pi,360); % degs
    s = s.*exp(-1j*phase_adjust/180*pi);
    
    s = s*param.radar.rx_gain(rx);
    
    % Apply channel imbalances
    s = s.*param.errors.chan_equal(rx);
    
    % Remove imaginary portion of the data in order to prepare for
    % ADC quantization
    s = real(s);
    
    % Detect clipped signal
    pos_clip_idxs = find(s(:) > param.radar.Vpp/2);
    neg_clip_idxs = find(s(:) < -param.radar.Vpp/2);
    
    if ~isempty(pos_clip_idxs) || ~isempty(neg_clip_idxs)
      warning('Signal is too large and has been clipped by ADC!!!')
      dbstop
    end
    s(pos_clip_idxs) = param.radar.Vpp/2;
    s(neg_clip_idxs) = -param.radar.Vpp/2;
    
    % Data quantization
    which_bits = max(log2(param.radar.wfs{wf}.presums)-2,0);
    quantization_factor = (2^14/param.radar.Vpp) * (param.radar.wfs{wf}.presums/2^which_bits);
    s = round(s*quantization_factor + 2^param.radar.adc_bits);
    
    if param.debug_level > 1
      figure(1); clf;
      imagesc(x,depth,20*log10(abs(s)/max(max(abs(s)))),[-60 0]);
      %     imagesc(x,t*1e6,20*log10(abs(s)/max(max(abs(s)))),[-60 0]);
      xlabel('Along-Track Position (m)')
      ylabel('Depth in Ice (m)')
      title('Simulated SAR signal')
    end
    if 1
      % =======================================================================
      % Convert Data to CReSIS 1U-DAQ .dat form
      % =======================================================================
      fprintf('Creating simulated binary file...\n')
      out_fn = fullfile(sprintf('%s/chan%d/',param.data_path,rx),sprintf(param.data_fn,rx));
      fprintf('%s\n',out_fn)
      sprintf(param.data_path,rx);
      if ~exist(sprintf('%s/chan%d/',param.data_path,rx),'dir')
        mkdir(sprintf('%s/chan%d/',param.data_path,rx))
      end
      save(fullfile(param.data_path,'simulator_param.mat'),'param')
      fid = fopen(out_fn,'w+');
      for rline=1:size(s,2)
        bit_count = 0;
        for wf=1:length(param.radar.wfs)
          % Data Format
          % Byte 0, uint32: 0xDEADBEEF
          fwrite(fid,hex2dec('DEADBEEF'),'uint32',0,'b');
          bit_count = bit_count + 32;
          % Byte 4, uint32: Version
          fwrite(fid,1.0,'uint32',0,'b');
          bit_count = bit_count + 32;
          % Byte 8, uint32: UTC seconds of day
          fwrite(fid,utc_sod(rline),'uint32',0,'b');
          bit_count = bit_count + 32;
          % Byte 12, uint32: UTC fraction
          fwrite(fid,utc_frac(rline),'uint32',0,'b');
          bit_count = bit_count + 32;
          % Byte 16, uint32: ERPI
          fwrite(fid,rline,'uint32',0,'b');
          bit_count = bit_count + 32;
          % Byte 20, uint32: Number of Waveforms
          fwrite(fid,1.0,'uint32',0,'b');
          bit_count = bit_count + 32;
          % Byte 24, uint32: Bit Shifts
          fwrite(fid,0,'uint32',0,'b');
          bit_count = bit_count + 32;
          % Byte 28, uint32: Decimate configuration
          fwrite(fid,0,'uint32',0,'b');
          bit_count = bit_count + 32;
          % Byte 32+(N-1)*8, uint32: Waveform settings
          % Reserved
          fwrite(fid,0,'ubit18',0,'b');
          % Number of samples in waveform
          fwrite(fid,n,'ubit14',0,'b');
          bit_count = bit_count + 32;
          % Byte 36+(N-1)*8, uint32: Waveform settings
          % Reserved
          fwrite(fid,0,'ubit3',0,'b');
          % Number of logical right bit shifts
          fwrite(fid,which_bits,'ubit5',0,'b');
          % Number of sample clocks between PRF trigger and recording
          fwrite(fid,(11e-6+param.radar.wfs{wf}.t0)*param.radar.fs,'ubit14',0,'b');
          % Number of presums (N-1)
          fwrite(fid,param.radar.wfs{wf}.presums-1,'ubit10',0,'b');
          bit_count = bit_count + 32;
          % Fill in blank bytes
          bit_jump = 160-bit_count/8;
          fwrite(fid,zeros(1,bit_jump/4),'uint32',0,'b');
          % Byte 160, M1*uint16: Data for waveform 1
          fwrite(fid,s(:,rline),'uint16',0,'b');
        end
      end
      fclose(fid);
    end
    clear s
  end
end

return











