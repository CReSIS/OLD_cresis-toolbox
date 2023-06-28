function [Data, DCM, imp_resp, DCM_fd] = doa_wideband_data(param)
% [Data, DCM, imp_resp, DCM_fd] = doa_wideband_data(param)
%
% Funcion to simulate multichannel array wideband data.
%
% Usually called from sim.doa script
% Inputs:
% param
%       .src
%           .Nsnap = scalar value where user specifies number of snapshots.
%            Used to compute number of fast time bins, Nt, needed in
%            simulation to obtain the specified number of snapshots.  When
%            the widening factor is 0 (narrowband assumption), Nsnap = Nt.
%            When W > 0, Nsnap = W*Nt.
%           .SNR  = 1xQ vector specifying SNR per channel of each source.
%           .DOAs = 1xQ vector used to define manifold (degrees),
%           .y_pc = Px1 vector containing y coordinate of phase
%           centers (meters) obtained from lever_arm.m
%           .z_pc = Px1 vector containing z coordinate of phase
%           centers (meters) obtained from lever_arm.m
%           .fs = sampling frequency (Hz)
%           .f0 = chirp start frequency (Hz)
%           .f1 = chirp stop frequency (Hz)
%           .fc = carrier frequency (Hz)
%           .widening factor = maximum number of samples for  a signal to
%           propagate over the aperture,
%           .ft_wind = fast time windowing function (i.e. @boxcar or
%           @hanning),
%           .noise_en = boolean used to turn off noise (for debugging),
%
% See also:
%   lever_arm.m
% =========================================================================
c = 2.997924580003452e+08; % physical_constants too slow

if ~exist('param.src.noise.en','var')
  param.src.noise.en = 1;
end

% -------------------------------------------------------------------------
% Setup processing fields
% -------------------------------------------------------------------------
Psig_dB     = param.src.SNR;
Ps          = 10.^(Psig_dB./10);
sigma_n     = sqrt(1/2);                    % assumes unity noise power
Q           = numel(param.src.DOAs);     % number of sources
if isfield(param.src,'DOAs_MP')
  Q_MP = numel(param.src.DOAs_MP);
else
  param.src.DOAs_MP = [];
  Q_MP = 0;
end
Num_sens    = length(param.src.y_pc); % number of sensors
fs          = param.src.fs;               % sampling frequency
BW          = param.src.f1 - param.src.f0; % chirp BW
dt          = 1/fs;                         % sampling interval
fc          = param.src.fc;               % carrier frequncy 
W           = param.method.wb_td.widening_factor;  % widening factor
NB          = param.method.wb_fd.filter_banks;  % Number of narrowband filter banks
Nc          = length(param.src.y_pc);

% Transmit beamforming
src_params = param.src;
if isfield(src_params,'tx_weights') && ~isempty(src_params.tx_weights)
  comp_tx_weight = src_params.tx_weights;
else
  comp_tx_weight = ones(Q+Q_MP,1);
end
  
% Array calibration errors
if isfield(param,'error_params') && ~isempty(param.error_params)
  error_params   = param.error_params;
  error_ypc      = error_params.error_ypc;
  error_zpc      = error_params.error_zpc;
  error_phase    = error_params.error_phase;
  error_g_s      = error_params.error_g_s;
  error_g_p      = error_params.error_g_p;
  error_g_offset = error_params.error_g_offset;
else
  error_ypc      = zeros(Nc,1);
  error_zpc      = zeros(Nc,1);
  error_phase    = zeros(Nc,1);
  error_g_s      = zeros(Nc,1);
  error_g_p      = zeros(Nc,1);
  error_g_offset = zeros(Nc,1);
end

% Mutual coupling matrix
if isfield(param.src,'mutual_coup_mtx') && ~isempty(param.src.mutual_coup_mtx)
  C = param.src.mutual_coup_mtx;
else
  C = eye(Nc);
end

% Multipath delays and weights
if isfield(param.src,'tau_MP') && ~isempty(param.src.tau_MP)
  tau_MP = param.src.tau_MP;
else
  tau_MP = zeros(Q+Q_MP,1);
end

if isfield(param.src,'w_MP') && ~isempty(param.src.w_MP)
  w_MP = param.src.w_MP;
else
  w_MP = zeros(Q+Q_MP,1);
end

% Setup matrix of time delays of each channel for the given DOAs
% -------------------------------------------------------------------------
if ~isempty(param.src.DOAs)
  for idx = 1:Q
    doa(:,idx) = [param.src.DOAs(idx),param.src.DOAs_MP(1+(idx-1)*Q_MP:idx*Q_MP)];
  end
  doa = doa(:).';
  
  uy          = sin(doa.*(pi/180));
  uz          = sqrt(1-uy.^2);
  Tau_mtx     = (2/c).*(-1.*(param.src.z_pc+error_zpc)*uz + (param.src.y_pc+error_ypc)*uy);
end
% -------------------------------------------------------------------------
% Determine number of fast time samples needed
% -------------------------------------------------------------------------
dy_max          = abs(param.src.y_pc(1) - param.src.y_pc(end));
tau_max         = 2*dy_max/c;
wrap_bin        = ceil(abs(tau_max*fs));
Nt              = param.src.Nsnap + (W-1);
Nt              = Nt + (2*wrap_bin+1);

% -------------------------------------------------------------------------
% Setup baseband and passband frequency vectors for time delaying channels.
% -------------------------------------------------------------------------
df              = (1/(Nt*dt));
f_bb            = (-floor(Nt/2)*df : df : floor((Nt -1) /2)*df).';
f_pb            = ifftshift(f_bb + fc);
time            = (0:dt:(Nt-1)*dt).';

% -------------------------------------------------------------------------
% Setup a vector of complex SAR pixel values measured with respect to
% reference element
% -------------------------------------------------------------------------
S               = complex(zeros(Q,Nt));
% Ps              = Ps./(Num_sens);

tmp_idx = 1;
for idx = 1:Q;
  weight      = comp_tx_weight(tmp_idx);
  sigma_s     = sqrt(Ps(tmp_idx)/2);
  S(tmp_idx,:) = weight*sigma_s.*(randn(1,Nt) - 1i*randn(1,Nt));
  
  % MPCs associated with this target
  for idx_MP = 1:Q_MP
    sigma_MP = sqrt((10.^(w_MP(idx_MP+tmp_idx)./10))/2);
    
    S(idx_MP+tmp_idx,:) =  sigma_MP * S(tmp_idx,:);
  end
  tmp_idx = tmp_idx + Q_MP + 1;
end


% -------------------------------------------------------------------------
% Setup impulse response model and corresponding time series
% -------------------------------------------------------------------------
Ts          = 1/fs;
Ttot        = 10*(tau_max + W*Ts);
Lwin        = ceil(Ttot/Ts);
% Force number of samples in frequency domain window to be odd
if ~mod(Lwin,2)
  Lwin = Lwin + 1;
end
df          = 1/(Lwin*Ts);
wind_func   = param.src.ft_wind;

% Compute impulse response
Hwin        = wind_func(Lwin);
Hwin        = ifftshift(Hwin);
Hwin        = [Hwin(1:floor(Lwin/2)+1);zeros(10*Lwin,1);Hwin(floor(Lwin/2)+2:end)];
Nfft        = length(Hwin);
dt          = 1/(Nfft*df);
% Store impulse response time series and corresponding time
% vector in array param
imp_resp.time_vec   = dt.*[0:floor((Nfft-1)/2), -floor(Nfft/2):-1];
imp_resp.time_vec   = fftshift(imp_resp.time_vec);
imp_resp.vals       = fftshift(ifft(Hwin));
imp_resp.vals       = imp_resp.vals./max(imp_resp.vals);

% -------------------------------------------------------------------------
% Define freq domain window
% -------------------------------------------------------------------------
good_idxs       = find(f_bb > -1 * BW/2 & f_bb < BW/2);
Hwin            = zeros(Nt,1);
Hwin(good_idxs) = param.src.ft_wind(length(good_idxs));
Hwin            = ifftshift(Hwin);

% -------------------------------------------------------------------------
% Window sources to obtain modeled complex envelope of signal from each 
% direction
% -------------------------------------------------------------------------
for q_idx = 1:Q+Q_MP
   Sref_fd = fft(S(q_idx,:));
   Sref_fd = Sref_fd(:);
   Sref_fd = Sref_fd.*Hwin;   
   
   S(q_idx,:) = ifft(Sref_fd);      
end

% -------------------------------------------------------------------------
% Build MxNtxQ data cube of signal measurements across array over Nt
% samples
% -------------------------------------------------------------------------
tmp_dat     = complex(zeros(Num_sens,Nt,Q+Q_MP));

% Create delayed version of each windowed signal 
if ~isempty(param.src.DOAs)
    for q_idx   = 1:Q+Q_MP
        tau_vec  = Tau_mtx(:,q_idx);
        Sref_fd  = fft(S(q_idx,:));
        Sref_fd  = Sref_fd(:);
        
        tmp_doa = doa(q_idx).*(pi/180);
        gain_error_exp = error_g_s.*(sin(tmp_doa)-sin(error_g_p)).^2 + error_g_offset;
        gain_error  = exp(-gain_error_exp./2);
        %   gain_error  = 10.^(gain_error_exp./20);
        %   gain_error  = 1+gain_error_exp;
        %   gain_error  = 10.^(-gain_error_exp./20);
        phase_error = error_phase;
        pg_error = gain_error .* exp(1i*phase_error);
        phase_delay_MP = exp(-1i*2*pi*f_pb*tau_MP(q_idx));
        for m_idx = 1:Num_sens
            tau = tau_vec(m_idx);
            tmp_dat(m_idx,:,q_idx) = ifft(phase_delay_MP.*Sref_fd.*exp(1i*2*pi*f_pb.*tau)*pg_error(m_idx));
            %         tmp_dat(m_idx,:,q_idx) = ifft(Sref_fd.*exp(1i*2*pi*f_pb.*tau));
            %         tmp_dat(m_idx,:,q_idx) = ifft(Sref_fd.*exp(1i*2*pi*f_pb.*tau)*exp(-1i*2*pi*fc.*tau));
        end
    end
end
% =========================================================================
% Create MxNt matrix of additive white Gaussian noise
% =========================================================================
AWGN            = complex(zeros(Num_sens,Nt));
for idx = 1:Num_sens
    noise       = sigma_n.*(randn(1,Nt) - 1i*randn(1,Nt));     
    noise_fd    = fft(noise);
    noise_fd    = noise_fd(:);
    
    AWGN(idx,:) = ifft(noise_fd.*Hwin);
end

% =========================================================================
% Sum plane waves incident on each element by summing over Q dimension
% =========================================================================
array_data      = sum(tmp_dat,3);

% Account for mutual coupling effect
array_data = C*array_data;

% Add noise
% -------------------------------------------------------------------------
if param.src.noise.en
    array_data  = array_data + AWGN;
end

% Truncate Data
% Data = M * Nt
%   M = number of sensors
%   Nt = number of snapshots
% -------------------------------------------------------------------------
Data        = array_data(:,wrap_bin+1:end-wrap_bin-1);

% =========================================================================
% Estimate space time covariance matrix (time domain)
% =========================================================================
[Num_sens,Nt] = size(Data);
Nsnap         = Nt - (W-1);

DCM = complex(zeros(Num_sens*W,Num_sens*W));
for idx = 1:Nsnap
  x   = reshape(Data(:,idx:idx+W -1),Num_sens*W,1);
  x   = x(:);
  DCM = DCM + x*x';
end
DCM = (1/Nsnap)*DCM;

% =========================================================================
% Estimate space time covariance matrix (frequency domain)
% =========================================================================

% Perform DFT for each set of data samples
DCM_fd = complex(zeros(Num_sens*NB,Num_sens));
for idx = 1:Nsnap/NB
  x_nb = fft(Data(:,(idx-1)*NB+(1:NB)),[],2);
  for nb = 1:NB
    DCM_fd((nb-1)*Num_sens+(1:Num_sens),:) = DCM_fd((nb-1)*Num_sens+(1:Num_sens),:) + x_nb(:,nb)*x_nb(:,nb)';
  end
end
% DCM_fd = (1/Nsnap)*DCM_fd; % John: divide by total number of snapshots
DCM_fd = (1/(Nsnap/NB))*DCM_fd; % Mohanad: divide by the number of f-domain snapshots

return

