function [Data, DCM, imp_resp, DCM_fd] = doa_wideband_data(param)
% [Data, DCM, imp_resp, DCM_fd] = doa_wideband_data(param)
%
% Funcion to simulate multichannel array wideband data.
%
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
physical_constants

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
Nc    = length(param.src.y_pc); % number of sensors
fs          = param.src.fs;               % sampling frequency
BW          = param.src.f1 - param.src.f0; % chirp BW
dt          = 1/fs;                         % sampling interval
fc          = param.src.fc;               % carrier frequncy 
W           = param.method.wb_td.widening_factor;  % widening factor
NB          = param.method.wb_fd.filter_banks;  % Number of narrowband filter banks

if ~isempty(param.src.DOAs)

% Setup matrix of time delays of each channel for the given DOAs
% -------------------------------------------------------------------------
uy          = sin(param.src.DOAs.*(pi/180));
uz          = sqrt(1-uy.^2);
Tau_mtx     = (2/c).*(-1.*param.src.z_pc*uz + param.src.y_pc*uy);


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
% Ps              = Ps./(Nc);

for idx = 1:Q;
    sigma_s     = sqrt(Ps(idx)/2);
    S(idx,:)    =  sigma_s.*(randn(1,Nt) - 1i*randn(1,Nt));
end

% -------------------------------------------------------------------------
% Setup impulse response model and corresponding time series
% -------------------------------------------------------------------------
Ts          = 1/fs;
Ttot        = 10*(tau_max + W*Ts);
Lwin        = ceil(Ttot/Ts);
% Force number of samples in frequency domain window to be
% odd
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
for q_idx = 1:Q
   Sref_fd = fft(S(q_idx,:));
   Sref_fd = Sref_fd(:);
   Sref_fd = Sref_fd.*Hwin;
   S(q_idx,:) = ifft(Sref_fd);   
   
end

% -------------------------------------------------------------------------
% Build MxNtxQ data cube of signal measurements across array over Nt
% samples
% -------------------------------------------------------------------------
tmp_dat     = complex(zeros(Nc,Nt,Q));

if ~isempty(param.src.DOAs)
% Create delayed version of each windowed signal 
for q_idx   = 1:Q
   tau_vec  = Tau_mtx(:,q_idx);
   Sref_fd  = fft(S(q_idx,:));
   Sref_fd  = Sref_fd(:);
  
    for m_idx = 1:Nc
        tau = tau_vec(m_idx);
        tmp_dat(m_idx,:,q_idx) = ifft(Sref_fd.*exp(1i*2*pi*f_pb.*tau));
%         tmp_dat(m_idx,:,q_idx) = ifft(Sref_fd.*exp(1i*2*pi*f_pb.*tau)*exp(-1i*2*pi*fc.*tau));
    end      
end

end
% =========================================================================
% Create MxNt matrix of additive white Gaussian noise
% =========================================================================
AWGN            = complex(zeros(Nc,Nt));
for idx = 1:Nc
    noise       = sigma_n.*(randn(1,Nt) - 1i*randn(1,Nt));     
    noise_fd    = fft(noise);
    noise_fd    = noise_fd(:);
    
    AWGN(idx,:) = ifft(noise_fd.*Hwin);
end

% =========================================================================
% Sum plane waves incident on each element by summing over Q dimension
% =========================================================================
array_data      = sum(tmp_dat,3);

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
[Nc,Nt] = size(Data);
Nsnap         = Nt - (W-1);
Nsnap = 100  ; % simulation as the paper


DCM = complex(zeros(Nc*W,Nc*W));
for idx = 1:Nsnap
  x   = reshape(Data(:,idx:idx+W -1),Nc*W,1);
  x   = x(:);
  DCM = DCM + x*x';
end
DCM = (1/Nsnap)*DCM;
% =========================================================================
% Estimate space time covariance matrix (frequency domain)
% =========================================================================

% Perform DFT for each set of data samples
DCM_fd = complex(zeros(Nc*NB,Nc));
for idx = 1:Nsnap/NB
  x_nb = fft(Data(:,(idx-1)*NB+(1:NB)),[],2);
  for nb = 1:NB
    DCM_fd((nb-1)*Nc+(1:Nc),:) = DCM_fd((nb-1)*Nc+(1:Nc),:) + x_nb(:,nb)*x_nb(:,nb)';
  end
end
DCM_fd = (1/Nsnap)*DCM_fd;

return

