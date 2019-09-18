clear
close all
tic
taup=1; % uncompressed pulse width
BW=100; % chirp bandwidth
f0=-50;   
% time bandwidth product
time_B_product = BW * taup;
if(time_B_product < 5 )
    fprintf('************ Time Bandwidth product is TOO SMALL ***********')
    fprintf('\n Change b and or taup')
    return
end
%% (1) LFM single pulse
% =========================================================================

Nt = 500; fs = 500; dt = 1/fs; tp = dt*(0:Nt-1);
alpha=BW/taup;
LFM = exp(1i*2*pi*f0*tp + 1i * pi * alpha .* tp.^2);
figure
subplot(3,1,1)
plot(tp,real(LFM)); ylabel('Amplitude'); xlabel('Time (sec)'); 
title('Real part of LFM'); grid
subplot(3,1,2)
plot(tp,imag(LFM)); ylabel('Amplitude'); xlabel('Time (sec)'); 
title('Imaginary part of LFM'); grid
LFMFFT = fftshift(fft(LFM)); df=fs/Nt;
if mod(Nt,2)==0    % even number
    f=(-Nt/2:1:Nt/2-1)*df;
else               % Odd number
    f=(-(Nt-1)/2:1:(Nt-1)/2)*df;
end
subplot(3,1,3) % Spectrum
plot(f, 20*log10(abs(LFMFFT)./max(abs(LFMFFT)))); xlabel('Frequency (Hz)');
ylabel('Normalized power,dB'); title('Spectrum for LFM waveform'); 
axis([-100 100 -40 1]); grid
%% (2) Autocorrelation LFM pulse
% =========================================================================

Rxx= zeros(1,Nt); % initialize
for L = 1:Nt
    Rxx(L) = (1/Nt).*LFM(L:Nt) * LFM(1:Nt-L+1)';
    % The autocorrelation of a real, stationary signal x(t) is defined by
    % Rxx(tau) = E[x(t) x(t+tau)]
end
Rxx_n=conj(Rxx(end:-1:2)); % By symmetry property, Rxx(-tau) = Conj(Rxx(tau))
Rxx=[Rxx_n Rxx]; % combine the Rxx for -tau and +tau
Nr=size(Rxx,2);
lag = (-Nr/2:1:Nr/2-1);
% Plotting from own MATLAB autocorrelation code
figure
subplot(2,2,1)
plot(lag,abs(Rxx)); ylabel('R_{xx}'); xlabel('Delay index'); 
title({'Autocorrelation of LFM pulse','using own code'});grid on;grid minor
% Plotting autocorrelation using MATLAB inbuilt XCORR function
subplot(2,2,2)
plot(lag,abs(xcorr(LFM)/size(LFM,2))); ylabel('R_{xx}'); 
xlabel('Delay index');title({'Autocorrelation of LFM pulse','using XCORR'}); 
grid  on; grid minor
% Plotting Autocorrelation on dB scale
subplot(2,2,3)
plot(lag,20*log10(abs(Rxx))); ylabel('R_{xx} (dB)'); xlabel('Delay index'); 
title('Autocorrelation of LFM pulse');  grid on; grid minor; 
axis([-500 500 -60 1]); % axis([-1 1 -60 1]) %
dF=fs/Nr;
if mod(Nr,2)==0    % even number
    fr=(-Nr/2:1:Nr/2-1)*dF;
else               % Odd number
    fr=(-(Nr-1)/2:1:(Nr-1)/2)*dF;
end
% Plotting PSD of LFM from Autocorrelation
% The Fourier transform of Rxx(tau) is the Power Spectral Density (PSD) Sx(f)
subplot(2,2,4)
PSDRxx=abs(fftshift(fft(Rxx,[],2))); PSDRxx=PSDRxx./max(PSDRxx);
plot(fr,10*log10(PSDRxx)); title('Power spectrum from R_{xx} of LFM'); 
ylabel('Normalized power,dB'); xlabel('Frequency (Hz)');
axis([-100 100 -40 1]); grid on; grid minor

%% (3) Ambiguity Function LFM pulse
%%
% =========================================================================
% the ambiguity fn is defined as: a(t,f) = abs(sumi(u(k)*u'(i-t)*exp(j*2*pi*f*i)))
eps = 0.000001;
ii = 0;
BW=BW;
alpha =  BW/2 /taup; % upchirp, if down chirp pre multiply by -1,  mu =  BW / 2. / taup;
Del=.01;
for tau = -1.1*taup:Del:1.1*taup
    ii = ii + 1;
    jj = 0;
    for fd = -BW:Del:BW
        jj = jj + 1;
        val1 = 1. - abs(tau) / taup;
        val2 = pi * taup * (1.0 - abs(tau) / taup);
        val3 = (fd + alpha * tau);
        val = val2 * val3;
        x(jj,ii) = abs( val1 * (sin(val+eps)/(val+eps))).^2;
    end
end
taux = -1.1*taup:Del:1.1*taup;
fdy = -BW:Del:BW;
figure
mesh(taux,fdy,x)
xlabel ('Delay (sec)')
ylabel ('Doppler (Hz)')
zlabel ('Ambiguity function')
figure
contour(taux,fdy,x)
xlabel ('Delay (sec)')
ylabel ('Doppler (Hz)')
grid on
grid minor
% plotting for fd=0 cut
Del2=0.00001;
taux = -1.5*taup:Del2:1.5*taup;
fd = 0;
alpha = BW /2/  taup;
ii = 0;
for tau = -1.5*taup:Del2:1.5*taup
    ii = ii + 1;
    val1 = 1 - abs(tau) / taup;
    val2 = pi * taup * (1 - abs(tau) / taup);
    val3 = (fd + alpha * tau);
    val = val2 * val3;
    xx(ii) = abs( val1 * (sin(val+eps)/(val+eps)));
end
figure
plot(taux,10*log10(xx.^2)) ; % plot(taux,xx.^2)
grid on
grid minor
xlabel ('Delay (sec)')
ylabel ('Ambiguity (dB)')
axis([-1.500 1.500 -60 1])
%% (4a) Scattering / target Model
% =========================================================================

%----------------------------------------------------------
rng(4);
tau=[0.3  0.67];
N_tgt=size(tau,2); % # of targets
k=1; n=N_tgt;A=3;
GG = A*complex(rand(k,n), rand(k,n)); %random complex scattering coeff
Tgt = zeros(N_tgt,size(tp,2));
da=.00005;
for n=1:N_tgt
    DeltaFF=sinc((tp - tau(n))./da); % lim da-> 0 , (1/a)sinc(x/a) = delta(x)
    tgt=GG(n).*DeltaFF;
    Tgt(n,:)=tgt;
end
targetZ=sum(Tgt);
figure
subplot(3,1,1)
plot(tp,real(targetZ))
xlabel('Delay (sec)'); ylabel('RCS - amplitude'); 
title('Point targets (real part)'); grid on; grid minor
subplot(3,1,2)
plot(tp,imag(targetZ))
xlabel('Delay (sec)'); ylabel('RCS - amplitude'); 
title('Point targets (img part)'); grid on; grid minor
subplot(3,1,3)
plot(tp,abs(targetZ))
xlabel('Delay (sec)'); ylabel('RCS - amplitude'); 
title('Point targets (abs)- Non zero Doppler'); grid on; grid minor
%% (4b) Tx on one pulse
% =========================================================================

S1=LFM;
Ax=[S1 zeros(1,50)]; % Zero padding
Bx=[targetZ zeros(1,50)];
s_received1 = conv(Ax,Bx);
nrx1=size(s_received1,2);
tc = dt*(0:nrx1-1);
figure
plot(tc,abs(s_received1))
title('Transmit waveform (one pulse) convolved with channel response')
xlabel('Time (sec)')
ylabel('Amplitude')
axis([.5 1.5 1.2 7])
grid
%% (5) MF output
% =========================================================================

s_transmitted1 = S1;
s_mixed1 = conj(fft(s_transmitted1,nrx1,2)) .* fft(s_received1,[],2);
xdft1 = ifft(s_mixed1);
psdx1 = (1/(fs*nrx1)) * abs(xdft1).^2;
psdx1 = psdx1./max(psdx1);
figure
plot(tc,10*log10(abs(psdx1)))
xlabel('Delay (Sec) ');
ylabel('Normalized power (dB)');
title('MF output- single pulse,  zero doppler')
grid on; grid minor
axis([0 1 -50 0])
%% (6a) LFM Multiple pulses
% =========================================================================

S1=LFM.';
PRI=2; PRF=1/PRI; fD_max=.5*PRF;
Npls=25;% number of pulses, odd
for rline=1:Npls
    S_train(:,rline)=S1;
end
%% (6b) TX on multiple pulses (Zero Doppler)
% =========================================================================

for nx=1:Npls
    A2=[S_train(:,nx).' zeros(1,50)]; % Zero padding
    B2=[targetZ zeros(1,50)];
    s_received2(:,nx) = conv(A2,B2.');
end
nrx2=size(s_received2,1);
s_transmitted2 = S1; % S_train(:,1);
s_mixed2 = bsxfun(@times,conj(fft(s_transmitted2,nrx2,1)) , fft(s_received2,[],1));
xdft2 = ifft(s_mixed2);
psdx2 = (1/(fs*nrx2)) * abs(xdft2).^2;
psdx2 = psdx2./ max(max(psdx2));
t2 = dt*(0:nrx2-1);
figure
plot(t2,10*log10(abs(psdx2)))
xlabel('Delay (Sec) ');
ylabel('Normalized power (dB)'); %Normalized if psdx2 = psdx2 ./ max(psdx2); is used
title('Multiple pulses, zero doppler')
legend('Range line 1','Range line 2','Range line 3','...','...','...','Location','northeast')
axis([0 1 -45 1])
grid on; grid minor
XX2=fftshift(fft(psdx2,[],2),2);
f1=linspace(-fD_max,fD_max,Npls);
figure
imagesc(t2,f1,10*log10(abs(XX2.')))
xlabel('Delay (Sec)')
ylabel('Doppler frequency (Hz)')
h = colorbar;
ylabel(h, 'Normalized power (dB)')
title('Echogram for multiple pulses, targets with zero doppler')
caxis([-30 0])

[M1,N1]=size(XX2.');
doppler1=linspace(-fD_max,fD_max,M1);
contour(t2,doppler1,abs(XX2.'));
xlabel('Delay (Sec)');
ylabel('Doppler frequency (Hz)');
title('Contour for multiple pulses, targets with zero doppler')
grid on
grid minor

%% Non Zero Doppler for one target
% =========================================================================

fD = [0.1 0]; % assign Doppler
Npls=Npls+30;
% number of pulses, odd, increased to get better Doppler resolution
EX1 = []; EX2 = [];
for ii=1:Npls
    A3=[S1.' zeros(1,10)]; % Zero padding
    B3i=[Tgt(1,:) zeros(1,10)];
    B3ii=[Tgt(2,:) zeros(1,10)];
    s1_received3(:,ii) = conv(A3,B3i);
    s2_received3(:,ii) = conv(A3,B3ii);
    Ns=size(s1_received3,1);
    t3 = dt*(0:Ns-1);
    tss=((ii-1)*PRI + t3).';
    EX1(:,ii)=exp(1i*2*pi.*fD(1).*tss);
    EX2(:,ii)=exp(1i*2*pi.*fD(2).*tss);
    s1_received3(:,ii)=bsxfun(@times ,s1_received3(:,ii),EX1(:,ii));
    s2_received3(:,ii)=bsxfun(@times ,s2_received3(:,ii),EX2(:,ii));
end
s_received3x=s1_received3 + s2_received3; %
s_transmitted3 = S1;
s_mixed3x = bsxfun(@times,conj(fft(s_transmitted3,Ns,1)) , fft(s_received3x,[],1));
% s_mixed3x = conj(fft(s_transmitted3,Ns,1)) .* fft(s_received3,[],1); % run on R2018
Ns = size(s_mixed3x,1);
xdft3x = ifft(s_mixed3x,[],1);
psdx3x = (1/(fs*Ns)) * abs(xdft3x).^2;
psdx3x = psdx3x./max(max(psdx3x));
Ns=size(psdx3x,1);
t3x = dt*(0:Ns-1);
figure
plot(t3x,10*log10(abs(psdx3x))) %,'DisplayName','range line =%')
xlabel('Delay (Sec)')
ylabel('Normalized power (dB)'); %Normalized if psdx2 = psdx2 ./ max(psdx2); is used
title('Multiple pulses,  non zero doppler')
legend('Range line 1','Range line 2','Range line 3','...','...','...','Location','northeast')
axis([0 1 -45 1])
grid on; grid minor

%-----------------------------------------------------------
X3=fftshift(fft(xdft3x,[],2),2);
N6 = size(X3,2);
t4 = dt*(0:N6-1);
dFf=PRF/N6;
if mod(N6,2)==0    % even number
    fd=(-N6/2:1:N6/2-1)*dFf; %*.001;
else               % Odd number
    fd=(-(N6-1)/2:1:(N6-1)/2)*dFf; %*.001;
end
figure
imagesc(t3x,fd,10*log10(abs(X3).'))
h = colorbar;
ylabel(h, 'Relative power (dB)')
title('Echogram for multiple pulses, one target with non zero Doppler')
ylabel('Doppler frequency (Hz)'); xlabel('Delay (Sec)')
clims = caxis;
caxis(clims(2) + [-20 0])
grid on;
ax = gca; % Get handle to current axes.
ax.XColor = 'k';
ax.YColor = 'k';
ax.GridAlpha = 0.7;  % Make grid lines less transparent.
ax.GridLineStyle = '--';
ax.GridColor = [0, 0.5, 0.5];
ax.Layer = 'top';
grid minor

[Mn,Nn]=size(X3.');
doppler2=linspace(-fD_max,fD_max,Mn);
contour(t3x,doppler2,abs(X3.'));
xlabel('Delay (Sec)');
ylabel('Doppler frequency (Hz)');
title('Contour for multiple pulses, one target with non zero Doppler')
grid on
grid minor
%% Along track position vs range plot 
% =========================================================================
xp= -2000:1:2000; yp=0; hp=1000; % Aircraft position (radar)
xt=0; yt=-2000; ht=0; % target ground offset
h=1000;  % height of straight and level flight path
R=sqrt((xp-xt).^2 + (yp-yt).^2 + (hp-ht).^2);  % Calculation of Range
plot(xp,R)
title('Plot of along track range')
xlabel('Aircraft position, X (meters)')
ylabel('Range, R (meters)')
toc