
format long
fn = 'C:\Users\dangermo\OneDrive - University of Kansas\Desktop\UTIG\X53b\RADnh5\bxds.10000';
% fn = 'C:\Users\dangermo\OneDrive - University of Kansas\Desktop\UTIG\X53a_first_flight_clean\RADnh5\bxds.10000';

% ===================================================================
%% Check input arguments
% ===================================================================
if ~exist('param','var')
  param = struct();
end
if ~isfield(param,'recs') || isempty(param.recs)
  param.recs = [0 inf];
end
if ~isfield(param,'first_byte') || isempty(param.first_byte)
  param.first_byte = 0;
end

% Reset/clear hdr struct
hdr = [];

% ===============================================================
%% Open file big-endian for reading
% ===============================================================
% [fid,msg] = fopen(fn,'r','ieee-le');
[fid,msg] = fopen(fn,'r','ieee-be');
if fid < 1
  fprintf('Could not open file %s\n', fn);
  error(msg);
end

%% Slow Reader

%% Read first header
% Bytes: 0-1
nsamp = fread(fid,1,'uint16');
% Bytes: 2
nchan = fread(fid,1,'uint8');
% Bytes: 3
vr0 = fread(fid,1,'uint8');
% Bytes: 4
vr1 = fread(fid,1,'uint8');
% Bytes: 5
choff = fread(fid,1,'uint8');
% Bytes: 6
ver = fread(fid,1,'uint8');
% Bytes: 7
resvd = fread(fid,1,'uint8');
% Bytes: 8
absix = fread(fid,1,'double');
% Bytes: 16
relix = fread(fid,1,'double');
% Bytes: 24
xinc = fread(fid,1,'single');
% Bytes: 28
rseq = fread(fid,1,'uint32');
% Bytes: 32
scount = fread(fid,1,'uint16');
% Bytes: 34
tscount = fread(fid,1,'uint16');
% Bytes: 36
rtime = fread(fid,tscount,'double');
% Bytes: 36 + 8*tscount

header_rec_size = 36 + 8*tscount + 17*2;
data_rec_size = 2*nsamp*nchan;
rec_size = header_rec_size + data_rec_size;
if nchan ~= 2
  error('This file has nchan=%d. Only nchan=2 supported.', nchan);
end
fseek(fid,0,1);
file_size = ftell(fid);
num_rec_in_file = floor(file_size/rec_size/2)-1;
if param.recs(1) + param.recs(2) > num_rec_in_file
  param.recs(2) = num_rec_in_file - param.recs(1);
end

%% Read data
fseek(fid,rec_size*param.recs(1),-1);

Nc = 4;
data = cell(Nc,1);
hdr = cell(Nc,1);
rec = zeros(Nc,1);
for chan = 1:Nc
  data{chan} = zeros(nsamp,param.recs(2));
  hdr{chan}.nsamp = zeros(1,param.recs(2));
  hdr{chan}.choff = zeros(1,param.recs(2));
  hdr{chan}.tscount = zeros(1,param.recs(2));
end

while any(rec < param.recs(2))
  if fread(fid,1,'uint16') ~= 3200
    error('Bad record');
  end
  fseek(fid,3,0);
  choff = fread(fid,1,'uint8');
  choff = bitand(choff,bin2dec('00011111'));
  if rec(choff+1) >= param.recs(2)
    fseek(fid,-6 + header_rec_size+data_rec_size,0);
    continue
  end
  rec(choff+1) = rec(choff+1) + 1;
  rec(choff+2) = rec(choff+2) + 1;
  fseek(fid,-6,0);

  % Bytes: 0-1
  hdr{choff+1}.nsamp(rec(choff+1)) = fread(fid,1,'uint16');
  % Bytes: 2
  hdr{choff+1}.nchan(rec(choff+1)) = fread(fid,1,'uint8');
  % Bytes: 3
  hdr{choff+1}.vr0(rec(choff+1)) = fread(fid,1,'uint8');
  % Bytes: 4
  hdr{choff+1}.vr1(rec(choff+1)) = fread(fid,1,'uint8');
  % Bytes: 5
  hdr{choff+1}.choff(rec(choff+1)) = fread(fid,1,'uint8');
  % Bytes: 6
  hdr{choff+1}.ver(rec(choff+1)) = fread(fid,1,'uint8');
  % Bytes: 7
  hdr{choff+1}.resvd(rec(choff+1)) = fread(fid,1,'uint8');
  % Bytes: 8
  hdr{choff+1}.absix(rec(choff+1)) = fread(fid,1,'double');
  % Bytes: 16
  hdr{choff+1}.relix(rec(choff+1)) = fread(fid,1,'double');
  % Bytes: 24
  hdr{choff+1}.xinc(rec(choff+1)) = fread(fid,1,'single');
  % Bytes: 28
  hdr{choff+1}.rseq(rec(choff+1)) = fread(fid,1,'uint32');
  % Bytes: 32
  hdr{choff+1}.scount(rec(choff+1)) = fread(fid,1,'uint16');
  % Bytes: 34
  hdr{choff+1}.tscount(rec(choff+1)) = fread(fid,1,'uint16');
  % Bytes: 36
  hdr{choff+1}.rtime{rec(choff+1)} = fread(fid,tscount,'double');
  % Odd stuff
  fread(fid,17,'uint16');

  data{choff+1}(:,rec(choff+1)) = fread(fid,hdr{choff+1}.nsamp(rec(choff+1)),'int16');
  data{choff+2}(:,rec(choff+2)) = fread(fid,hdr{choff+1}.nsamp(rec(choff+1)),'int16');

  if 0
    figure(1); clf;
    plot(data{choff+1}(:,rec(choff+1)));
    hold on
    plot(data{choff+1}(:,rec(choff+1)));
    ylim([-1000 1000]);
    keyboard
  end
end

fclose(fid);

% 50 MHz
% 52.5-67.5 MHz
% 2.5-17.5 MHz
fs = 50e6;
Nt = 3200;
dt = 1/fs;
df = fs/Nt;
time = dt*(0:Nt-1).';
freq = df*(-floor(Nt/2):floor((Nt-1)/2));

Nt_512 = 512;
df_512 = fs/Nt_512;
freq_512 = df_512*(-floor(Nt_512/2):floor((Nt_512-1)/2));

dt_x = 32*160e-6;
Nx = param.recs(2);
time_x = dt_x*(0:Nx-1);
BW_x = 1/dt_x;
df_x = BW_x/Nx;
freq_x = df_x*(-floor(Nx/2):floor((Nx-1)/2));

% figure(1); clf;
for chan = 1:4
  if 0
    h_fig = figure(chan); clf;
    set(h_fig,'WindowStyle','docked');
  else
    subplot(2,2,chan);
  end
  if 0
    imagesc(db(data{chan} - mean(data{chan},2)));
  elseif 0
    % range-doppler
    imagesc(db(fft(data{chan},[],2)));
  elseif 0
    % frequency-space
    imagesc([],freq_512/1e6,fftshift(db(fft(data{chan}(2500+(1:Nt_512),:))),1));
  elseif 1
    % frequency-space PSD
    plot(freq_512/1e6, db(mean(abs(fft(data{chan}(2500+(1:Nt_512),:),[],1)).^2,2),'power'));
    grid on;
    xlabel('Frequency (MHz)')
    ylabel('Relative power (dB)')
  elseif 0
    % range-doppler noise/late-record
    imagesc(fftshift(db(fft(data{chan}(2500+(1:Nt_512),:),[],2)),1));
  elseif 1
    % range-doppler noise/late-record PSD
    plot(freq_x, fftshift(db(mean(abs(fft(data{chan}(2500+(1:Nt_512),:),[],2).^2),1),'power'),2));
    grid on;
    xlabel('Spatial frequency (Hz)')
    ylabel('Relative power (dB)')
  elseif 1
    % f-k
    imagesc(freq_x,freq_512/1e6,fftshift(db(fft2(data{chan}(2500+(1:Nt_512),:)))));
    xlabel('Spatial frequency (Hz)')
    ylabel('Frequency (MHz)')
  end
  if chan == 1 || chan == 3
    title(sprintf('Chan %d (low gain)',chan-1))
  else
    title(sprintf('Chan %d (high gain)',chan-1))
  end
  colormap(1-gray(256))
  caxis([0 150])
  hold on;
end
link_figures;
xlim([-25 25]);
