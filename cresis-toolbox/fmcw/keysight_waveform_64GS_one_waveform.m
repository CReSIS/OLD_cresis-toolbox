close all; clear all;

fs = 64e9; %512 samples corresponds to 1 cycle at 125MHz
dt = 1/fs;
tau = 240e-6;
max_del = 10e-6;
nt_tau = round(tau/dt);
nt = round((tau+max_del)/dt);
tt = (0:dt:(nt-1)*dt)';

f0 = 2e9;
bw = 16e9;
k = bw/tau;

chirp1 = [tukeywin(nt_tau,.1);zeros(nt-nt_tau,1)].*sin(2*pi*(f0*tt + 0.5*k*(tt.*tt)));

m_clk = zeros(nt,1);

%wave_id is used to uniquely identify the current waveform
%it is written to the header of the data file every epri through marker 2
%64e9/4096 is 125/8 (using 3rd bit (number 2) of the counter
wave_id = uint64(hex2dec('464D4357')); 
for index = 1:32,
    if mod(bitshift(wave_id,1-index),2),
      m_clk((index-1)*4096+2*2048:(index-1)*4096+2*2048+4095) = ones(4096,1);
    end;
end;
wave_id = uint64(hex2dec('4F49425F'));
for index = 33:64,
    if mod(bitshift(wave_id,33-index),2),
      m_clk((index-1)*4096+2*2048:(index-1)*4096+2*2048+4095) = ones(4096,1);
    end;
end;

data = zeros(nt,4);
data(:,1) = chirp1;

%marker 1 is used as the pri trigger which is sent from the AWG to the DAQ
%the location of the pri trigger can be adjusted to account for setup/hold
%instabilities due to cable lengths
%64e9/4096 is 125/8
data(0+(4097:8192),3) = ones(4096,1); %PRI Trigger
data(:,4) = m_clk(1:nt);

%write the file
fid1 = fopen('fmcw_1ch_64GSPS_2to18GHz_240us_4khz.csv','w');
fprintf(fid1,'SampleRate=64000000000\r\n');
fprintf(fid1,'SetConfig=true\r\n');
fprintf(fid1,'Y1,Y2,SampleMarker1,SampleMarker2\r\n');
fprintf(fid1,'%12.10f,%12.10f,%i,%i\r\n',data');
fclose(fid1);