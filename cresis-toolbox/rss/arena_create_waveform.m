% script arena_create_waveform
%
% Creates a waveform file for an arena arbitrary waveform generator/DAC.
% This is just a file with N double values in it. Waveform full scale
% is +/- 1.0 (double).

out_fn = 'C:\Temp\test_waveform.bin';
fs = 2480e6;
Tpd = 200e-6;
Nt = round(Tpd*fs);
f0 = (2e9+20e9)/16
f1 = (4e9+20e9)/16
BW = (f1 - f0)
K = -BW/Tpd

dt = 1/fs;
time = dt*(0:Nt-1).';

wf = 0.5*tukeywin(Nt,0) .* cos(2*pi*(fs-f0)*time + pi*K*time.^2);

fid = fopen(out_fn,'w');
fwrite(fid,wf,'double');
fclose(fid);

