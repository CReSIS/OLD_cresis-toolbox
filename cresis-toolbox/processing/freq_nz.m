function freq = freq_nz(freq,fs,nz)
% freq = freq_nz(freq,fs,nz)
%
% Function for handling undersampling by knowing which nyquist zone the is
% sampled in.
%
% Example:
% fs = 10;
% f_ddc = 13;
% nz = 2;
% freq = -5:0.1:4.9;
% plot(freq,freq_nz(freq_alias(freq+f_ddc,fs),fs,nz),'.'); grid on;

freq = freq + (2*(freq>=0)-1)*(-1)^nz*round(nz/2)*fs;
