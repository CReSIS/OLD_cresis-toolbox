function waveform = gen_randn_wave(N,corr_len,scale)
% waveform = gen_randn_wave(N,corr_len,scale)
% 
% Generates a real waveform from a Gaussian random process that is filtered
% with a 3 dB bandwidth of 1/corr_length and is scaled by scale. Some
% effort is taken to remove transients caused by the filtering process.
% 
% N = number of samples in output waveform
% corr_len = approx. correlation length of filter (used to compute the
%   filter's 3 dB bandwidth as 1/corr_len)
% scale = output is scaled by this amount
%
% Author: John Paden

[B_distort,A_distort] = butter(2,1/corr_len);
Nt_distort = N + round(corr_len*2);
waveform = filtfilt(B_distort,A_distort,scale*randn(Nt_distort,1));
waveform = waveform(1+round(corr_len*2):end);

return;
