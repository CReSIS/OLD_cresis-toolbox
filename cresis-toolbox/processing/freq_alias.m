function freq = freq_alias(freq,fs)
% freq = freq_alias(freq,fs)
%
% Function for handling aliasing due to a DDC shift of digital data sampled
% at fs.

freq = mod(freq + fs/2,fs)-fs/2;
