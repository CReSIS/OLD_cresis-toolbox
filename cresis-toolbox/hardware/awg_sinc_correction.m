function corr = sinc_correction(freq, fs, max_f, correction)

freq(freq>max_f) = max_f;

correction.power = -correction.power;
correction.power = correction.power - max(interp1(correction.freq,correction.power,max_f,'linear','extrap'));

if 0
  % Debug/test
  plot(freq,interp1(correction.freq,correction.power,freq,'linear','extrap'));
  keyboard
end

corr = sinc(max_f/fs) ./ sinc(freq/fs) .* 10.^(interp1(correction.freq,correction.power,freq,'linear','extrap')/20);

return;
