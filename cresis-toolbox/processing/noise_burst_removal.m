noise_window = [3 1];
background_window = [11 51];
zero_window = [11 1];
threshold = 10^(13/10);

for idx=1:size(data,3)
  
  background = fir_dec(abs(data(:,:,idx)).^2,ones(1,background_window(2))/background_window(2),1);
  background = fir_dec(background.',ones(1,background_window(1))/background_window(1),1).';
  
  noise = fir_dec(abs(data(:,:,idx)).^2,ones(1,noise_window(2))/noise_window(2),1);
  noise = fir_dec(noise.',ones(1,noise_window(1))/noise_window(1),1).';
  
  zero = noise > background*threshold;
  zero = fir_dec(zero,ones(1,zero_window(2))/zero_window(2),1);
  zero = fir_dec(zero.',ones(1,zero_window(1))/zero_window(1),1).';
  zero(zero>0) = 1;
  
  fprintf('%.0f %.f\n', idx, sum(zero(:)));
  data(:,:,idx) = data(:,:,idx) .* ~zero;
end

