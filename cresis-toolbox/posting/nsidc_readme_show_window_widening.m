% nsidc_show_window_widening
%
% Demonstrates pulse compression with time domain and frequency domain
% windowing for the read me.
%
% Author: John Paden
%
% See also: type "nsidc_help.m"

rx_window_guards = 1:0.05:1;
for guard_idx = 1:length(rx_window_guards)
  rx_window_guard = rx_window_guards(guard_idx);
  
  fs = 120e6;
  dt = 1/fs;
  Tpd = 1e-6;
  BW = 20e6;
  f0 = -BW/2;
  alpha = BW / Tpd;
  Nt = round(Tpd/dt * 3);
  time = (0:Nt-1).' * dt;
  
  Hwin = tukeywin(Tpd*fs,0.);
  Hwin = Hwin(2:end-1);
  
  signal = exp(j*2*pi*f0*time + j*pi*alpha*time.^2);
  signal(1:length(Hwin)) = signal(1:length(Hwin)) .* Hwin;
  signal(length(Hwin)+1:end) = 0;
  
  ref = conj(fft(signal));
  T = dt*Nt;
  df = 1/T;
  freq = (-floor(Nt/2) : floor((Nt-1)/2)) * df;
  Nt_win = length(find(abs(freq) < rx_window_guards*BW/2));
  Hwin_freq = hanning(Nt_win);
  start_idx = find(freq >= -rx_window_guards*BW/2,1);
  ref = fftshift(ref);
  ref(start_idx - 1 + (1:Nt_win)) = ref(start_idx - 1 + (1:Nt_win)) .* Hwin_freq;
  ref(1:start_idx-1) = 0;
  ref(start_idx + Nt_win:end) = 0;
  ref = ifftshift(ref);
  
  Mt = 50;
  output = ifft( fft(signal) .* ref);
  output = interpft(output,Mt*length(output));
  time_interp = (0:Mt*Nt-1).' * dt/Mt;
  output = output / max(output);
  
  figure(1); clf;
  plot(freq/1e6,fftshift(lp(ref)));
  hold on;
  plot(freq/1e6,fftshift(lp(fft(signal))),'r:');
  hold off;
  
  figure(2); clf;
  plot(time_interp*1e6, lp(output));
  xlabel('Time (us)');
  ylim([-150 0]);
  % xlim([1 12000]);
  title(sprintf('rx window guard %.2f',rx_window_guards));
  grid on;
  
  time_resolution_ideal = 1/BW*1e9;
  fprintf('Ideal resolution %f ns\n', time_resolution_ideal);
  time_resolution_16dB = time_interp(find(lp(output)<-16,1))*1e9;
  fprintf('16 dB resolution %f ns\n', time_resolution_16dB);
  fprintf('Widening %.0f%%\n', 100*time_resolution_16dB / time_resolution_ideal);
  fprintf('Expected sidelobe level %.1f dB\n', 20*log10(Tpd*BW))
  
  if guard_idx < length(rx_window_guards)
    pause(1)
  end
  
end

for Ny = [5 100 1000]
beam_pattern = ifft(hanning(Ny),Ny*100);
beam_pattern = beam_pattern/max(beam_pattern);
plot(lp(beam_pattern))
Ny
find(lp(beam_pattern)<-16,1)
end


return;
