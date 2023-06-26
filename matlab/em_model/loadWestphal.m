% Loads data from Westphal_Jiracek_ice_data.txt
%
% The e' and e'' are high, so I would recommend not correcting
% for density (e.g. tuto density corrections are commented out).

dataWJ = load('Westphal_Jiracek_ice_data.txt');

ward_glacial.freq = [150 300 500 1000 2700]*1e6;
ward_glacial.density = dataWJ(1,1);
ward_glacial.T = dataWJ(1:8,2);
ward_glacial.e_real = dataWJ(1:8,3:7);
ward_glacial.e_imag = dataWJ(1:8,8:12) .* ward_glacial.e_real;
ward_glacial.e_real = ward_glacial.e_real ...
  .* (1 + 1.7*0.917 + 0.7*0.917^2) ...
  ./ (1 + 1.7*ward_glacial.density + 0.7*ward_glacial.density^2);
ward_glacial.e_imag = ward_glacial.e_imag ...
  .* (0.52*0.917 + 0.62*0.917.^2) ...
  ./ (0.52*ward_glacial.density + 0.62*ward_glacial.density.^2);

tuto.freq = [150 300 500 1000 2700]*1e6;
tuto.density = dataWJ(9,1);
tuto.T = dataWJ(9:16,2);
tuto.e_real = dataWJ(9:16,3:7);
tuto.e_imag = dataWJ(9:16,8:12) .* tuto.e_real;
tuto.e_real = tuto.e_real ...
  .* (1 + 1.7*0.917 + 0.7*0.917^2) ...
  ./ (1 + 1.7*tuto.density + 0.7*tuto.density^2);
tuto.e_imag = tuto.e_imag ...
  .* (0.52*0.917 + 0.62*0.917.^2) ...
  ./ (0.52*tuto.density + 0.62*tuto.density.^2);

little.freq = [150 300 500 1000 2700]*1e6;
little.density = dataWJ(17,1);
little.T = dataWJ(17:23,2);
little.e_real = dataWJ(17:23,3:7);
little.e_imag = dataWJ(17:23,8:12) .* little.e_real;
little.e_real = little.e_real ...
  .* (1 + 1.7*0.917 + 0.7*0.917^2) ...
  ./ (1 + 1.7*little.density + 0.7*little.density^2);
little.e_imag = little.e_imag ...
  .* (0.52*0.917 + 0.62*0.917.^2) ...
  ./ (0.52*little.density + 0.62*little.density.^2);

arctic.freq = [150 300 500 1000 2700]*1e6;
arctic.density = dataWJ(25,1);
arctic.T = dataWJ(25:32,2);
arctic.e_real = dataWJ(25:32,3:7);
arctic.e_imag = dataWJ(25:32,8:12) .* arctic.e_real;
arctic.e_real = arctic.e_real ...
  .* (1 + 1.7*0.917 + 0.7*0.917^2) ...
  ./ (1 + 1.7*arctic.density + 0.7*arctic.density^2);
arctic.e_imag = arctic.e_imag ...
  .* (0.52*0.917 + 0.62*0.917.^2) ...
  ./ (0.52*arctic.density + 0.62*arctic.density.^2);

ward_sea.freq = [150 300 500 1000 2700]*1e6;
ward_sea.density = dataWJ(33,1);
ward_sea.T = dataWJ(33:end,2);
ward_sea.e_real = dataWJ(33:end,3:7);
ward_sea.e_imag = dataWJ(33:end,8:12) .* ward_sea.e_real;
ward_sea.e_real = ward_sea.e_real ...
  .* (1 + 1.7*0.917 + 0.7*0.917^2) ...
  ./ (1 + 1.7*ward_sea.density + 0.7*ward_sea.density^2);
ward_sea.e_imag = ward_sea.e_imag ...
  .* (0.52*0.917 + 0.62*0.917.^2) ...
  ./ (0.52*ward_sea.density + 0.62*ward_sea.density.^2);

clear('dataWJ');

return

freq = tuto.freq;
figure(1); clf;
% Plot -1, -10, and -30
plot(freq/1e6,tuto.e_real(1,:),'rx-');
hold on;
plot(freq/1e6,tuto.e_real(3,:),'gx-');
plot(freq/1e6,tuto.e_real(5,:),'bx-');

w = 2*pi*freq;
e0 = 8.854e-12;
u0 = 4e-7*pi;
sigma = 0;
for index=1:length(freq)
  er(:,index) = ice(273.15+[-1 -10 -30],0.917,freq(index),3e-6,0.22);
  atten(:,index) = 20*log10(exp(1)) ...
    *real(j*sqrt(-j*w(index)*u0*(sigma+j*w(index)*er(:,index)*e0)));
  tuto.atten(:,index) = 20*log10(exp(1)) ...
    *real(j*sqrt(-j*w(index)*u0*(sigma+j*w(index) ...
      *(tuto.e_real(:,index)-j*tuto.e_imag(:,index))*e0)));
end

plot(freq/1e6,real(er(1,:)),'rx--');
plot(freq/1e6,real(er(2,:)),'gx--');
plot(freq/1e6,real(er(3,:)),'bx--');
hold off;
title('e'', Westphal (solid lines), Model (dashed lines)');
xlabel('Frequency (MHz)');
ylabel('er (real)');

figure(2); clf;

plot(freq/1e6,tuto.atten(1,:),'rx-');
hold on;
plot(freq/1e6,tuto.atten(3,:),'gx-');
plot(freq/1e6,tuto.atten(5,:),'bx-');

plot(freq/1e6,atten(1,:),'rx--');
plot(freq/1e6,atten(2,:),'gx--');
plot(freq/1e6,atten(3,:),'bx--');
hold off;
title('{\alpha}, Westphal (solid lines), Model (dashed lines)');
xlabel('Frequency (MHz)');
ylabel('dB/m');

