function alpha = er_to_extinction(freq,er)

% For e0 and u0
physical_constants;

w = 2*pi*freq;

% gamma is the propagation constant
gamma = j*sqrt(-j*w.*u0.*(j*w.*er(1,:)*e0));

% alpha is the attenuation constant
alpha = real(gamma);


return;

er = iceCond(273.15+[-10 0],0.917,195e6,2e-5,0.22,273.15-15);
freq = 195e6;

alpha = er_to_extinction(freq,er);

two_way_loss_per_km_dB = lp(exp(-4*alpha*1000))

round(15./lp(exp(-4*alpha*1000))*1000)
