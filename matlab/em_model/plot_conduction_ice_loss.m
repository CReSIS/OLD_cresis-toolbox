% Plot ice loss from Matsuoka
%
% This file is for illustration purposes.  Not very good...

physicalConstants;

T = [225:273];
er = ice(T,0.917,500e6,0e-6,0.22);
w = 2*pi*500e6;
atten = 20*log10(exp(1))*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
hand1 = plot(T,atten*1000,'r-');

xlabel('Temperature (Kelvin)');
ylabel('Attenuation (dB/km)');
hold on;

T = [225:273];
er = ice(T,0.917,1000e6,0e-6,0.22);
w = 2*pi*1000e6;
atten = 20*log10(exp(1))*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
hand2 = plot(T,atten*1000,'y-');

T = [225:273];
er = ice(T,0.917,500e6,1e-6,0.22);
w = 2*pi*500e6;
atten = 20*log10(exp(1))*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
hand3 = plot(T,atten*1000,'g-');

T = [225:273];
er = ice(T,0.917,1000e6,1e-6,0.22);
w = 2*pi*1000e6;
atten = 20*log10(exp(1))*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
hand4 = plot(T,atten*1000,'b-');

T = [225:273];
er = ice(T,0.917,500e6,2e-6,0.22);
w = 2*pi*500e6;
atten = 20*log10(exp(1))*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
hand5 = plot(T,atten*1000,'m-');

T = [225:273];
er = ice(T,0.917,1000e6,2e-6,0.22);
w = 2*pi*1000e6;
atten = 20*log10(exp(1))*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
hand6 = plot(T,atten*1000,'k-');

axis([223.15 273.15 0 60]);
title('Ice with Impurities (MHz-Micromolarity)');
legend([hand1 hand2 hand3 hand4 hand5 hand6],'Matsuoka-500-0e-6',...
    'Matsuoka-1000-0e-6','Matsuoka-500-1e-6',...
    'Matsuoka-1000-1e-6','Matsuoka-500-2e-6','Matsuoka-1000-2e-6');
hold off;



