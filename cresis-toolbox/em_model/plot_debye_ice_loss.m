% Plot different ice losses that have been measured
%
% This file is for illustration purposes.  Not very good...

physicalConstants;

% Plot ice loss from Johari-Charette data
joh = load('johari_charette_table.txt');

T = 273.15 + joh(:,1);
er = joh(:,2) - j*joh(:,2).*joh(:,3);
w = 2*pi*35e6;
atten35 = 20*log10(exp(1))*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
er = joh(:,4) - j*joh(:,4).*joh(:,5);
w = 2*pi*60e6;
atten60 = 20*log10(exp(1))*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));

hand1 = plot(T,atten35*1000,'bx');
xlabel('Temperature (Kelvin)');
ylabel('Attenuation (dB/km)');
hold on;
hand2 = plot(T,atten60*1000,'bo');

% Plot ice loss from Johari data
joh = load('johari_table.txt');
w = 2*pi*joh(:,1)*1e6;
T = 273.15-5;
er = joh(:,2) - j*joh(:,3);
atten = 20*log10(exp(1))*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));

T = ones(size(w))*T;
hand3 = plot(T,atten*1000,'rx');

% Plot ice loss from Westphal data
load_westphal
T = 273.15 + tuto.T;
symb = {'gd','rd','bd','yd','cd'};
for c_ind = 1:size(tuto.e_real,2)
  w = 2*pi*tuto.freq(c_ind);
  er = tuto.e_real(:,c_ind) - j*tuto.e_imag(:,c_ind);
  atten = 20*log10(exp(1))*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
  hand4(c_ind) = plot(T,atten*1000,symb{c_ind});
end

% Plot ice loss from Debye (Auty/Cole) data
auty = load('auty_ice_prop.txt');
T = 273.15 + auty(:,1);
w = 2*pi*60e6;
er = auty(:,3) + (auty(:,2)-auty(:,3)) ./ (1 + j*w*auty(:,4)/1e5);
atten = 20*log10(exp(1))*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
hand5 = plot(T,atten*1000,'k+');

% Plot ice loss from Matsuoka

T = [225:273];
er = ice(T,0.917,60e6,0,0.22);
w = 2*pi*60e6;
atten = 20*log10(exp(1))*real(j*sqrt(-j*w.*u0.*(j*w.*er*e0)));
hand6 = plot(T,atten*1000,'r-');

axis([223.15 273.15 0 60]);
title('Pure Ice Loss (MHz Data)');
legend([hand1 hand2 hand3 hand4 hand5 hand6],'Johari-35',...
    'Johari-60','Johari-0.5 thru 100',...
    'Westphal-150','Westphal-300','Westphal-500','Westphal-1000','Westphal-2700',...
    'Auty-Cole-60','Matsuoka-60');
hold off;



