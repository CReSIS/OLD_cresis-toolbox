function makeIcePropTable.m
% Makes tables for prop_ice.doc

% Gas constant Joules per mole per kelvin (from internet)
R = 8.3143;

ice_prop = load('auty_ice_prop.txt');
T = ice_prop(:,1);

% Make the time-constant (Debye relaxation time) table
% e = e_inf + (e_stat-e_inf)/(1+j*w*tao)

% Table and values are from Auty and Cole 1952
% tao = A*exp(B/(R*T)), A = 5.3e-16 sec, B = 13200 calorie/mole
A = 5.3e-16;
B = 13200;
% Conversion from calories to joules (from internet)
B = B * 4.19;

tao = A*exp(B./(R*(T+273.15)));

hand1 = semilogy(T+273.15,ice_prop(:,4)/1e5,'r');
hold on;
hand2 = semilogy(T+273.15,tao,'b');
hold off;
xlabel('Temperature (Kelvin)');
ylabel('Debye time-constant (sec/rad)');
legend([hand1 hand2], 'Auty Table','Function');
title('Comparison of table and functions (small unexplained error)');

pause;
clear all;
clf;

ice_prop = load('matsuoka_ice_prop.txt');
T = ice_prop(:,1);

% Gas constant Joules per mole per kelvin (from internet)
R = 8.3143;
% T0 (Curie-Weiss temperature) and beta from Kawada 1978 & Matsuoka 1996
T0 = 15;
beta = 23700;
% From Matsuoka 1996
tao0 = 5.3e-16;

% Activation energy (from Matsuoka 1996)
E = [T>=223]*55300 + [T<223]*22600;

% Resonant frequency (correcting sign error in exp)
tao = tao0*exp(E./(R*T));
% Correction to match data... not sure why???
tao(find(T<223)) = tao(find(T<223))/2e-8
res_w = 1./tao;
res_freq = res_w/(2*pi);
% res_freq = (2*pi*tao0*exp(E./(R*T))).^-1

% Resonant frequency using tao table from Auty and Cole 1952
auty = load('auty_ice_prop.txt');
tao_auty = exp(interp1(auty(:,1),log(auty(:,4)/1e5),T-273.15,'spline','extrap'));
res_freq_auty = 1./tao_auty/(2*pi);

% From Matsuoka 1996 (including /1e9 correction)
A_d = beta./(T-T0).*res_freq/1e9;

% A_d_auty = beta./(T-T0).*res_freq_auty/1e9;
A_d_auty = interp1(auty(:,1),auty(:,2)-auty(:,3),T-273.15,'spline','extrap') ...
   .*res_freq_auty/1e9;

hand1 = semilogy(T,ice_prop(:,2)/1e4,'r');
hold on;
hand2 = semilogy(T,A_d,'b');
hand3 = semilogy(T,A_d_auty,'g');
hold off;
xlabel('Temperature (Kelvin)');
ylabel('A');
legend([hand1 hand2 hand3], 'Table','Function','Auty Table');
title('Comparison of function and tables (small unexplained error)');

return;
