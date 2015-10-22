%Returns dielectric of wet snow according to Hallikainen
%freq = frequency in Hz
%mv = % of liquid water content
% roh = density in g/cc

% mv=0:0.1:12;
% f=6e9;
% roh=0.3;
% [er_r,er_i]=Debye(mv,freq,roh);
% delta=er_r-1-1.832*roh;
% figure(1); clf;
% f1=plot(mv,delta,'k-');
% ylim([0,1.6]);
% xlim([0.0,13]);
% title('Debye-like model')
% xlabel('Liquid water content, mv(%)');
% ylabel('Incrementa dielectric constant')
% figure(2); clf;
% f2=plot(mv,er_i,'k-');
% ylim([0,1.1]);
% xlim([0.0,13]);
% xlabel('Liquid water content, mv(%)');
% ylabel('Dielectric loss factor')
% title('Debye-like model')
%%Debye like model
function [er_r,er_i]=Debye(mv,freq,roh)
fo=9.07e9;
A1=1;
B1=0;
A2=1;
A=1+(1.83*roh)+(0.02*A1.*(mv.^1.015))+B1;
B=0.073*A1;
C=0.073*A2;
er_r=A+B*mv.^1.31/(1+(freq/fo)^2);
er_i=C*mv.^1.31*(freq/fo)/(1+(freq/fo)^2);
return;
