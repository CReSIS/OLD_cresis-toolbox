%Modified Debye like model
% Returns dielectric of wet snow according to Hallikainen
% freq = frequency in GHz
% clear;
% clc;
% f=3:1:37;
% mv=2:2:12;
% roh=0.25;
% for r = 1:1:length(mv)
% [er_r(:,r),er_i(:,r)]=modified_debyelike_model(mv(r),f,roh);
% end
% figure(1); clf;
% f1=plot(f,er_r,'k-');
% title('Modified Debye-like model, snow desnity used = 0.25g/cc')
% xlabel('frequency f (GHz)');
% ylabel('Permitivity of wet snow')
% ylim([1.2,3.4]);
% figure(2); clf;
% f2=plot(f,er_i,'k');
% xlabel('frequency f (GHz)');
% ylabel('Dielectric loss factor')
% title('Modified Debye-like model, snow desnity used = 0.25g/cc')

function [er_r,er_i]=modified_debyelike_model(mv,f,roh)
fo=8.8;
A1=0.78+(0.03.*f)-(0.58e-3.*(f.^2));
A2=0.97-(0.39e-2.*f)+(0.39e-3.*(f.^2));
B1=0;
%B1= 0.31-(0.05.*f)+(0.87e-3.*(f.^2));
A=1+(1.83*roh)+(0.02.*A1.*(mv.^1.015))+B1;
B=0.073.*A1;
C=0.073.*A2;
x=1.31;
N=B*(mv^x);
D=1+((f./fo).^2);
er_r=A+(N./D);
er_i=(C.*mv^x.*(f/fo))./(1+(f/fo).^2);
return;

