function er_r = drysnow(ei,vi,roh_ds)
% Permitivity of DrySnow using the 2 eqns on page 2062 MRS
% figure 5 in the report
% 
% roh_ds=0:0.01:1;
% vi=roh_ds/0.916;
% ei=3.15;
% er_r = drysnow(ei,vi,roh_ds)
% figure(3); clf;
% f1=plot(roh_ds,er_r,'k-');
% xlabel('Snow Density (g cm^-3)');
% ylabel('Dry snow permittivity eds (real part)');
%
% rho_ds = g/cm^3 with 0.916 g/cm^3 being solid ice

N=3*vi*(ei-1);
D=(2+ei)-(vi*(ei-1));
T=N/D;
eds=1+T;
r = find(roh_ds<=0.5);
er_r(r)=(1+1.9*roh_ds(r));
r = find(roh_ds>=0.5);
er_r(r)=(0.51+2.88*roh_ds(r));
%  er_r(r)=(1+0.47*vi(r)).^3;

% Another permitivity relation using only rho_ds from:
% Ulaby, F. T., R. K. Moore, and A. K. Fung (1986), Microwave Remote Sensing: Active
% and Passive, Vol. II { Radar Remote Sensing and Surface Scattering and Emission
% Theory, Microwave Remote Sensing (Book 2), vol. 2, Artech House, 685 Canton Street,
% Norwood, MA 02062.
% er_r = (1 + 0.51*rho_ds).^3;

return;
