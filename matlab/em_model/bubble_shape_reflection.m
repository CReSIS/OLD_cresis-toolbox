function bubble_shape_reflection
% From: S.F. Ackley, T.E. Keliher, "Ice sheet internal radio-echo reflections
% and associated physical property changes with depth," Journal of Geophysical
% Research, Vol. 84, No. B10, (Sep. 10, 1979), pp. 5675-5680.
% -- Specifically page 5676.  We have hard copies of this article.
%
% Equations assume air bubble inclusions are prolate spheroids.
% These are ellipsoids with two of the magnitudes of the ellipsoidal axes
% equal to b and the third equal to a.

% Parameters
e1 = 1.78^2;
e2 = 1;
ab_ratio = 10;
v2 = linspace(0,0.5,21);

% Equations
Aa = -1./(ab_ratio.^2-1) + ab_ratio./(ab_ratio.^2-1).^1.5 ...
  .* log(ab_ratio+sqrt(ab_ratio.^2-1));
Ab = 0.5*(1-Aa);
ea = e1*(1 + v2*(e2-e1)/(e1+Aa*(e2-e1)));
eb = e1*(1 + v2*(e2-e1)/(e1+Ab*(e2-e1)));

eab = 1/3*ea + 2/3*eb;

es = e1*(1 + 3*v2*(e2-e1)/(2*e1 + e2));

hand1 = plot(v2,ea,'k-');
hold on;
hand2 = plot(v2,eb,'k:');
hand3 = plot(v2,eab,'kd-');
hand4 = plot(v2,es,'k+-');
hold off;
title(sprintf('Ratio a/b = %.0f',ab_ratio));
xlabel('Volume fraction of air');
ylabel('Dielectric');
legend([hand1 hand2 hand3 hand4],'e_a','e_b','e_{ab}','e_s');
pause;

% Parameters
ab_ratio = 1:32;
v2 = 0.0025;
v2 = 0.2;

% Equations
Aa = -1./(ab_ratio.^2-1) + ab_ratio./(ab_ratio.^2-1).^1.5 ...
  .* log(ab_ratio+sqrt(ab_ratio.^2-1));
Ab = 0.5*(1-Aa);
ea = e1*(1 + v2*(e2-e1)./(e1+Aa*(e2-e1)));
eb = e1*(1 + v2*(e2-e1)./(e1+Ab*(e2-e1)));

eab = 1/3*ea + 2/3*eb;

es = e1*(1 + 3*v2*(e2-e1)/(2*e1 + e2)) * ones(1,length(eab));

hand1 = plot(ab_ratio,ea,'k-');
hold on;
hand2 = plot(ab_ratio,eb,'k:');
hand3 = plot(ab_ratio,eab,'kd-');
hand4 = plot(ab_ratio,es,'k+-');
hold off;
title(sprintf('Volume fraction of air: %.4f',v2));
xlabel('Ratio a/b');
ylabel('Dielectric');
legend([hand1 hand2 hand3 hand4],'e_a','e_b','e_{ab}','e_s');
pause;

PRC1 = (1./sqrt(ea)-1./sqrt(eb))./(1./sqrt(ea)+1./sqrt(eb));
PRC2 = (1./sqrt(ea)-1./sqrt(eab))./(1./sqrt(ea)+1./sqrt(eab));
hand1 = plot(ab_ratio,20*log10(abs(PRC1)),'k-');
hold on;
hand2 = plot(ab_ratio,20*log10(abs(PRC2)),'k-');
hold off;
title(sprintf('Volume fraction of air: %.4f',v2));
xlabel('Ratio a/b');
ylabel('Power Reflection Coefficient (dB)');
legend([hand1 hand2],'PRC e_a w/ e_b','PRC e_a w/ e_{ab}');


