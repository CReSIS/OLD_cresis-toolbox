

%% 1. Adapt number of snapshots so that flat earth distance traveled covers
% theta = 0.5 deg (or some value that makes sense from looking at the MUSIC
% imagery)

% Not sure if you will have time to do this, so later on...

%% 2. Set the standard deviation according to the expected max elevation
% change from the last bin as we discussed. This max elevation angle change
% would be found from looking at the variation in theta in the MUSIC
% imagery and trying to tie this to the geometric calculations that we
% discussed.
%
% This theta_std should be related to K/(R*sin(theta)) and not to the
% transition model. The following shows geometric calculations.
theta_std = 1;
figure(1); clf;
theta = linspace(0,89,101)/180*pi;
h = 20; % Typical height deviation
er_ice = 3.15;
T = 2000; % ice thickness
H_air = 500;
theta_ice = asin(sin(theta)/sqrt(er_ice));
R = H_air./cos(theta) + T./cos(theta_ice);
bound = pi/2 - theta_ice - asin( (R-h./cos(theta_ice))./R.*sin(pi/2+theta_ice));
plot(theta*180/pi, bound*180/pi)
hold on;
plot(theta*180/pi, h./(R.*sin(theta_ice)) * 180/pi);
ylim([0 30]);
grid on;
xlabel('Incidence angle (deg)');
ylabel('Upper bound (deg)');
title('Two geometric ways to set the upper bound or std of Gaussian')
legend('Exact','Plane wave')

%% 3. Set the prior pdf using something like this:

theta = linspace(-5,5,101);

f_prior = zeros(size(theta));

% W: widening factor beyond the standard deviation
W = 2;
% C: compression factor for negative values
C = 2;
f_prior(C*theta<-3*theta_std) = 0;
% mm: Logical mask
mm = C*theta>=-3*theta_std & theta<0; 
f_prior(mm) = exp(-(C*theta(mm)).^2/(theta_std*W).^2);
% mm: Logical mask
mm = theta >= 0 & theta>-3*theta_std;
f_prior(mm) = exp(-theta(mm).^2/(theta_std*W).^2);
f_prior(theta>3*theta_std) = 0;

figure(2); clf;
plot(theta, f_prior);
grid on;
xlabel('\theta (deg)');
ylabel('f_{prior}(\theta)');
title('Prior pdf centered on mean from transition model');


%% 3. Alternative to geometric bound
% This is the fixed 3*dtheta bound that you proposed, but I think it should
% be continuous and monotonically increasing so here is one option that is
% also continuously differentiable
theta = linspace(0,5,101);
K = 3;
theta_transition_deg = 0.8;
transition_speed = 3.5;

figure(3); clf;
plot(theta, theta.*(1 + (K-1)./(1+exp((theta-theta_transition_deg)*transition_speed)) ));
grid on; xlabel('\theta'); ylabel('Bounds (deg)'); hold on;
theta3 = theta;

K = 3;
theta_transition_deg = 0.5;
theta3(theta<=theta_transition_deg) = K*theta(theta<=theta_transition_deg);
plot(theta,theta3);
title('Upper bound replacement for 3*theta (continuously differentiable)');
legend('Continous','Original');

% Here is another option that just satisfies continuous and monotonically
% increasing
theta = linspace(0,5,101);
K = 3;
theta_transition_deg = 0.5;
transition_arc_deg = 1.5;
bound = zeros(size(theta));
mm = theta<theta_transition_deg;
bound(mm) = theta(mm)*K;
mm = theta>=theta_transition_deg & theta<=theta_transition_deg+transition_arc_deg;
bound(mm) ...
  = linspace(K*theta(find(mm,1)),theta_transition_deg+transition_arc_deg,sum(mm));
mm = theta>theta_transition_deg+transition_arc_deg;
bound(mm) = theta(mm);

figure(4); clf;
plot(theta, bound);
grid on; xlabel('\theta'); ylabel('Bounds (deg)'); hold on;
theta3 = theta;
theta3(theta<=theta_transition_deg) = K*theta(theta<=theta_transition_deg);
plot(theta,theta3);
title('Upper bound replacement for 3*theta (monotonically increasing)');
legend('Continous','Original');

