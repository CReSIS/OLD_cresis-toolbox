function [c,ceq] = doa_nonlcon(theta,theta_guard)
% [c,ceq] = sim.doa_nonlcon(theta,theta_guard)
%
% Nonlinear constraint for DOA optimization. To be used with fmincon.m
%
% theta = vector of theta estimates (radians)
% theta_guard = minimum distance between sources when optimizing
%
% c = nonlinear constraint that must satisfy <=0
% ceq = nonlinear constraint that must satisfy == 0
%
% Author: John Paden, Theresa Stumpf
%
% See also: doa.m

c = -1*ones(size(theta));
if length(theta) >= 2
  % If more than or equal to 2 sources, then enforce nonlinear constraint
  for theta_idx = 1:length(theta)
    c(theta_idx) = -(min(abs(theta(theta_idx) - theta([1:theta_idx-1 theta_idx+1:end]))) - theta_guard);
  end
end

ceq = zeros(size(theta));

end
