function g = h2g(h, epsilon);

% G = h2g(H,EPSILON)
%
% Hybrid-H to Hybrid-G transformation
% H and G are matrices of size [2,2,F]
% where F is the number of frequencies
% (the number of ports is always 2)
% 
% EPSILON is a limit used in finding correspondent Hybrid-H matrices in the
% vicinity of singularities; by default 1e-12, should be enough for most
% realistic problems; could be increased for a gain in speed
%
% written by tudor dima, tudima@zahoo.com, change the z into y

if nargin < 2 epsilon = 1e-12; end;
d = 1 - h(1,2,:).*h(2,1,:);
[n,i] = min(abs(d));
exact_g = 1;
while n <= epsilon
    exact_g = 0; csi = round(rand);
    h(1+csi,1+~csi,i) = h(1+csi,1+~csi,i)+(rand-0.5)*epsilon;
    d = 1 - h(1,2,:).*h(2,1,:);
    [n,i] = min(abs(d));
end;

if exact_g == 0
    fprintf(1,'%s\n%s\n', 'h2g: correspondent hybrid-G matrix non-existent', 'an approximation is produced');
end;

g(1,1,:) = h(2,2,:)./d;
g(1,2,:) = -h(1,2,:)./d;
g(2,1,:) = -h(2,1,:)./d;
g(2,2,:) = h(1,1,:)./d;
