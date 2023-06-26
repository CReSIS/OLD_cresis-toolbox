function y = g2y(g);

% Y = g2y(G)
%
% Hybrid G to Admittance transformation
% only for 2-by-2 matrices
%
% martie 26

if g(2,2) == 0
  disp('correspondent admittance matrix non-existent');
else
y(1,1) = -g(2,1)/g(2,2);
y(1,2) = 1/g(2,2);
y(2,1) = g(1,1) - g(1,2)*g(2,1)/g(2,2);
y(2,2) = g(1,2)/g(2,2);
end;