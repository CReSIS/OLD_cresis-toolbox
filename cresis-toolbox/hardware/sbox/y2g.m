function g = y2g(y);

% G = y2g(Y)
%
% Admittance to Hybrid G transformation
% only for 2-by-2 matrices
%
% martie 26

if y(2,2) == 0
  disp('correspondent admittance matrix non-existent');
else
g(1,1) = y(1,1) - y(1,2)*y(2,1)/y(2,2);
g(1,2) = y(1,2)/y(2,2);
g(2,1) = -y(2,1)/y(2,2);
g(2,2) = 1/y(2,2);
end;