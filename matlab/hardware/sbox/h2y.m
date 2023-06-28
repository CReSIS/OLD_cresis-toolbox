function y = h2y(h);

% Y = h2y(H)
%
% Hybrid to Admittance transformation
% only for 2-by-2 matrices
%
% martie 27

if h(1,1) == 0
  disp('correspondent admittance matrix non-existent');
else
y(1,1) = 1/h(1,1);
y(1,2) = -h(1,2)/h(1,1);
y(2,1) = h(2,1)/h(1,1);
y(2,2) = h(2,2) - h(1,2)*h(2,1)/h(1,1);
end;