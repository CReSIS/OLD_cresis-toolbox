function h = y2h(y);

% H = y2h(Y)
%
% Admittance Hybrid transformation
% only for 2-by-2 matrices
%
% martie 27

if y(1,1) == 0
  disp('correspondent admittance matrix non-existent');
else
y(1,1) = 1/y(1,1);
y(1,2) = -y(1,2)/y(1,1);
y(2,1) = y(2,1)/y(1,1);
y(2,2) = y(2,2) - y(1,2)*y(2,1)/y(1,1);
end;