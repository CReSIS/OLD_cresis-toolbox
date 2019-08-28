function c = lcm_vector(a)
% c = lcm_vector(a)
%
% Find the least common multiple of a vector of numbers using a recursive
% algorithm based on lcm.m

if length(a) == 1
  c = a;
elseif length(a) == 2
  c = lcm(a(1),a(2));
else
  c = lcm(lcm_vector(a(1:end-1)),a(end));
end
