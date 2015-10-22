function [coeff_output, A] = detrending_coeff(win_len, order)
%win_len = 51;
%order = 2;
n = (win_len - 1) / 2;
A = ones(win_len, order + 1); 
x = -n : n;
for j = 0 : order
    A(:, j + 1) = x.^j;
end

coeff_output = inv(A' * A) * A'; 