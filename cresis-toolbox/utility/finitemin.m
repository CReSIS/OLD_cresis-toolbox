function min_val = finitemin(A)

min_val = min(A(isfinite(A)));

end
