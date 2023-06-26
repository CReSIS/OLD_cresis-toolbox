function max_val = finitemax(A)

max_val = max(A(isfinite(A)));

end
