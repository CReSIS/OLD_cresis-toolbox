function idxs = monotonic_indexes(A,strict)
% idxs = monotonic_indexes(A,strict)
%
% idxs = monotonic_indexes([1 2 3 3 4 5],true)
% idxs = monotonic_indexes([1 2 3 3 4 5],false)
% idxs = monotonic_indexes([1 2 3 2 3 4 5],true)
% idxs = monotonic_indexes([-inf 2 3 2 3 4 inf 5],true)
% idxs = monotonic_indexes([NaN -inf 2 3 2 3 nan inf 5],true)

mask = false(size(A));
cur_value = nan;
if strict
  for idx = 1:length(A)
    if A(idx)>cur_value || ~isnan(A(idx)) && isnan(cur_value)
      cur_value = A(idx);
      mask(idx) = true;
    end
  end
else
  for idx = 1:length(A)
    if A(idx)>=cur_value || ~isnan(A(idx)) && isnan(cur_value)
      cur_value = A(idx);
      mask(idx) = true;
    end
  end
end
idxs = find(mask);

