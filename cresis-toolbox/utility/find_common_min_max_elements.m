function [min_idx max_idx] = find_common_min_max_elements(A)
% [min_idx max_idx] = find_common_min_max_elements(A)
%
% A = cell vector of 1-D vectors, vectors can be of different lengths
%   VECTORS MUST BE SORTED
% min_idx: index of the smallest element that is contained in everyone
%   of the vectors in A (index is relative to A{1})
% max_idx: index of the largest element that is contained in everyone
%   of the vectors in A (index is relative to A{1})
%
% Example:
% A = {};
% A{1} = [1 4 5 6 8];
% A{2} = [2 4 5 6 8];
% A{3} = [1 3 6 8];
% [min_idx max_idx] = find_common_min_max_elements(A);
% A{1}(min_idx)
% A{1}(max_idx)
%
% Author: John Paden

min_idx = 1;
while min_idx <= length(A{1})
  % Look for this element in all the other vectors
  found = true;
  for vec_idx = 2:length(A)
    found_idx = find(A{vec_idx} <= A{1}(min_idx),1,'last');
    if isempty(found_idx) || A{vec_idx}(found_idx) ~= A{1}(min_idx)
      found = false;
      break;
    end
  end
  if found
    break;
  else
    min_idx = min_idx + 1;
  end
end
if ~found
  error('No common elements found');
end

max_idx = length(A{1});
while max_idx >= 1
  % Look for this element in all the other vectors
  found = true;
  for vec_idx = 2:length(A)
    found_idx = find(A{vec_idx} >= A{1}(max_idx),1);
    if isempty(found_idx) || A{vec_idx}(found_idx) ~= A{1}(max_idx)
      found = false;
      break;
    end
  end
  if found
    break;
  else
    max_idx = max_idx - 1;
  end
end

return;

