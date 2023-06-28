function A = nansum(A,dim)
% A = nansum(A,dim)
%
% Replacement function for Matlab nansum function since it requires a
% toolbox that people often do not have.

if ~exist('dim','var')
  sizeA = size(A);
  if any(sizeA > 1)
    dim = find(sizeA>1,1);
  else
    dim = 1;
  end
end

% Take the mean over all non-NaN
A(isnan(A)) = 0;
A = sum(A,dim);

end
