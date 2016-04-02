function Amean = nanmean(A,dim)
% Amean = nanmean(A,dim)
%
% Replacement function for Matlab nanmean function since it requires a
% toolbox that people often do not have.

if ~exist('dim','var')
  sizeA = size(A);
  if any(sizeA > 1)
    dim = find(sizeA>1,1);
  else
    dim = 1;
  end
end

% Count the number of non-NaN values
Nnonnan = sum(~isnan(A),dim);
Nnonnan(Nnonnan==0) = NaN;

% Take the mean over all non-NaN
A(isnan(A)) = 0;
Amean = sum(A,dim) ./ Nnonnan;

end
