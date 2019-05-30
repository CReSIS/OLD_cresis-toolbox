function Astd = nanstd(A,ignore,dim)
% Astd = nanstd(A,ignore,dim)
%
% Replacement function for Matlab nanstd function since it requires a
% toolbox that people often do not have.

if ~exist('dim','var')
  sizeA = size(A);
  if any(sizeA > 1)
    dim = find(sizeA>1,1);
  else
    dim = 1;
  end
end

A = permute(A,[dim 1:dim-1 dim+1:ndims(A)]);
sizeA = size(A);

sizeAstd = size(A);
sizeAstd(1) = 1;
Astd = zeros(sizeAstd,class(A));

for n = 1:prod(sizeA(2:end))
  vals = A(:,n);
  Astd(1,n) = std(vals(~isnan(vals)));
end

Astd = permute(Astd,[2:dim 1 dim+1:ndims(A)]);

end
