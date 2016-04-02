function Amed = nanmedian(A,dim)
% Amed = nanmedian(A,dim)
%
% Replacement function for Matlab nanmedian function since it requires a
% toolbox that people often do not have.

if ~exist('dim','var')
  sizeA = size(A);
  if any(sizeA > 1)
    dim = find(sizeA>1,1);
  else
    dim = 1;
  end
end

sizeAmed = size(A);
sizeAmed(dim) = 1;
Amed = zeros(sizeAmed,class(A));

A = permute(A,[dim 1:dim-1 dim+1:ndims(A)]);
sizeA = size(A);

for n = 1:prod(sizeA(2:end))
  vals = A(:,n);
  Amed(1,n) = median(vals(isnan(vals)));
end

Amed = permute(Amed,[2:dim 1 dim+1:ndims(A)]);

end
