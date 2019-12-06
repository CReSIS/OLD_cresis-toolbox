function y = medfilt1(x,n,max_isnan)
% y = medfilt1(x,n)
%
% Replacement of Matlab Signal Processing Toolbox medfilt1 function.
%
% Better handling of edges.
% Handles NaN.
% Only operates on first dimension, except for row vectors.
% Only allows n to be odd.

if nargin < 3 || isempty(max_isnan)
  max_isnan = n-1;
end

size_x = size(x);

y = zeros(size_x);

if size(x,1) == 1
  y = medfilt1_1D(x.',n,max_isnan).';
else
  for col = 1:prod(size_x(2:end))
    y(:,col) = medfilt1_1D(x(:,col),n,max_isnan);
  end
end

end

function x = medfilt1_1D(x,n,max_isnan)

N = length(x);

% Since multiple calls to median would be very slow for long vectors, we
% construct a matrix of all the vectors that we need to take a median of
% and call median one time.
xx = nan(n,N);
num_isnan = zeros(1,N);
rel_n = -(n-1)/2 : (n-1)/2;
for k=1:n
  if rel_n(k) < 0
    xx(k,1-rel_n(k):end) = x(1:end+rel_n(k));
  elseif rel_n(k) > 0
    xx(k,1:end-rel_n(k)) = x(1+rel_n(k):end);
  else
    xx(k,:) = x;
  end
  num_isnan = num_isnan + isnan(xx(k,:));
end
xx = sort(xx);
row_idxs = (n+1)/2 - round(num_isnan/2);
row_idxs(num_isnan>max_isnan) = n;
idxs = row_idxs + (0:N-1)*n;
x(:) = xx(idxs);
idxs = (n+1)/2 - round((num_isnan-0.5)/2) + (0:N-1)*n;
x(:) = x(:) + xx(idxs).';
x = x/2;

end
