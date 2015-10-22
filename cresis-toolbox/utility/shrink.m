function B = shrink(A, order)
% B = shrink(A, order)
%
% Takes a logical matrix and makes every one that is adjacent to a zero
% into a zero also.
%
% Author: John Paden
%
% See also: grow

if ~exist('order','var') || isempty(order)
  order = 1;
end

B = A;
while order > 0
  A = B;
  B(:,1:end-1) = B(:,1:end-1) & A(:,2:end);
  B(:,2:end) = B(:,2:end) & A(:,1:end-1);
  B(1:end-1,:) = B(1:end-1,:) & A(2:end,:);
  B(2:end,:) = B(2:end,:) & A(1:end-1,:);
  order = order - 1;
end

return;
