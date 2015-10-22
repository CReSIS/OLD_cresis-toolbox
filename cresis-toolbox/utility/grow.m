function B = grow(A, order)
% B = grow(A, order)
%
% Takes a logical matrix and makes every zero that is adjacent to a one
% into a one also.
%
% Author: John Paden
%
% See also: shrink

if ~exist('order','var') || isempty(order)
  order = 1;
end

B = A;
while order > 0
  A = B;
  B(:,1:end-1) = B(:,1:end-1) | A(:,2:end);
  B(:,2:end) = B(:,2:end) | A(:,1:end-1);
  B(1:end-1,:) = B(1:end-1,:) | A(2:end,:);
  B(2:end,:) = B(2:end,:) | A(1:end-1,:);
  order = order - 1;
end

return;
