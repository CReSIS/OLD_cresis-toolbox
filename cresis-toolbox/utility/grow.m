function B = grow(A, order, conn)
% B = grow(A, order, conn)
%
% Takes a logical matrix and makes every zero that is adjacent to a one
% into a one also.
%
% A: logical matrix
% order: nonnegative scalar integer which specifies the number of times to
% run the grow operation (default is 1)
% conn: conn may be 4 or 8 (default is 4)
%   Same definition as bwboundaries. 4 only grows north, south, east, west.
%   8 grows northeast, southeast, northwest, southwest too.
%
% Author: John Paden
%
% See also: grow, shrink

if ~exist('order','var') || isempty(order)
  order = 1;
end
if ~exist('conn','var') || isempty(conn)
  conn = 4;
end

B = A;
while order > 0
  A = B;
  B(:,1:end-1) = B(:,1:end-1) | A(:,2:end);
  B(:,2:end) = B(:,2:end) | A(:,1:end-1);
  B(1:end-1,:) = B(1:end-1,:) | A(2:end,:);
  B(2:end,:) = B(2:end,:) | A(1:end-1,:);
  if conn == 8
    B(1:end-1,1:end-1) = B(1:end-1,1:end-1) | A(2:end,2:end);
    B(1:end-1,2:end) = B(1:end-1,2:end) | A(2:end,1:end-1);
    B(2:end,1:end-1) = B(2:end,1:end-1) | A(1:end-1,2:end);
    B(2:end,2:end) = B(2:end,2:end) | A(1:end-1,1:end-1);
  end
  order = order - 1;
end

return;
