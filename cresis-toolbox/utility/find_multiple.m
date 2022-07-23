function [varargout] = find_multiple(varargin)

% function [varargout] = find_multiple(varargin)
% [I,J,V] = find_multiple(X, matrix_cond, R, find_str_extra)
%
% Author: Hara Madhav Talasila

switch nargin
  case 3
    X = varargin{1};
    matrix_cond = varargin{2};
    R = varargin{3};
    find_str_extra = [];
  case 5
    X = varargin{1};
    matrix_cond = varargin{2};
    R = varargin{3};
    N_vals = varargin{4};
    find_str_extra = varargin{5};
  otherwise
    warning('find_multiple: argin not supported !!!');
    return;
end
    
% find for multiple values in a vector "vals" (row/column) in "X" matrix
% [I,J,V] = find_multiple(varargin)
% add varargin support later
X_len = numel(X);
R_len = length(R);

I = [];
J = [];
V = [];


for idx = 1:length(R)
  if isempty(find_str_extra)
    eval_this = sprintf('[I{%d}, J{%d}, V{%d}] = find(X %s R(%d)); ', idx, idx, idx, matrix_cond, idx);
  else
    eval_this = sprintf('[I{%d}, J{%d}, V{%d}] = find(X %s R(%d), %d, ''%s''); ', idx, idx, idx, matrix_cond, idx, N_vals, find_str_extra);
  end
  eval(eval_this);
end

switch nargout
  case {0,1}
    varargout{1} = J;
  case 2
    varargout{1} = I;
    varargout{2} = J;
  case 3
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = V;
  otherwise
    warning('find_multiple: argout not supported !!!');
end
    
return;

%% test

a = reshape( 1:9, 3,3)';
a = 1:10;

b = [3 4];
[M N O] = find(a>=b)

[I J V] = find_multiple(a,'>=',b)

[I J V] = find_multiple(a,'>=',b, 1,'first')

tt = find_multiple(a,'>=',b, 1,'first')

find_multiple(a,'>=',b, 1,'first')