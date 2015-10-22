function [out] = read_param_xls_boolean(row,col,num,txt)
% [out] = read_param_xls_boolean(row,col,num,txt)
%
% Boolean Cell Read Function
%
% Support function for read_param_xls_*
%
% Author: John Paden, Brady Maasen

if size(num,1) < row || size(num,2) < col || isnan(num(row,col))
  if size(txt,1) < row || size(txt,2) < col || isempty(txt{row,col})
    out = false;
  else
    out = logical(eval(txt{row,col}));
  end
else
  out = logical(num(row,col));
end

return
