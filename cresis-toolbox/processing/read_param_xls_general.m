function [out] = read_param_xls_general(row,col,num,txt)
% [out] = read_param_xls_general(row,col,num,txt)
% 
% General Cell Read Function
%
% Support function for read_param_xls_*
%
% Author: John Paden, Brady Maasen

if size(num,1) < row || size(num,2) < col
  if size(txt,1) >= row && size(txt,2) >= col && ~isempty(txt{row,col})
    out = eval(txt{row,col});
  else
    out = [];
  end
elseif isnan(num(row,col))
  if size(txt,1) < row || size(txt,2) < col || isempty(txt{row,col})
    out = [];
  else
    out = eval(txt{row,col});
  end
else
  out = num(row,col);
end


return
