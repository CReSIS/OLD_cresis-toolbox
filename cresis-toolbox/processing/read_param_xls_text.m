function [out] = read_param_xls_text(row,col,num,txt)
% [out] = read_param_xls_text(row,col,num,txt)
%
% Text Cell Read Function
%
% Support function for read_param_xls_*
%
% Author: John Paden, Brady Maasen

if size(num,1) >= row && size(num,2) >= col && ~isnan(num(row,col))
  error('Cell (%d,%d)/(%d,%s) must be text (use '' before the number)\n', row, col, row, 'A'+col-1);
end
if size(txt,1) < row || size(txt,2) < col
  out = '';
else
  out = txt{row,col};
end


return
