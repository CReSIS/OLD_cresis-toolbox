function mask = sidelobe_mask(data, sl_rows, sl)
% mask = sidelobe_mask(data, sl_rows, sl)
%
% A = single(randn(10));
% mask = sidelobe_mask_mex(A,int32([-2:2]),single(ones(1,5)))
% mask = sidelobe_mask(A,int32([-2:2]),single(ones(1,5)))
%
% Author: John Paden

mask = false(size(data));
num_rows = size(data,1);
num_cols = size(data,2);

for col = 1:num_cols
  for row = 1:num_rows
    for sl_rows_idx = 1:length(sl_rows)
      cur_row = row+sl_rows(sl_rows_idx);
      if cur_row >= 1 && cur_row <= num_rows && data(cur_row,col)+sl(sl_rows_idx) < data(row,col)
        mask(cur_row,col) = true;
      end
    end
  end
end
