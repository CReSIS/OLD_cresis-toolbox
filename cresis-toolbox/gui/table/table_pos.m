function pos = table_pos(table, row, col, in_pos, height)
% pos = table_pos(table, row, col, in_pos)
%
% Convert from natural units to Matlab units

pos(1) = table.offset(1)+table.width_margin(row,col)+in_pos(1);
pos(2) = -table.offset(2)+table.height_margin(row,col)+table.size(2)-in_pos(2)-height;

return;
