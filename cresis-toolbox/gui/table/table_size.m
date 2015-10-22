function Tsize = table_size(table)
% Tsize = table_size(table)
%
% Size of table's user interface control. Called by table_draw.m.

set(table.ui,'ResizeFcn',[]);
set(table.ui,'Units','Points');
tmp = get(table.ui,'Position');
Tsize = tmp(3:4);

% Strange and unexplained correction factor:
%Tsize(2) = Tsize(2) * get(0,'ScreenPixelsPerInch')/72;

set(table.ui,'Units','Normalized');
set(table.ui,'ResizeFcn',@table_resize);


return;