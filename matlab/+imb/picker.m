% For debugging purposes, it is nice to have the handle (gpicker) available
global gpicker;
if isempty(gpicker)
  gpicker = imb.mapwin();
else
  gpicker(end+1) = imb.mapwin();
end

return;
