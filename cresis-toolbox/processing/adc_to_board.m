function board = adc_to_board(radar_name,adc)
% board = adc_to_board(radar_name,adc)
%
% Support function for determining which file grouping or "board"
% a particular ADC is associated with.  E.g. used with mcords2.
%
% radar_name is string containing radar name
% adc is one indexed
% board is zero indexed
%
% Author: John Paden

if any(strcmpi(radar_name,{'mcords2','mcords3'}))
  board = unique(floor((adc-1)/4));
else
  board = unique(adc);
end

return;
