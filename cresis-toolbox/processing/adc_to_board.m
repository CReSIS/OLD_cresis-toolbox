function [board,board_idx] = adc_to_board(radar_name,adc)
% [board,board_idx] = adc_to_board(radar_name,adc)
%
% Support function for determining which file grouping or "board"
% a particular ADC is associated with.  E.g. used with mcords2.
%
% Inputs:
% radar_name: is string containing radar name
% adc: adc is one indexed
%
% Outputs:
% board: can be one or zero indexed
% board_idx: is one indexed
%
% Author: John Paden

[output_dir,radar_type,radar_name] = ct_output_dir(radar_name);

if any(strcmpi(radar_name,{'mcords2','mcords3'}))
  board = unique(floor((adc-1)/4));
  board_idx = board+1;
else
  board = unique(adc);
  board_idx = board;
end

return;
