function [output_dir,radar_type,radar_name] = ct_output_dir(radar_name)
% [output_dir,radar_type,radar_name] = ct_output_dir(radar_name)
%
% Returns output directory, general radar type and name for each radar
% when given the specific radar_name. Used with the parameter spreadsheets
% param.radar_name field.
%
% radar_name: Specific radar string (e.g. accum2, snow8, mcords5, mcords5-accum)
%
% output_dir: Standard output directory (e.g. accum, snow, rds, accum)
% radar_type: 'deramp', 'stepped', or 'pulsed'
% radar_name: Digital format radar name (e.g. (accum, snow, mcords5,
%   mcords5)
%
% Example:
%   param.radar_name = 'mcords5-accum';
%   [output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name)
%
% Author: John Paden

if ~isempty(find(radar_name=='-'))
  [radar_name,output_dir_override] = strtok(radar_name,'-');
  output_dir_override = output_dir_override(2:end);
else
  output_dir_override = '';
end

if any(strcmpi(radar_name,{'kuband','kuband2','kuband3'}))
  radar_type = 'deramp';
  output_dir = 'kuband';
elseif any(strcmpi(radar_name,{'kaband3'}))
  radar_type = 'deramp';
  output_dir = 'kaband';
elseif any(strcmpi(radar_name,{'snow','snow2','snow3','snow5','snow8','snow9'}))
  radar_type = 'deramp';
  output_dir = 'snow';
elseif any(strcmpi(radar_name,{'accum0'}))
  radar_type = 'deramp';
  output_dir = 'accum';
elseif any(strcmpi(radar_name,{'accum'}))
  radar_type = 'stepped';
  output_dir = 'accum';
elseif any(strcmpi(radar_name,{'accum2','accum3'}))
  radar_type = 'pulsed';
  output_dir = 'accum';
elseif any(strcmpi(radar_name,{'rds','cords','acords','hfrds','hfrds2','icards','mcrds','mcords','mcords2','mcords3','mcords4','mcords5','wise'}))  radar_type = 'pulsed';
  output_dir = 'rds';
else
  error('Invalid radar type %s\n', radar_name);
end

if ~isempty(output_dir_override)
  output_dir = output_dir_override;
end

return;

