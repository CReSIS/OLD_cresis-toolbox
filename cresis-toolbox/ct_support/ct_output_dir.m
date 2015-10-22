function [output_dir,radar_type] = ct_output_dir(radar_name)

if any(strcmpi(radar_name,{'kuband','kuband2','kuband3'}))
  radar_type = 'fmcw';
  output_dir = 'kuband';
elseif any(strcmpi(radar_name,{'kaband3'}))
  radar_type = 'fmcw';
  output_dir = 'kaband';
elseif any(strcmpi(radar_name,{'snow','snow2','snow3','snow5'}))
  radar_type = 'fmcw';
  output_dir = 'snow';
elseif any(strcmpi(radar_name,{'accum0'}))
  radar_type = 'fmcw';
  output_dir = 'accum';
elseif any(strcmpi(radar_name,{'accum'}))
  radar_type = 'stepped';
  output_dir = 'accum';
elseif any(strcmpi(radar_name,{'accum2'}))
  radar_type = 'pulsed';
  output_dir = 'accum';
elseif any(strcmpi(radar_name,{'rds','cords','icards','acords','mcrds','mcords','mcords2','mcords3','mcords4','mcords5','wise'}))
  radar_type = 'pulsed';
  output_dir = 'rds';
else
  error('Invalid radar type %s\n', radar_name);
end

return;

