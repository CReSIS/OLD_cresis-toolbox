function default = default_radar_params_settings_match(defaults,settings)
% default = default_radar_params_settings_match(defaults,settings)
%
% Function for matching settings parameters with default parameters
%
% Author: John Paden

found = false;
if isfield(settings,'XML_File_Path')
  settings_fn = settings.XML_File_Path{1}.values{1};
  for default_idx = 1:length(defaults)
    if ~isempty(regexp(settings_fn, defaults{default_idx}.config_regexp))
      default = defaults{default_idx};
      found = true;
      break;
    end
  end
else
  settings_fn = '';
end

while ~found
  warning('Did not find a matching set of default parameters for %s.', settings_fn);
  for default_idx = 1:length(defaults)
    fprintf(' (%d) %s\n', default_idx, defaults{default_idx}.name);
  end
  default_idx = input('Choose default param to associate with this file: ');
  if length(default_idx) == 1
    try
      default = defaults{default_idx};
      found = true;
    end
  end
end

return