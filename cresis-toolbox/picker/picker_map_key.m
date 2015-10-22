function picker_map_key(src,event)

% Check to make sure that a key was pressed and not
% just a modifier (e.g. shift, ctrl, alt)
if ~isempty(event.Character)

  global hui; % hui: user interface handles
  global gCtrl; % Geobase: contains all the geographic info

  [x,y,but] = get_mouse_info(hui.mapfig.handle,hui.mapfig.axes.handle);
  fprintf('Map: x = %.2f, y = %.2f, but = %d\n', x, y, but);

  if 1 % This code for debugging
    if ischar(event.Key)
      fprintf('x = %.2f, y = %.2f, key = %s %d\n', x, y, event.Key, event.Character);
    else
      fprintf('x = %.2f, y = %.2f, key = %d %d\n', x, y, event.Key, event.Character);
    end
    if ~isempty(event.Modifier)
      fprintf('  Modifiers ');
      for ind = 1:length(event.Modifier)
        fprintf('%s ', event.Modifier{ind});
      end
      fprintf('\n');
    end
  end

  % see event.Modifier for modifiers
  switch event.Character
    case 8 % Backspace
    case 127 % Delete
    case 31 % Down-arrow
    case 30 % Up-arrow
    case {'z','Z'} % zoom mode
  end
end
