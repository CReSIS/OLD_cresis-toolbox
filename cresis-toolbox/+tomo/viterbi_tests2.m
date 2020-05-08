% Standalone test environment for viterbi implementation
% Author: Reece Mathews

function viterbi_tests2()
  global matrix layer;
  
  % CONSTANTS
  rows = 20;
  cols = 20;
  surf = 5;
  mult = 10;
  grnd = 17;
  matrix = zeros(rows, cols);
  for i = 1:cols
      matrix(floor(abs(sin(i)/4)*rows+1) + 8, i) = 30; % Surface
  end
  matrix(mult, :) = ones(1, cols) * 20; % Mult 1
  matrix(grnd, :) = ones(1, cols) * 10; % Bottom
  
  surf_bins = ones(1, cols) * surf;

  % Viterbi params
  along_track_weight = 1;
  along_track_slope = round(diff(surf_bins(1, :)));
  upper_bounds = ones(1, cols)*NaN;
  lower_bounds = ones(1, cols)*NaN;
  
  upper_bounds(6) = 1;
  lower_bounds(6) = 5;
  
  
  matrix = echo_norm(matrix,struct('scale',[-40 90]));
  
  % RUN
  layer = tomo.viterbi2(matrix, along_track_slope, along_track_weight, upper_bounds, lower_bounds);
  hfig = setup();
  resize(hfig);
end

function plot_viterbi()
  global layer;
  
  hold on;
  
  plot(layer, 'g');
      
  hold off;
end


function hfig = setup() 
  if ishandle(1)
    hfig = figure(1);
  else
    hfig = figure('SizeChangedFcn', @resize);
  end
end


function resize(src,~)
  global matrix;
  
  r = size(matrix, 1);
  c = size(matrix, 2);
  
  % Plot matrix
  clf;
  colormap(1-gray);
  imagesc(matrix);
      
  % X-axis
  xticks(.5:(c + .5));
  xlim([.5 (c + .5)]);
  xticklabels('');
  
  % Y-axis
  yticks(.5:(r + .5));
  ylim([.5 (r + .5)]);
  yticklabels('');
  grid;
  
  w = src.Position(3);
  h = src.Position(4);
  axes_pos = get(gca, 'Position');

  left = axes_pos(1) * w;
  bottom = axes_pos(2) * h;

  xmargin = axes_pos(3)*w/c;
  ymargin = axes_pos(4)*h/r;
  dimensions = [r c];
  for axis = 1:2
    for i = 1:dimensions(axis)
      char = uicontrol('Style', 'text', 'String', sprintf('%d | %d', i-1, i), 'FontSize', 10);
      char_pos = get(char, 'Extent');
      char_width = char_pos(3);
      char_height = char_pos(4);
      if axis == 2
        set(char, 'Position', [left + xmargin*i - xmargin/2 - char_width/2, bottom - char_height * 1.25, char_width, char_height]);
      else
        set(char, 'Position', [left - char_width * 1.25, bottom + ymargin*(r-i) + ymargin/2 - char_height/2, char_width, char_height]);
      end
    end
  end
 
  plot_viterbi();
end
