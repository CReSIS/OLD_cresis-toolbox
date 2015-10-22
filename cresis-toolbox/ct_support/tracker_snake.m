function layer = tracker_snake(A,dataPnts,offset_abs)
% A = image
% dataPnts = structure containing information about the manually
%   selected datapoints.
%   dataPnts.row
%   dataPnts.col
%   dataPnts.method: 's' == snake, 't' == threshold, 'c' == constant
%   dataPnts.snake.search_range
%   dataPnts.CFAR.noise_range
%   dataPnts.start = optional (all points after this not tracked)
%   dataPnts.stop = optional (all points after this not tracked)
% layer = output layer
%
% Layer tracking called from tracker_snake_manual_gui and can be
% called directly.

layer = zeros(1,size(A,2));
method = zeros(1,size(A,2));
if numel(dataPnts) > 0
  
  %% Prep Inputs
  if ~exist('offset_abs') || isempty(offset_abs)
    offset_abs = 0;
  end
  default_search_range = -3:3;
  noise_range = -15;
  threshold = 10;
  for ind = 1:length(dataPnts)
    row = dataPnts(ind).row;
    if ~isfield(dataPnts(ind),'snake')
      search_range = default_search_range;
    else
      search_range = dataPnts(ind).snake.search_range;
    end
    if dataPnts(ind).method == 'c'
      layer(1:end) = dataPnts(ind).row;
      return;
    end
    if row < 1-search_range(1)
      row = 1-search_range(1);
    elseif row > size(A,1)-search_range(end)
      row = size(A,1)-search_range(end);
    end
    startCol = dataPnts(ind).col - offset_abs;
    if startCol < 1
      startCol = 1;
    elseif startCol > size(A,2)
      startCol = size(A,2);
    end
    stopCol = dataPnts(ind).col + offset_abs;
    if stopCol < 1
      stopCol = 1;
    elseif stopCol > size(A,2)
      stopCol = size(A,2);
    end
    method(startCol:stopCol) = dataPnts(ind).method;
    %dataPnts(ind).threshold = A(row,dataPnts(ind).col) - threshold;
    if isfield(dataPnts(ind),'start') && dataPnts(ind).start
      layer(1:stopCol) = row;
    elseif isfield(dataPnts(ind),'stop') && dataPnts(ind).stop
      layer(startCol:end) = row;
    else
      layer(startCol:stopCol) = row;
    end
  end
  
  %% Apply the basic snake
  while sum(layer ~= 0) < size(A,2)
    for ind = 1:length(dataPnts)
      col = dataPnts(ind).col+offset_abs;
      if col >= 1 && col <= size(A,2) && layer(col) == 0
        if dataPnts(ind).method == 's'
          [val row] = max(A(layer(col-1)+search_range,col));
          row = layer(col-1)+search_range(1) - 1 + row;
        elseif dataPnts(ind).method == 't'
          row = track_layer_threshold(A(layer(col-1)+search_range,col), ...
            dataPnts(ind).thresh.value, dataPnts(ind).thresh.num_lower);
          if isempty(row)
            % Keep old value
            row = layer(col-1);
          else
            row = layer(col-1)+search_range(1) - 1 + row;
          end
          %           row = find(A(layer(col-1)+search_range,col) > dataPnts(ind).thresh.value,1);
          %           if isempty(row)
          %             [val row] = max(A(layer(col-1)+search_range,col));
          %           end
          %           row = layer(col-1)+search_range(1) - 1 + row
        end
        if row < 1-search_range(1)
          row = 1-search_range(1);
        elseif row > size(A,1)-search_range(end)
          row = size(A,1)-search_range(end);
        end
        layer(col) = row;
        method(col) = dataPnts(ind).method;
      end
    end
    
    for ind = 1:length(dataPnts)
      col = dataPnts(ind).col-offset_abs;
      if col >= 1 && col <= size(A,2) && layer(col) == 0
        if dataPnts(ind).method == 's'
          [val row] = max(A(layer(col+1)+search_range,col));
          row = layer(col+1)+search_range(1) - 1 + row;
        elseif dataPnts(ind).method == 't'
          row = track_layer_threshold(A(layer(col+1)+search_range,col), ...
            dataPnts(ind).thresh.value, dataPnts(ind).thresh.num_lower);
          if isempty(row)
            % Keep old value
            row = layer(col+1);
          else
            row = layer(col+1)+search_range(1) - 1 + row;
          end
        end
        if row < 1-search_range(1)
          row = 1-search_range(1);
        elseif row > size(A,1)-search_range(end)
          row = size(A,1)-search_range(end);
        end
        layer(col) = row;
        method(col) = dataPnts(ind).method;
      end
    end
    offset_abs = offset_abs + 1;
  end
end

return

function row = track_layer_threshold(A_rline,threshold,num_lower)

row = [];

% Find the leading edge
leading_edge_state = 0;
for bin = 1:length(A_rline)-num_lower
  if leading_edge_state == 1
    val_lower = true;
    for bin_low = bin:bin+num_lower-1
      if A_rline(bin_low) >= A_rline(bin-1)
        val_lower = false;
        break;
      end
    end
    if val_lower
      row = bin - 1;
      break;
    end
  end
  if leading_edge_state == 0 && A_rline(bin) > threshold
    leading_edge_state = 1;
  end
end

return
