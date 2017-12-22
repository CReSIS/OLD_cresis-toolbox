function [px_out,py_out,idxs] = remove_intersections(px,py,buffer,single_polygon_flag)
% [px_out,py_out,idxs] = tomo.remove_intersections(px,py,buffer,single_polygon_flag)
%
% Useful for gridding data that comes from self-intersecting swaths. This
% finds the bounding polygon with no self-intersections. This code is not
% tested very well and may have bugs. There is also an assumption/HACK of
% removing intersecting segments with any matching coordinates.
%
% px,py: polygon (may or may not have self intersections)
% buffer: HACK with inpolygon which allows us to check if a point is
%   inside of a polygon and not just on the edge. inpolygon normally returns
%   true for any point on the edge of the polygon. We do this by checking
%   (roughly) that all points within buffer distance from the point are inside
%   the polygon. If they are, we assume that the point is in the polygon and
%   not just on the edge of the polygon.
% single_polygon_flag: join all polygons into one single polygon (joins 
%   the polygons at the closest points between each polygon and ensures
%   no intersections).
%
% px_out,py_out: polygon without self intersections (column vectors)
% idxs: indexes to the original points for each entry in the polygon
%   For polygon breaks, a NaN is inserted. px_out can be constructed from
%   px and idxs:
%    px_out = zeros(size(idxs));
%    px_out(isnan(idxs)) = NaN;
%    px_out(~isnan(idxs)) = px(idxs);
%
% Example:
%   [px_out,py_out,idxs] = tomo.remove_intersections(px,py,0.1,true);
%    px_check = zeros(size(idxs));
%    px_check(isnan(idxs)) = NaN;
%    px_check(~isnan(idxs)) = px(idxs(~isnan(idxs)));
%    any(px_out(~isnan(px_out)) ~= px_check(~isnan(px_out)));
%
% Author: John Paden

% debug_state_check: Set to true to enable state check debugging
debug_state_check = 0;
% debug_plots: Set to true to enable plots for debugging
debug_plots = 0;
% debug_poly: Set to true to enable polygon debugging
debug_poly = 0;

% Find all the segments that intersect
%   (each row of segs contains a pair of indexes to crossing segments)
[x0,y0,segs] = tomo.selfintersect(px,py);

% Remove segments that include overlapping points (i.e. if either end of
% the segment has the same position as either end of the intersecting
% segment, remove that segment). HACK: The assumption here is that
% any overlapping point means that the segment overlaps another segment...
% a reasonable assumption for single/double precision 3D surfaces.
overlap_mask = logical(zeros(size(segs,1),1));
for seg_idx=1:size(segs,1)
  if any(px(segs(seg_idx,1)+[0 1 0 1]) == px(segs(seg_idx,2)+[0 1 1 0]) ...
      & py(segs(seg_idx,1)+[0 1 0 1]) == py(segs(seg_idx,2)+[0 1 1 0]))
    overlap_mask(seg_idx) = true;
  end
end
segs = segs(~overlap_mask,:);

% Plots for debugging
if 0
  clf;
  plot(px,py);
  hold on
  for row = 1:size(segs,1)
    % Draw segment
    h_plot = plot(px(segs(row,1)+[0 1]),py(segs(row,1)+[0 1]),'-','LineWidth',2);
    % Draw intersecting segment
    plot(px(segs(row,2)+[0 1]),py(segs(row,2)+[0 1]),'.','MarkerSize',20,'color',get(h_plot,'color'));
    fprintf('Intersection %d\n',row);
    pause
    for idx_offset = [-5:5]
      fprintf('  Point offset %d\n', idx_offset);
      idxs = segs(row,1)+idx_offset;
      idxs = idxs(idxs>=1 & idxs<=length(px));
      h_idx_plot = plot(px(idxs),py(idxs),'.','MarkerSize',40);
      pause;
      delete(h_idx_plot);
    end
  end
  keyboard
end

%% Inside/Border determination
segs_sort = sort(segs(:));

% Start at the first segment (idx = 1)
idx = 1;
% 0: border/keep, 1: inside/remove
state = inpolygon(px(idx)+buffer,py(idx),px,py) ...
  & inpolygon(px(idx)-buffer,py(idx),px,py) ...
  & inpolygon(px(idx),py(idx)+buffer,px,py) ...
  & inpolygon(px(idx),py(idx)-buffer,px,py);
if debug_state_check
  h_state_check = figure();
  plot(px,py)
  hold on;
end
mask = logical(zeros(size(px)));
for seg_idx = 1:length(segs_sort)
  if state == 1
    % Transitioning from inside to border
    seg = [idx segs_sort(seg_idx)];
    idx = segs_sort(seg_idx)+1;
    mask(seg(1):seg(2)) = true;
    if debug_state_check
      plot(px(seg(1):seg(2)),py(seg(1):seg(2)),'r.');
    end
    state = 0;
  else
    % Transitioning from border to inside
    seg = [idx segs_sort(seg_idx)];
    idx = segs_sort(seg_idx)+1;
    if debug_state_check
      plot(px(seg(1):seg(2)),py(seg(1):seg(2)),'g.');
    end
    state = 1;
  end
  
  % Debug code for state check
  if debug_state_check
    state_check = inpolygon(px(idx)+buffer,py(idx),px,py) ...
      & inpolygon(px(idx)-buffer,py(idx),px,py) ...
      & inpolygon(px(idx),py(idx)+buffer,px,py) ...
      & inpolygon(px(idx),py(idx)-buffer,px,py);
    if state ~= state_check
      h_state_check1 = plot(px(idx),py(idx),'x');
      h_state_check2 = plot(px(idx-1),py(idx-1),'<');
      if state
        state_str = 'Inside/Remove';
        state_check_str = 'Border/Keep';
      else
        state_str = 'Border/Keep';
        state_check_str = 'Inside/Remove';
      end
      title(sprintf('State Expected: %s, State Check: %s',state_str,state_check_str));
      warning('Check: "State Expected" should be "Inside/Remove" if the point "x" should be excluded). If not run "state = ~state;"');
      keyboard
      try; delete(h_state_check1); end;
      try; delete(h_state_check2); end;
    end
  end
end
if state
  seg = [idx length(px)];
  mask(seg(1):seg(2)) = true;
  if debug_state_check
    plot(px(seg(1):seg(2)),py(seg(1):seg(2)),'r.');
    hold on;
  end
else
  seg = [idx length(px)];
  if debug_state_check
    plot(px(seg(1):seg(2)),py(seg(1):seg(2)),'g.');
    hold on;
  end
end
if debug_state_check
  keyboard
  try; delete(h_state_check); end;
end

%% Construct polygon from outside segments
idxs = {[]};
px_out = {[]};
py_out = {[]};
cur_idx = find(~mask,1);
direction = 1;
skip_idx = NaN;
if debug_plots
  h_plots = figure();
  clf;
  plot(px,py,'b.-');
  hold on;
  h_polygon = plot(NaN,NaN,'k.','markersize',26);
  h_mask = plot(NaN,NaN,'r.','markersize',13);
  h_cur = plot(NaN,NaN,'g-','linewidth',4);
  h_cur(2) = plot(NaN,NaN,'gx','markersize',8,'linewidth',2);
  h_cur(3) = plot(NaN,NaN,'g>','markersize',8,'linewidth',2);
  set(h_mask, 'XData',px(mask), 'YData',py(mask));
  set(h_polygon, 'XData',px_out{end}, 'YData',py_out{end});
end
while any(~mask)
  % Handle the situation where we reach the start/end of the polygon
  if cur_idx < 1 || cur_idx > length(px)
    cur_idx = find(~mask,1);
    if debug_plots
      fprintf('Jump from %d to %d\n', old_cur_idx, cur_idx);
      set(h_mask, 'XData',px(mask), 'YData',py(mask));
      set(h_polygon, 'XData',px_out{end}, 'YData',py_out{end});
      keyboard
    end
    idxs{end+1} = [];
    px_out{end+1} = [];
    py_out{end+1} = [];
  end
  
  if debug_plots
    % Plot current segment
    if cur_idx+direction >= 1 && cur_idx+direction <= length(px)
      set(h_cur(1), 'XData',px(cur_idx+[0 direction]), 'YData',py(cur_idx+[0 direction]));
      set(h_cur(3), 'XData',px(cur_idx+direction), 'YData',py(cur_idx+direction));
    else
      set(h_cur(1), 'XData',NaN, 'YData',NaN);
      set(h_cur(3), 'XData',NaN, 'YData',NaN);
    end
    set(h_cur(2), 'XData',px(cur_idx), 'YData',py(cur_idx));
    set(h_mask, 'XData',px(mask), 'YData',py(mask));
    set(h_polygon, 'XData',px_out{end}, 'YData',py_out{end});
    old_cur_idx = cur_idx;
  end
  
  % Look to see if the current segment crosses a boundary
  match_idxs = find(segs(:) == cur_idx);
  
  % Choose the crossing that is closest and also ensure that the crossing
  % has points that still need to be traversed (mask == false).
  min_dist = inf;
  min_match_idx_idx = [];
  for match_idx_idx = 1:length(match_idxs)
    row = mod(match_idxs(match_idx_idx)-1,size(segs,1)) + 1;
    col = ~floor((match_idxs(match_idx_idx)-1)/size(segs,1)) + 1;
    dist = (px(cur_idx-direction)-x0(row)).^2 + (px(cur_idx-direction)-x0(row)).^2;
    if dist < min_dist && (~mask(segs(row,col)) || ~mask(segs(row,col)+1))
      min_dist = dist;
      min_match_idx_idx = match_idx_idx;
    end
  end
  match_idx = match_idxs(min_match_idx_idx);

  if ~isempty(match_idx)
    % There is at least one crossing to follow, get the indices for the
    % crossing segment the segment. (These are in the opposite col to the
    % current segment.)
    row = mod(match_idx-1,size(segs,1)) + 1;
    col = ~floor((match_idx-1)/size(segs,1)) + 1; % Switch columns

    if ~mask(cur_idx)
      mask(cur_idx) = true;
      idxs{end}(end+1) = cur_idx;
      px_out{end}(end+1) = px(cur_idx);
      py_out{end}(end+1) = py(cur_idx);
    end

    if ~mask(segs(row,col)) && ~mask(segs(row,col)+1)
      % Both directions of crossing are open, don't change direction
      cur_idx = segs(row,col);
      if direction
        cur_idx = segs(row,col)+1;
      else
        cur_idx = segs(row,col);
      end
    elseif ~mask(segs(row,col)+1)
      % Forward on the crossing is open
      direction = +1;
      cur_idx = segs(row,col)+1;
    elseif ~mask(segs(row,col))
      % Backward on the crossing is open
      direction = -1;
      cur_idx = segs(row,col);
    end

    if debug_plots
      fprintf('Jump from %d to %d\n', old_cur_idx, cur_idx);
      keyboard
    end

  end
  
  if mask(cur_idx)
    % We have closed a polygon, look for any other polygons
    cur_idx = find(~mask,1);
    
    if debug_plots
      fprintf('Jump from %d to %d\n', old_cur_idx, cur_idx);
      keyboard
    end
      
    idxs{end+1} = [];
    px_out{end+1} = [];
    py_out{end+1} = [];
    
  else
    % Continue along the current path
    mask(cur_idx) = true;
    idxs{end}(end+1) = cur_idx;
    px_out{end}(end+1) = px(cur_idx);
    py_out{end}(end+1) = py(cur_idx);
    cur_idx = cur_idx + direction;
  end
end
if debug_plots
  set(h_polygon,'XData',NaN,'YData',NaN);
  set(h_mask,'XData',NaN,'YData',NaN);
  set(h_cur,'XData',NaN,'YData',NaN);
  set(h_cur(2),'XData',NaN,'YData',NaN);
  set(h_cur(3),'XData',NaN,'YData',NaN);
  
  for idx = 1:length(px_out)
    h_final(idx) = plot(px_out{idx},py_out{idx},'.-','markersize',10);
    hold on;
  end
  keyboard;
end

% Remove zero area and trim self-intersecting polygons
if debug_poly
  clf;
end
for idx = length(px_out):-1:1

  if polyarea(px_out{idx},py_out{idx}) == 0
    px_out(idx) = [];
    py_out(idx) = [];
    idxs(idx) = [];
  else
    if length(px_out{idx}) >= 4
      px = [px_out{idx} px_out{idx}(1)];
      py = [py_out{idx} py_out{idx}(1)];
      [~,~,segs] = tomo.selfintersect(px,py);
      
      % Remove segments that include overlapping points (i.e. if either end of
      % the segment has the same position as either end of the intersecting
      % segment, remove that segment). HACK: The assumption here is that
      % any overlapping point means that the segment overlaps another segment...
      % a reasonable assumption for single/double precision 3D surfaces.
      overlap_mask = logical(zeros(size(segs,1),1));
      for seg_idx=1:size(segs,1)
        if any(px(segs(seg_idx,1)+[0 1 0 1]) == px(segs(seg_idx,2)+[0 1 1 0]) ...
            & py(segs(seg_idx,1)+[0 1 0 1]) == py(segs(seg_idx,2)+[0 1 1 0]))
          overlap_mask(seg_idx) = true;
        end
      end
      segs = segs(~overlap_mask,:);

      remove_mask = zeros(size(px_out{idx}));
      new_start = max(segs(segs<length(px_out{idx})/2));
      if isempty(new_start)
        new_start = 1;
      end
      new_end = min(segs(segs>=length(px_out{idx})/2));
      if isempty(new_end)
        new_end = length(px_out{idx});
      end
      px_out{idx} = px_out{idx}(new_start:new_end);
      py_out{idx} = py_out{idx}(new_start:new_end);
      idxs{idx} = idxs{idx}(new_start:new_end);
    end
    if polyarea(px_out{idx},py_out{idx}) == 0
      px_out(idx) = [];
      py_out(idx) = [];
      idxs(idx) = [];
    else
      if debug_poly
        h_plot = plot(px_out{idx},py_out{idx},'.-','linewidth',2,'markersize',20);
        hold on;
        keyboard
        delete(h_plot);
        plot(px_out{idx},py_out{idx},'.-','linewidth',2,'markersize',20);
      end
    end
  end
end

if debug_poly
  keyboard
end
if single_polygon_flag
  while length(px_out) > 1
    % Sort polygons based on distances from each other
    min_dists = inf*ones(numel(px_out),numel(px_out));
    min_p1_idxs = zeros(numel(px_out),numel(px_out));
    min_p2_idxs = zeros(numel(px_out),numel(px_out));
    for p1_idx = 1:numel(px_out)
      for p2_idx = p1_idx+1:numel(px_out)
        [min_dists(p1_idx,p2_idx),min_p1_idxs(p1_idx,p2_idx),min_p2_idxs(p1_idx,p2_idx)] ...
          = polygon_distance(px_out{p1_idx},py_out{p1_idx},px_out{p2_idx},py_out{p2_idx});
      end
    end
    
    % Join the two closest polynomials
    dist_idxs = find(isfinite(min_dists(:)));
    [~,dist_min_idx] = min(min_dists(dist_idxs));
    dist_min_idx = dist_idxs(dist_min_idx);
    
    p1_idx = mod(dist_min_idx-1,size(min_dists,1))+1;
    p2_idx = floor((dist_min_idx-1)/size(min_dists,1))+1;
    
    min_p1_idx = min_p1_idxs(dist_min_idx);
    min_p2_idx = min_p2_idxs(dist_min_idx);
    
    % Is p2 inside p1?
    in_poly = inpolygon(px_out{p2_idx}(1),py_out{p2_idx}(1),px_out{p1_idx},py_out{p1_idx});
    
    % If p2 is inside of p1, then it needs to be opposite cw/ccw direction as p1
    % If p2 is outside of p1, then is needs to have the same cw/ccw direction as p1
    cw1 = ispolycw(px_out{p1_idx},py_out{p1_idx});
    cw2 = ispolycw(px_out{p2_idx},py_out{p2_idx});
    if ~in_poly && cw1 || in_poly && ~cw1
      if ~cw2
        [px_out{p2_idx},py_out{p2_idx}] = poly2cw(px_out{p2_idx},py_out{p2_idx});
        min_p2_idx = numel(px_out{p2_idx}) - min_p2_idx + 1;
        idxs{p2_idx} = idxs{p2_idx}(end:-1:1);
      end
    else
      if cw2
        [px_out{p2_idx},py_out{p2_idx}] = poly2ccw(px_out{p2_idx},py_out{p2_idx});
        min_p2_idx = numel(px_out{p2_idx}) - min_p2_idx + 1;
        idxs{p2_idx} = idxs{p2_idx}(end:-1:1);
      end
    end
    
    % Join polygon 2 to polygon 1 at the closest point
    px_out{p1_idx} = [px_out{p1_idx}(1:min_p1_idx), px_out{p2_idx}([min_p2_idx:end,1:min_p2_idx]), px_out{p1_idx}(min_p1_idx:end)];
    py_out{p1_idx} = [py_out{p1_idx}(1:min_p1_idx), py_out{p2_idx}([min_p2_idx:end,1:min_p2_idx]), py_out{p1_idx}(min_p1_idx:end)];
    idxs{p1_idx} = [idxs{p1_idx}(1:min_p1_idx), idxs{p2_idx}([min_p2_idx:end,1:min_p2_idx]), idxs{p1_idx}(min_p1_idx:end)];
    
    % Remove polygon 2 from the list
    px_out(p2_idx) = [];
    py_out(p2_idx) = [];
    idxs(p2_idx) = [];
  end
  
  px_out = px_out{1};
  py_out = py_out{1};
  idxs = idxs{1};
  
else
  [px_out,py_out] = polyjoin(px_out,py_out);
  [idxs] = polyjoin(idxs,idxs);
end
if debug_plots
  try; delete(h_final); end;
  plot(px_out,py_out,'m.-','markersize',10);
  keyboard;
  try; delete(h_plots); end;
end

