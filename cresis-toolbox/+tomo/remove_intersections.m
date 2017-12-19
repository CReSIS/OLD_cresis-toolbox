function [px_out,py_out,idxs] = remove_intersections(px,py,buffer)
% [px_out,py_out,idxs] = tomo.remove_intersections(px,py,buffer)
%
% Useful for gridding data that comes from self-intersecting swaths. This
% finds the bounding polygon with no self-intersections. This code is not
% tested very well and may have bugs. There is also an assumption/HACK of
% removing intersecting segments with any matching coordinates.
%
% px,py: polygon (may or may not have self intersections)
% buffer: HACK with inpolygon which allows us to check if a point is
% inside of a polygon and not just on the edge. inpolygon normally returns
% true for any point on the edge of the polygon. We do this by checking
% (roughly) that all points within buffer distance from the point are inside
% the polygon. If they are, we assume that the point is in the polygon and
% not just on the edge of the polygon.
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
%   [px_out,py_out,idxs] = tomo.remove_intersections(px,py,0.1);
%    px_check = zeros(size(idxs));
%    px_check(isnan(idxs)) = NaN;
%    px_check(~isnan(idxs)) = px(idxs(~isnan(idxs)));
%    any(px_out(~isnan(px_out)) ~= px_check(~isnan(px_out)));
%
% Author: John Paden

% Find all the segments that intersect
[~,~,segs] = tomo.selfintersect(px,py);

% Remove segments that include overlapping points (i.e. if either end of
% the segment has the same position as either end of the intersecting
% segment, remove that segment). HACK: The assumption here is that a single
% overlapping point only occurs when a segment overlaps... a reasonable
% assumption for single/double precision 3D surfaces.
overlap_mask = logical(zeros(size(segs,1),1));
for seg_idx=1:size(segs,1)
  if any(px(segs(seg_idx,1)+[0 1 0 1]) == px(segs(seg_idx,2)+[0 1 1 0]) ...
      & py(segs(seg_idx,1)+[0 1 0 1]) == py(segs(seg_idx,2)+[0 1 1 0]))
    overlap_mask(seg_idx) = true;
  end
end
segs = segs(~overlap_mask,:);

if 0
  clf;
  plot(px,py);
  hold on
  for row = 1:size(segs,1)
    % Draw segment
    h_plot = plot(px(segs(row,1)+[0 1]),py(segs(row,1)+[0 1]),'-','LineWidth',2);
    % Draw intersecting segment
    plot(px(segs(row,2)+[0 1]),py(segs(row,2)+[0 1]),'.','MarkerSize',20,'color',get(h_plot,'color'));
    row
    pause
  end
  keyboard
end

segs_sort = sort(segs(:));

idx = 1;
state = inpolygon(px(idx)+buffer,py(idx),px,py) ...
  & inpolygon(px(idx)-buffer,py(idx),px,py) ...
  & inpolygon(px(idx),py(idx)+buffer,px,py) ...
  & inpolygon(px(idx),py(idx)-buffer,px,py);
clf;
mask = logical(zeros(size(px)));
for seg_idx = 1:length(segs_sort)
  if state == 1
    % Transitioning from inside to outside
    seg = [idx segs_sort(seg_idx)];
    idx = segs_sort(seg_idx)+1;
    mask(seg(1):seg(2)) = true;
    %plot(px(seg(1):seg(2)),py(seg(1):seg(2)),'r.');
    %hold on;
    state = 0;
  else
    % Transitioning from outside to inside
    seg = [idx segs_sort(seg_idx)];
    idx = segs_sort(seg_idx)+1;
    %plot(px(seg(1):seg(2)),py(seg(1):seg(2)),'g.');
    %hold on;
    state = 1;
  end
  
  state_check = inpolygon(px(idx)+buffer,py(idx),px,py) ...
    & inpolygon(px(idx)-buffer,py(idx),px,py) ...
    & inpolygon(px(idx),py(idx)+buffer,px,py) ...
    & inpolygon(px(idx),py(idx)-buffer,px,py);
  if state ~= state_check
    figure;
    plot(px,py)
    hold on;
    plot(px(idx),py(idx),'x');
    hold off;
    if state
      state_str = 'Inside';
      state_check_str = 'Outside';
    else
      state_str = 'Outside';
      state_check_str = 'Inside';
    end
    title(sprintf('State Expected: %s, State Check: %s',state_str,state_check_str));
    warning('Identify the correct state (state is true if the point is inside the polygon).');
    keyboard
  end
end
if state
  seg = [idx length(px)];
  mask(seg(1):seg(2)) = true;
  %plot(px(seg(1):seg(2)),py(seg(1):seg(2)),'r.');
  %hold on;
end
%keyboard

%plot(px(~mask),py(~mask),'.')
%keyboard
%clf

idxs = [];
px_out = {[]};
py_out = {[]};
cur_idx = find(~mask,1);
direction = 1;
clf;
while any(~mask)
  if mask(cur_idx)
    match_idx = find(segs(:) == cur_idx-direction);
    row = mod(match_idx-1,size(segs,1)) + 1;
    col = ~floor(match_idx/size(segs,1)) + 1;
    cur_idx = segs(row,col);
    if mask(cur_idx)
      if cur_idx > 1 && ~mask(cur_idx-1)
        direction = -1;
        cur_idx = cur_idx + direction;
      elseif cur_idx < length(px) && ~mask(cur_idx+1)
        direction = 1;
        cur_idx = cur_idx + direction;
      else
        idxs(end+1,1) = NaN;
        px_out{end+1} = [];
        py_out{end+1} = [];
        while mask(cur_idx)
          cur_idx = mod(cur_idx,length(mask)) + direction;
        end
      end
    end
    %plot(px(good_idxs{end}),py(good_idxs{end}));
  else
    last_state = mask(cur_idx);
    mask(cur_idx) = true;
    idxs(end+1,1) = cur_idx;
    px_out{end}(end+1) = px(cur_idx);
    py_out{end}(end+1) = py(cur_idx);
    cur_idx = cur_idx + direction;
  end
end

[px_out,py_out] = polyjoin(px_out,py_out);
%plot(px_out,py_out);
