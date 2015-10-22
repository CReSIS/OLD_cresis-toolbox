function plot_crossovers(obj)
% echowin.plot_crossovers(obj)
%
% Update crossover plot handles. Assumes that obj.eg.image* fields are up
% to date.

physical_constants;

if isempty(obj.eg.crossovers.gps_time)
  obj.eg.crossovers.x_curUnit = [];
else
  obj.eg.crossovers.x_curUnit = interp1(obj.eg.image_gps_time,obj.eg.image_xaxis,obj.eg.crossovers.gps_time,'nearest','extrap');
end

% ======================================================================
%% Convert the data along the y-axis according to the units
% perform y-axis conversion (from twtt)
yaxis_choice = get(obj.left_panel.yaxisPM,'Value');
if yaxis_choice == 1 % TWTT
  % convert crossover values too
  % calculate elevation offsets for crossovers first
  % raw elevation_1/elevation_2 in meters, elev_offset in twtt
  elev_offset = 2*(obj.eg.crossovers.cross_elev-obj.eg.crossovers.source_elev)./c.*1e6;
  % calculate crossover twtts with elevation differences compensated
  obj.eg.crossovers.y_curUnit = obj.eg.crossovers.twtt*1e6 - elev_offset;
  
elseif yaxis_choice == 2 % WGS_84 Elevation
  % convert crossovers too
  elevation = obj.eg.crossovers.source_elev;
  elev_offset = 2*(obj.eg.crossovers.cross_elev-obj.eg.crossovers.source_elev)./c;
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surface,...
    obj.eg.crossovers.gps_time,'linear');
  % correct surface using elevation offset info (since we are using surface
  % for frame1 and not frame2 and the twtt of the crossover is from frame2)
  range = min(obj.eg.crossovers.twtt-elev_offset,surface).*c/2 ...
        + max(0,obj.eg.crossovers.twtt-elev_offset-surface).*c/(sqrt(er_ice)*2);
  obj.eg.crossovers.y_curUnit = elevation - range;
  
elseif yaxis_choice == 3 % Depth/Range
  % convert crossovers too
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surface,...
    obj.eg.crossovers.gps_time,'linear');
  elev_offset = 2*(obj.eg.crossovers.cross_elev-obj.eg.crossovers.source_elev)./c;
  % correct surface using elevation offset info (since we are using surface
  % for frame1 and not frame2 and the twtt of the crossover is from frame2)
  % doesn't work right now for some reason
  obj.eg.crossovers.y_curUnit = min(obj.eg.crossovers.twtt-elev_offset,surface).*c/2 ...
        + max(0,obj.eg.crossovers.twtt-elev_offset-surface).*c/(sqrt(er_ice)*2);
      
elseif yaxis_choice == 4 % Range bin
  % convert crossovers too
  elev_offset = 2*(obj.eg.crossovers.cross_elev-obj.eg.crossovers.source_elev)./c;
  % calculate crossover twtts with elevation differences compensated
  corrected_twtt = obj.eg.crossovers.twtt - elev_offset;
  obj.eg.crossovers.y_curUnit = interp1(obj.eg.time,1:length(obj.eg.time),corrected_twtt);
end

%% Plot crossover data
delete(obj.eg.crossovers.h);
obj.eg.crossovers.h = [];
for idx = 1:length(obj.eg.crossovers.x_curUnit)
  obj.eg.crossovers.h(2*(idx-1)+1) = plot(obj.right_panel.axes.handle, ...
    obj.eg.crossovers.x_curUnit(idx),obj.eg.crossovers.y_curUnit(idx),'o', ...
    'markersize',25);
  obj.eg.crossovers.h(2*(idx-1)+2) = plot(obj.right_panel.axes.handle, ...
    obj.eg.crossovers.x_curUnit(idx),obj.eg.crossovers.y_curUnit(idx),'x', ...
    'markersize',25);
end

end

