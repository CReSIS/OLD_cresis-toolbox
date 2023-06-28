function plot_crossovers(obj)
% echowin.plot_crossovers(obj)
%
% Update crossover plot handles. Assumes that obj.eg.image* fields are up
% to date.

physical_constants;
vel_air = c/2;
slowness_air = 1/vel_air;
vel_ice = c/(2*sqrt(er_ice));
slowness_ice = 1/vel_ice;

if isempty(obj.crossovers.gps_time)
  obj.crossovers.x_curUnit = [];
else
  obj.crossovers.x_curUnit = interp1(obj.eg.image_gps_time,obj.eg.image_xaxis,obj.crossovers.gps_time,'nearest','extrap');
end

% ======================================================================
%% Convert the data along the y-axis according to the units
% perform y-axis conversion (from twtt)
yaxis_choice = get(obj.left_panel.yaxisPM,'Value');
if yaxis_choice == 1 % TWTT
  % convert crossover values too
  % calculate elevation offsets for crossovers first
  % raw elevation_1/elevation_2 in meters, elev_offset in twtt
  elev_offset = (obj.crossovers.cross_elev-obj.crossovers.source_elev)*slowness_air*1e6;
  % calculate crossover twtts with elevation differences compensated
  obj.crossovers.y_curUnit = obj.crossovers.twtt*1e6 - elev_offset;
  
elseif yaxis_choice == 2 % WGS_84 Elevation
  % convert crossovers too
  elevation = obj.crossovers.source_elev;
  elev_offset = (obj.crossovers.cross_elev-obj.crossovers.source_elev)*slowness_air;
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surf_twtt,...
    obj.crossovers.gps_time,'linear');
  % correct surface using elevation offset info (since we are using surface
  % for frame1 and not frame2 and the twtt of the crossover is from frame2)
  range = min(obj.crossovers.twtt-elev_offset,surface).*vel_air ...
        + max(0,obj.crossovers.twtt-elev_offset-surface).*vel_ice;
  obj.crossovers.y_curUnit = elevation - range;
  
elseif yaxis_choice == 3 % Depth/Range
  % convert crossovers too
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surf_twtt,...
    obj.crossovers.gps_time,'linear');
  elev_offset = (obj.crossovers.cross_elev-obj.crossovers.source_elev)*slowness_air;
  % correct surface using elevation offset info (since we are using surface
  % for frame1 and not frame2 and the twtt of the crossover is from frame2)
  % doesn't work right now for some reason
  obj.crossovers.y_curUnit = min(obj.crossovers.twtt-elev_offset,surface).*vel_air ...
        + max(0,obj.crossovers.twtt-elev_offset-surface).*vel_ice;
      
elseif yaxis_choice == 4 % Range bin
  % convert crossovers too
  elev_offset = (obj.crossovers.cross_elev-obj.crossovers.source_elev)*slowness_air;
  % calculate crossover twtts with elevation differences compensated
  corrected_twtt = obj.crossovers.twtt - elev_offset;
  obj.crossovers.y_curUnit = interp1(obj.eg.time,1:length(obj.eg.time),corrected_twtt);
  
elseif yaxis_choice == 5 % Surface flat
  % convert crossovers too
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surf_twtt,...
    obj.crossovers.gps_time,'linear');
  elev_offset = (obj.crossovers.cross_elev-obj.crossovers.source_elev)*slowness_air;
  % correct surface using elevation offset info (since we are using surface
  % for frame1 and not frame2 and the twtt of the crossover is from frame2)
  % doesn't work right now for some reason
  obj.crossovers.y_curUnit = min(0,obj.crossovers.twtt-elev_offset-surface).*vel_air ...
        + max(0,obj.crossovers.twtt-elev_offset-surface).*vel_ice;
end

%% Plot crossover data
delete(obj.crossovers.h);
obj.crossovers.h = [];
for idx = 1:length(obj.crossovers.x_curUnit)
  obj.crossovers.h(2*(idx-1)+1) = plot(obj.h_axes, ...
    obj.crossovers.x_curUnit(idx),obj.crossovers.y_curUnit(idx),'o', ...
    'markersize',25);
  obj.crossovers.h(2*(idx-1)+2) = plot(obj.h_axes, ...
    obj.crossovers.x_curUnit(idx),obj.crossovers.y_curUnit(idx),'x', ...
    'markersize',25);
end

end

