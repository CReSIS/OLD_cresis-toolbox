function picker_update_xaxis(hObj,event)
% This function updates XTickLabel for along-track and rangeline modes when
% matlab zoom tool, picker zoom tool, and left/right with and without
% modifier.

global gCtrl;
global hui;

gCtrl.source.along_track = geodetic_to_along_track(gCtrl.source.geo{...
        gCtrl.source.cur_pick}.Latitude,gCtrl.source.geo{...
        gCtrl.source.cur_pick}.Longitude,gCtrl.source.geo{...
        gCtrl.source.cur_pick}.Elevation);
      
% Find current figure name
cur_fig = get(gcf,'Name');

if strcmp(cur_fig,'3: Pick') % for pickfig

  if gCtrl.source.rlines_flag(1) == 0 % in along-track mode
    % Update XTick to contain accurate range line indexes
    set(hui.pickfig.axes.handle,'XTickLabel',get(hui.pickfig.axes.handle,'XTick'));

    % Take updated XTick and find corresponding along track indexes
    cur_xticklabel = get(hui.pickfig.axes.handle,'XTickLabel');      
    new_xticklabel = zeros(size(cur_xticklabel,1),1);
    for label_idx = 1:length(cur_xticklabel)
      new_xticklabel(label_idx) = roundn(gCtrl.source.along_track(round(str2double(cur_xticklabel(label_idx,:))))*1e-3,-1);
    end
    set(hui.pickfig.axes.handle,'XTickLabel',new_xticklabel);
    set(get(gca,'XLabel'),'String','Along Track Distance (km)')
  elseif gCtrl.source.rlines_flag(1) == 1 % In range line mode
    % Reset ranglines
    set(hui.pickfig.axes.handle,'XTickLabel',get(hui.pickfig.axes.handle,'XTick'));
    set(get(gca,'XLabel'),'String','Range line')
  end
  
elseif strcmp(cur_fig,'4: View') % for viewfig
  
  if gCtrl.source.rlines_flag(2) == 0 % in along-track mode
    % Update XTick to contain accurate range line indexes
    set(hui.viewfig.axes.handle,'XTickLabel',get(hui.viewfig.axes.handle,'XTick'));

    % Take updated XTick and find corresponding along track indexes
    cur_xticklabel = get(hui.viewfig.axes.handle,'XTickLabel');      
    new_xticklabel = zeros(size(cur_xticklabel,1),1);
    for label_idx = 1:length(cur_xticklabel)
      new_xticklabel(label_idx) = roundn(gCtrl.source.along_track(round(str2double(cur_xticklabel(label_idx,:))))*1e-3,-1);
    end
    set(hui.viewfig.axes.handle,'XTickLabel',new_xticklabel);
    set(get(gca,'XLabel'),'String','Along Track Distance (km)')
  elseif gCtrl.source.rlines_flag(2) == 1 % In range line mode
    % Reset ranglines
    set(hui.viewfig.axes.handle,'XTickLabel',get(hui.viewfig.axes.handle,'XTick'));
    set(get(gca,'XLabel'),'String','Range line')
  end
  
end

return