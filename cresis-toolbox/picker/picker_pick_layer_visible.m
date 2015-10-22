function picker_pick_layer_visible
% picker_pick_layer_visible
%
% Function supporting pick figure. Called from picker_pick
% and picker_pick_key.
%

global hui; % hui: user interface handles
global gCtrl; % Geobase: contains all the geographic info

% =================================================================
% Notes
% =================================================================

% gCtrl.tool.layer_visible_pick
%  element purpose
%  1       1: show manual/auto layers, 2: show quality by color
%     Shift-Space Bar toggles (which also forces master visible layer to on)
%  2       boolean flag (master layers visible on/off)
%     Space Bar toggles
%  3       boolean flag (layer 1 visible on/off)
%     Shift-1 toggles
%  4       boolean flag (layer 2 visible on/off)
%     Shift-2 toggles
%  5       boolean flag (layer 3 visible on/off)
%     Shift-3 toggles

% Layer 1 Manual: hui.pickfig.layer_h(1)
% Layer 1 Automated: hui.pickfig.layer_h(2)
% Layer 2 Manual: hui.pickfig.layer_h(3)
% Layer 2 Automated: hui.pickfig.layer_h(4)

% Layer 1 (good,moderate,derived): hui.pickfig.quality_h(1:3)
% Layer 2 (good,moderate,derived): hui.pickfig.quality_h(1:3)


if gCtrl.tool.layer_visible_pick(2) == 1 && gCtrl.tool.layer_visible_pick(3) == 1
  % =================================================================
  % Layer 1 Is Visible
  % =================================================================
  if gCtrl.tool.layer_visible_pick(1) == 1
    for idx = 1:2
      set(hui.pickfig.layer_h(idx),'Visible','on');
    end
    for idx = 1:3
      set(hui.pickfig.quality_h(idx),'Visible','off');
    end
  else
    for idx = 1:2
      set(hui.pickfig.layer_h(idx),'Visible','off');
    end
    for idx = 1:3
      set(hui.pickfig.quality_h(idx),'Visible','on');
    end
  end
else
  % =================================================================
  % Layer 1 Is Not Visible
  % =================================================================
  for idx = 1:2
    set(hui.pickfig.layer_h(idx),'Visible','off');
  end
  for idx = 1:3
    set(hui.pickfig.quality_h(idx),'Visible','off');
  end
end

if gCtrl.tool.layer_visible_pick(2) == 1 && gCtrl.tool.layer_visible_pick(4) == 1
  % =================================================================
  % Layer 2 Is Visible
  % =================================================================
  if gCtrl.tool.layer_visible_pick(1) == 1
    for idx = 3:4
      set(hui.pickfig.layer_h(idx),'Visible','on');
    end
    for idx = 4:6
      set(hui.pickfig.quality_h(idx),'Visible','off');
    end
  else
    for idx = 3:4
      set(hui.pickfig.layer_h(idx),'Visible','off');
    end
    for idx = 4:6
      set(hui.pickfig.quality_h(idx),'Visible','on');
    end
  end
else
  % =================================================================
  % Layer 2 Is Not Visible
  % =================================================================
  for idx = 3:4
    set(hui.pickfig.layer_h(idx),'Visible','off');
  end
  for idx = 4:6
    set(hui.pickfig.quality_h(idx),'Visible','off');
  end
end

% if gCtrl.tool.layer_visible_pick(2) == 1 && gCtrl.tool.layer_visible_pick(5) == 1
%   % =================================================================
%   % Layer 3 Is Visible
%   % =================================================================
%   if gCtrl.tool.layer_visible_pick(1) == 1
%     for idx = 5:6
%       set(hui.pickfig.layer_h(idx),'Visible','on');
%     end
%   else
%     for idx = 5:6
%       set(hui.pickfig.layer_h(idx),'Visible','off');
%     end
%   end
% else
%   % =================================================================
%   % Layer 3 Is Not Visible
%   % =================================================================
%   for idx = 5:6
%     set(hui.pickfig.layer_h(idx),'Visible','off');
%   end
% end

return;
