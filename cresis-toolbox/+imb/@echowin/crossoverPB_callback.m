function crossoverPB_callback(obj,hObj,event)
% crossoverPB_callback(obj,hObj,event)
%
% Shows the cross over window

if ~obj.crossovers.en
  obj.crossovers.en = true;
  obj.load_crossovers();
  obj.crossovers.gui.visibility_toggle(true);
end

obj.crossovers.gui.figure_visibility_toggle(true);

end
