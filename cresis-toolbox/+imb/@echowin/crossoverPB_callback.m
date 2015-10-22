function crossoverPB_callback(obj,hObj,event)
% crossoverPB_callback(obj,hObj,event)
%
% Shows the cross over window

if ~obj.crossovers_en
  obj.crossovers_en = true;
  obj.load_crossovers();
  obj.eg.crossovers.gui.visibility_toggle(true);
end

obj.eg.crossovers.gui.figure_visibility_toggle(true);

end
