function update_caxis(obj,varargin)

if ~ishandle(obj.img)
  return;
end

cmap = get(obj.h_gui.colormapPM,'Value');
if cmap == 1
  base_cmap = 1-gray(256);
else
  base_cmap = jet(256);
end
hist_exp = obj.h_gui.hist_slider.get_value();
cmap = interp1(linspace(0,1,256).', base_cmap, linspace(0,1,256).^(10.^(hist_exp/10)));

if ~get(obj.h_gui.caxis_autoCB,'Value')
  clims = obj.h_gui.min_slider.get_value();
  clims(2) = obj.h_gui.max_slider.get_value();
  clims = sort(clims);
  
  caxis(get(obj.img,'parent'),clims);
  colormap(get(obj.img,'parent'),cmap);
  
else
  xlims = xlim(get(obj.img,'parent'));
  ylims = ylim(get(obj.img,'parent'));
  X = get(obj.img,'XData');
  Y = get(obj.img,'YData');
  C = get(obj.img,'CData');
  
  if isempty(X)
    return;
  end
  xbins = sort(interp1(X([1 end]), [1 size(C,2)], xlims));
  xbins(1) = floor(xbins(1));
  xbins(2) = ceil(xbins(2));
  if isnan(xbins(1))
    xbins(1) = 1;
  end
  if isnan(xbins(2))
    xbins(2) = size(C,2);
  end
  
  ybins = sort(interp1(Y([1 end]), [1 size(C,1)], ylims));
  ybins(1) = floor(ybins(1));
  ybins(2) = ceil(ybins(2));
  if isnan(ybins(1))
    ybins(1) = 1;
  end
  if isnan(ybins(2))
    ybins(2) = size(C,1);
  end

  clims = [finitemin(finitemin(C(ybins(1):ybins(end),xbins(1):xbins(end)))) finitemax(finitemax(C(ybins(1):ybins(end),xbins(1):xbins(end))))];
  
  if numel(clims) == 2 && all(isfinite(clims))
    caxis(get(obj.img,'parent'),clims);
    colormap(get(obj.img,'parent'),cmap);
    
    obj.h_gui.min_slider.set_value(clims(1));
    obj.h_gui.max_slider.set_value(clims(end));
  end
end

return

