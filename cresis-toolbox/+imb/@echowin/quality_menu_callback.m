function quality_menu_callback(obj,source,event)

% callback function to give focus back to echowin after quality is changed

try
  warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
  javaFrame = get(obj.h_fig,'JavaFrame');
  javaFrame.getAxisComponent.requestFocus;
catch
  obj.status_text_set(sprintf('Focus error, click inside echogram window before using key shortcuts'),'replace');
end

return;