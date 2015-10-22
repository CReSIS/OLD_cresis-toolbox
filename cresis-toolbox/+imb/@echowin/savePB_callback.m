function savePB_callback(obj,hObj,event)
% savePB_callback(obj,hObj,event)

obj.status_text_set(sprintf('(%s) Saving...', datestr(now,'HH:MM:SS')),'replace');
drawnow;

imb.save_undo_stack(obj.undo_stack);

%% Hack to fix GUI object focus
try
  warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
  javaFrame = get(obj.h_fig,'JavaFrame');
  javaFrame.getAxisComponent.requestFocus;
% catch
%   obj.status_text_set(sprintf('Focus error, click inside echogram window before using key shortcuts'),'replace');
end

obj.status_text_set(sprintf('Done %s', datestr(now,'HH:MM:SS')),'append');

end
