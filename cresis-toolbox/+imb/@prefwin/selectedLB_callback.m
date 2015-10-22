function seasonPB_callback(obj, hObj, event)

if strcmpi(get(obj.h_fig,'SelectionType'),'open')
  removePB_callback(obj, hObj, event);
end

return
