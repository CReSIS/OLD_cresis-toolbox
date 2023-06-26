function mapsPM_callback(obj,status,event)

temp = get(obj.h_gui.mapsPM,'String');

map = temp{get(obj.h_gui.mapsPM,'Value')};

if strcmpi(map,'Connect to OPS')
  if obj.ops_connect()
    set(obj.h_gui.mapsPM,'Value',1);
    return;
  end
end

obj.season_update();
