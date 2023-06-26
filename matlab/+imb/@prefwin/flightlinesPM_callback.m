function flightlinesPM_callback(obj,status,event)

temp = get(obj.h_gui.flightlinesPM,'String');

flightline = temp{get(obj.h_gui.flightlinesPM,'Value')};

if strcmpi(flightline,'Connect to OPS')
  if obj.ops_connect()
    set(obj.h_gui.flightlinesPM,'Value',1);
    return;
  end
end

obj.season_update();
