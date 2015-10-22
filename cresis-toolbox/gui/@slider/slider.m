classdef slider < handle
  % slider: A class for creating an integrated label, slider, and edit
  
  properties
    h_text
    h_slider
    h_LE
    value_range
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
   
  events
    slider_changed
  end
  
  methods
    function obj = slider(parent,name,value_range,init_value)
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      
      %%% Post Initialization %%%
      % Any code, including access to object
      obj.value_range = value_range;
      
      obj.h_text = uicontrol('Parent',parent);
      set(obj.h_text,'Style','Text');
      set(obj.h_text,'String',name);
      
      obj.h_slider = uicontrol('Parent',parent);
      set(obj.h_slider,'Style','Slider');
      set(obj.h_slider,'Min',min(value_range));
      set(obj.h_slider,'Max',max(value_range));
      set(obj.h_slider,'Value',init_value);
      set(obj.h_slider,'Callback',@obj.slider_callback);
      
      obj.h_LE = uicontrol('Parent',parent);
      set(obj.h_LE,'Style','Edit');
      set(obj.h_LE,'String',sprintf('%g', init_value));
      set(obj.h_LE,'Callback',@obj.slider_callback);
    end
    
    function delete(obj)
    end

    function slider_callback(obj, h_obj, event)
      if h_obj == obj.h_slider
        val = get(obj.h_slider,'Value');
        set(obj.h_LE,'String',sprintf('%g', val));
      elseif h_obj == obj.h_LE
        val = str2double(get(obj.h_LE,'String'));
        if val < min(obj.value_range)
          val = min(obj.value_range);
          set(obj.h_LE,'String',sprintf('%g', val));
        elseif val > max(obj.value_range)
          val = max(obj.value_range);
          set(obj.h_LE,'String',sprintf('%g', val));
        end
        set(obj.h_slider,'Value',val);
      end
      notify(obj,'slider_changed');
    end
    
    function val = get_value(obj)
      val = get(obj.h_slider,'Value');
    end
    
    function set_value(obj, val)
      set(obj.h_slider,'Value',val);
      set(obj.h_LE,'String',sprintf('%g', val));
    end
    
    function val = get_min(obj)
      val = get(obj.h_slider,'Min');
    end
    
    function set_min(obj, val)
      set(obj.h_slider,'Min',val);
    end
    
    function val = get_max(obj)
      val = get(obj.h_slider,'Max');
    end
    
    function set_max(obj, val)
      set(obj.h_slider,'Max',val);
    end
    
  end
  
end


