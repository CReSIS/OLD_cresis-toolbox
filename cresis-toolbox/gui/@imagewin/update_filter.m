function update_filter(obj,varargin)

if ~ishandle(obj.img)
  return;
end

param1 = obj.h_gui.param1_slider.get_value();
param2 = obj.h_gui.param2_slider.get_value();

sig_filt = inline(obj.h_gui.signal_listbox.get_string());
noise_filt = inline(obj.h_gui.noise_listbox.get_string());

if nargin(sig_filt) == 1 && nargin(noise_filt) == 1
  set(obj.img,'CData',real(sig_filt(obj.CData) - noise_filt(obj.CData)));
elseif nargin(sig_filt) == 2 && nargin(noise_filt) == 1
  set(obj.img,'CData',real(sig_filt(obj.CData,param1) - noise_filt(obj.CData)));
elseif nargin(sig_filt) == 3 && nargin(noise_filt) == 1
  set(obj.img,'CData',real(sig_filt(obj.CData,param1,param2) - noise_filt(obj.CData)));
elseif nargin(sig_filt) == 1 && nargin(noise_filt) == 2
  set(obj.img,'CData',real(sig_filt(obj.CData) - noise_filt(obj.CData,param1)));
elseif nargin(sig_filt) == 2 && nargin(noise_filt) == 2
  set(obj.img,'CData',real(sig_filt(obj.CData,param1) - noise_filt(obj.CData,param1)));
elseif nargin(sig_filt) == 3 && nargin(noise_filt) == 2
  set(obj.img,'CData',real(sig_filt(obj.CData,param1,param2) - noise_filt(obj.CData,param1)));
elseif nargin(sig_filt) == 1 && nargin(noise_filt) == 3
  set(obj.img,'CData',real(sig_filt(obj.CData) - noise_filt(obj.CData,param1,param2)));
elseif nargin(sig_filt) == 2 && nargin(noise_filt) == 3
  set(obj.img,'CData',real(sig_filt(obj.CData,param1) - noise_filt(obj.CData,param1,param2)));
elseif nargin(sig_filt) == 3 && nargin(noise_filt) == 3
  set(obj.img,'CData',real(sig_filt(obj.CData,param1,param2) - noise_filt(obj.CData,param1,param2)));
end

obj.set_limits(false);
obj.update_caxis();

end
