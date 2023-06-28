function save_mat_callback(obj,h_obj,event)
% imagewin.save_mat_callback(obj,h_obj,event)
%
% Callback when save mat button is pushed

if ~isempty(obj.save_mat_fh)
  obj.save_mat_fh(h_obj,event);
else
  warning('Property save_mat_fh is not set');
end

end
