function struct_out = empty_struct_element(struct_in)

  field_names = fieldnames(struct_in);
  
  for f = 1:length(field_names)
    
    field_name = field_names{f};
  
    if isstruct(struct_in(1).(field_name))
      struct_out.(field_name) = empty_struct_element(struct_in(1).(field_name));
    else
      struct_out.(field_name) = [];
    end
    
  end

end