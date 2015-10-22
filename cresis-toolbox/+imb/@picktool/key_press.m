function key_press(obj,src,event)

if ~isempty(event.Key)
  switch event.Key
    
    case 'f1'
      % print out help for this tool
      fprintf('--------------------%s Tool Help--------------------\n\n',obj.tool_name);
      fprintf('%s',obj.help_string);
    case 'p'
      obj.close_win();      
  end
    
end

return