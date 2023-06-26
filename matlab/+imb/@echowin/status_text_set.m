function status_text_set(obj,str,type)
% echowin.status_text_set(obj, dest, str, type)
%
% Prints status messages to the echowin obj's status bar.
% str is the string to be printed
% type specifies properties of the printing.
% Recognized inputs are:
%   'replace': Replace the status display's existing string with str
%   'append': Append str to the status display's existing string
%
% Note that management of the mouseCoordText component of the status bar is
% maintained separately by the button_motion mouse movement callback.
%

switch type
  case 'replace'
    set(obj.right_panel.status_panel.statusText,'String',str);
  case 'append'
    old_str = get(obj.right_panel.status_panel.statusText,'String');
    set(obj.right_panel.status_panel.statusText,'String',strcat(old_str,str));
  case 'default'
    set(obj.right_panel.status_panel.statusText,'String',str);
end
%drawnow;

return;
