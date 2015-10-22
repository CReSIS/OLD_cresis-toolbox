function cmds_list = attach_echowin(obj,h_doc)
% cmds_list = attach_echowin(obj, h_doc)
%
% Member function of undo_stack. Attaches the specified document to
% the undo stack list. Returns a list of commands that need to be
% executed with the implied cmds_direction of 'redo'.

obj.doc_list{end+1} = h_doc;

cmds_list = obj.stack(1:obj.pointer);

end
