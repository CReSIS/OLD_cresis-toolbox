function remove_document(obj, h_doc)
% remove_document(obj, h_doc)
%
% Member function of undo_stack. Removes the specified document from
% the undo stack list

for idx = 1:length(obj.doc_list)
  if obj.doc_list{idx} == h_doc;
    obj.doc_list = obj.doc_list([1:idx-1 idx+1:end]);
    return
  end
end

return
