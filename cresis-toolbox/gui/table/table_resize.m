function table_resize(hObject, eventdata, handles)
% table_resize(hObject, eventdata, handles)

[hObject,figure] = gcbo;

table = get(hObject,'UserData');

table_draw(table);

return;
