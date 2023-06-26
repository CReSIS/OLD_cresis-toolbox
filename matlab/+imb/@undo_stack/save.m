function save(obj)
% save(obj)
%
% Notify all documents that a save has happened.

notify(obj,'save_event');

end
