function handle = get_parent_figure(handle)
% handle = get_parent_figure(handle)
%
% Recurses through the parents of the handle passed in and finds the first
% parent that is a figure or returns [] if no parents are of type figure.

while ~isempty(handle) & ~strcmp('figure', get(handle,'type'))
  handle = get(handle,'parent');
end

end
