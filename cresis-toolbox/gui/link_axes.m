function link_axes(fig_list,option)
% link_axes(fig_list,option)
%
% Simple function to link axes more quickly.
%
% fig_list: list of figure handles or numbers, default is all figures
% option: options to linkaxes, default is 'xy'
%
% Author: John Paden
%
% See also: 

h_axes = [];

if ~exist('fig_list','var') || isempty(fig_list)
  fig_list = get(0,'Children');
end
if ~exist('option','var') || isempty(option)
  option = 'xy';
end

for f_idx = 1:length(fig_list)
  h_children = get(fig_list(f_idx),'Children');
  
  for h_idx = 1:length(h_children)
    if isa(h_children(h_idx),'matlab.graphics.axis.Axes')
      h_axes(end+1) = h_children(h_idx);
    end
  end
end

linkaxes(h_axes,option);

end
