function copy_fig(fh,fig_num)
% function copy_fig(fh,fig_num)
%
% Copy all properties of a figure window into a new figure window.  
% 
% Inputs:
%   fh: figure handle of figure to be copied
%   fig_num (optional): figure number to be copied to
%
% Examples:
%   1) copy_fig(gcf)
%   2) fig1 = figure(1);
%      copy_fig(fig1,2)
%
% Author: Logan Smith

error(nargchk(1,2,nargin,'struct'));
h1=fh;
if nargin > 2
    h2=figure(fig_num); clf;
else
    h2=figure;
end
objects=allchild(h1);
copyobj(get(h1,'children'),h2);

end