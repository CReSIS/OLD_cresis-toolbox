function h_fig = get_figures_userdata(h_fig,user_data)
% get_figures_userdata(h_fig,user_data)
%
% Function for clearing or resetting figure's UserData field. See
% get_figures for more details
%
% h_fig: Leave empty to do all figures or specify 1) vector of figure
%   handles, or string containing UserData to match to.
% user_data: string to uniquely identify these new figures (Default is
%   empty)

if nargin<1
  h_fig = [];
end

if nargin<2
  user_data = '';
end

h_figs = get(0,'Children');

if isempty(h_fig)
  mask = true(size(h_figs));
elseif ischar(h_fig)
  mask = false(size(h_figs));
  for fig_idx = 1:length(h_figs)
    if ~isempty(regexpi(get(h_figs(fig_idx),'UserData'), h_fig))
      mask(fig_idx) = true;
    end
  end
else
  mask = false(size(h_figs));
  for fig_idx = 1:length(h_figs)
    if any(h_fig == h_figs(fig_idx).Number)
      mask(fig_idx) = true;
    end
  end
end

fig_idxs = find(mask);
for fig_idx = 1:length(fig_idxs)
  set(h_figs(fig_idxs(fig_idx)),'UserData',sprintf('%s:%d',user_data,fig_idx));
end


