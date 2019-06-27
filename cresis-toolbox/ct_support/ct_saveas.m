function ct_saveas(fig,fn,varargin)
% ct_saveas(fig,fn,varargin)
%
% Function calls saveas, but adds a creation function callback to force the
% "Visible" property to "on" when the figure is opened. This is useful when
% creating invisible figures and saving them.

% Set the create creation function to force the visibility to on when
% opening.
try
  set(fig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
end
saveas(fig,fn,varargin{:});
