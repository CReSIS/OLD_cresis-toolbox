function ct_saveas(fig,fn,varargin)

saveas(fig,fn,varargin{:});

try
  load(fn,'-mat','hgS_070000');
  hgS_070000.properties.Visible = 'on';
  save(fn,'-append','hgS_070000');
end
