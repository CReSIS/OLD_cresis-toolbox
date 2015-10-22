function update_frame_and_sourceLB(obj,frm)
% echowin.update_frame_and_sourceLB(obj,frm)
%
% Updates the frameLB with selection frm and updates the sourceLB for the
% sources that are available for this frame.

% Set the currently selected frame in frameLB
set(obj.left_panel.frameLB,'Value',frm);

% Set the currently selected frame in the figure's title bar
if strcmpi(class(obj.h_fig),'double')
  set(obj.h_fig,'name',sprintf('%d: %s', obj.h_fig, obj.eg.frame_names{frm}));
else
  set(obj.h_fig,'name',sprintf('%d: %s', obj.h_fig.Number, obj.eg.frame_names{frm}));
end

% Find the sources that are available for the currently selected frame
source_mask = any(obj.eg.source_fns_existence(frm,:,:),3);

% Get the current source so that we can try to keep the string selected after
% updating the source listbox
source_idx = get(obj.left_panel.sourceLB,'Value');
current_sources = get(obj.left_panel.sourceLB,'String');
source_idx = find(strcmp(current_sources{source_idx},obj.eg.sources(source_mask)));

% Keep matlab from complaining when current Value ends up being larger than
% the length of sources.
set(obj.left_panel.sourceLB,'Value',1);

% Update sourceLB with currently available list of sources
set(obj.left_panel.sourceLB,'String',obj.eg.sources(source_mask));

if ~isempty(source_idx)
  % Update with old source string
  set(obj.left_panel.sourceLB,'Value',source_idx);
end

end
