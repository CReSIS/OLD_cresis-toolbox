function comp_result = frame_id_comp(frm_id1,frm_id2)
% comp_result = frame_id_comp(frm_id1,frm_id2)
%
% Compares two frame IDs. Specifically built for building
% binary searches of lists of sorted frame IDs.
%
% If frm_id1 is lower than frm_id2, then returns -1
% If frm_id1 is equal to frm_id2, then returns 0
% If frm_id1 is more than frm_id2, then returns 1
%
% Example:
%  frame_id_comp('20110316_01_000','20110406_02_392')
%  frame_id_comp('20110406_02_392','20110406_02_392')
%  frame_id_comp('20110406_02_392','20110316_01_000')
%
% Author: John Paden
%
% See also: picker.m (used by picker's binary search)

comp = frm_id1 - frm_id2;
diff_idx = find(comp~=0,1);
if isempty(diff_idx)
  comp_result = 0;
else
  comp_result = comp(diff_idx);
end

return;
