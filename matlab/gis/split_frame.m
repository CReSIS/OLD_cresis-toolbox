function [yyyymmdd,segment,frame] = split_frame(frm)
% 
% Split FRAME into YYYYMMDD,SEGMENT,FRAME
%
% [yyyymmdd,segment,frame] = split_frame(frm);
% 
% frm = Standard FRAME variable (Cell)
% 
% Example:
%
%   FRAME = '2011051101001';
%
%   [yyyymmdd,segment,frame] = split_frame(FRAME);
%
%    yyyymmdd = '20110511'
%    segment = '01'
%    frame = '001'
%
% Author: Kyle Purdon

% Pre-Allocate New Variables
yyyymmdd = cell(length(frm),1);
segment = cell(length(frm),1);
frame = cell(length(frm),1);

% Fill New Variables
for frm_idx = 1:length(frm)
  yyyymmdd{frm_idx} = frm{frm_idx}(1:8);
  segment{frm_idx} = frm{frm_idx}(9:10);
  frame{frm_idx} = frm{frm_idx}(11:end);
end
end