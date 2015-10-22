function frm = repair_frame(yyyymmdd,segment,frame)
%
% Repair FRAME from YYYYMMDD,SEGMENT,FRAME
%
% [frm] = repair_frame(yyyymmdd,segment,frame)
%
%   yyyymmdd = YYYYMMDD variable (1:8) of FRAME
%   segment = SEG variable (9:10) of FRAME
%   frame = FRM variable (11:end) of FRAME
%
% Example:
%
%   yyyymmdd = '20110511';
%   segment = '01';
%   frame = '001';
%   
%   [FRAME] = repair_frame(yyyymmdd,segment,frame);
%
%   FRAME = '2011051101001';
%
% Author: Kyle Purdon

% Pre-Allocate New Variable
frm = cell(length(frame),1);

% Fill New Variable  
for idx = 1:length(frame)
  frm{idx} = strcat(yyyymmdd{idx},segment{idx},frame{idx});
end
end