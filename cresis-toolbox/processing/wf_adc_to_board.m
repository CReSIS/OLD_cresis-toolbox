function [board,board_idx,profile] = wf_adc_to_board(param,img)
% [board,board_idx,profile] = wf_adc_to_board(param,img)
%
% Support function for determining which file (referred to as board since
% each board usually generates its own file and may have multiple adcs) a
% particular wf-adc pair comes from. This also returns the profile for more
% complex data mappings like what the Arena systems use.
%
% Inputs:
% param: parameter structure usually from parameter spreadsheet
%   .records
%     .file.version
%     .data_map
% img: list of wf-adc pairs (N by 2 matrix where first column is the wf and
%   the second column is the adc. wf is one indexed. adc is one indexed.
%
% Outputs:
% board: usually one or zero indexed (refers to the naming convention in
%   the raw data filenames). Returns a list of unique boards required for
%   the wf-adc pairs passed in.
% board_idx: is one indexed. Returns a list of the unique board_idxs
%   required for the wf-adc pairs passed in.
% profile: N by M matrix where N is the number of separate profiles 
%   requested (i.e. the number of wf-adc pairs that are passed in) and M
%   is the number of fields required for a particular file version. 
%
% Author: John Paden

if any(param.records.file.version == [402 403])
  % NI systems with 4 adcs per board
  board = unique(floor((img(:,2).'-1)/4));
  board_idx = board+1;
  profile = [];
  
elseif any(param.records.file.version == [410])
  % MCRDS with all adcs on one board
  board = 1;
  board_idx = 1;
  profile = [];
  
elseif any(param.records.file.version == [9 10 103 412])
  % RSS which uses param.records.data_map
  board = zeros(1,size(img,1));
  profile = zeros(size(img,1),2);
  for wf_adc = 1:size(img,1)
    found = false;
    for board_idx = 1:length(param.records.data_map)
      for profile_idx = 1:size(param.records.data_map{board_idx},1)
        wf = param.records.data_map{board_idx}(profile_idx,3);
        adc = param.records.data_map{board_idx}(profile_idx,4);
        if img(wf_adc,1) == wf && img(wf_adc,2) == adc
          board(wf_adc) = board_idx;
          profile(wf_adc,1) = param.records.data_map{board_idx}(profile_idx,1); % mode
          profile(wf_adc,2) = param.records.data_map{board_idx}(profile_idx,2); % subchannel
          found = true;
        end
      end
    end
    if ~found
      error('Did not find wf-adc pair (%d,%d).', img(wf_adc,1), img(wf_adc,2));
    end
  end
  board = unique(board);
  board_idx = board;
  
else
  % All other systems
  board = unique(img(:,2).');
  board_idx = board;
  profile = [];
end

return;
