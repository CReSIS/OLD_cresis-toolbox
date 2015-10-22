function [comp_time radar_time] = basic_load_mcrds_time(filename,hdr)
% [comp_time radar_time] = basic_load_mcrds_time(filename,hdr)
%
% filename = string containing MCRDS raw file
% hdr = header returned from basic_load_mcrds_hdr
%
% comp_time = double vector containing computer time for each record
%   in the MCRDS file (ANSI-C standard, seconds since Jan 1, 1970
% radar_time = double vector containing 10 MHz clock counts for each
%   record in the MCRDS file (radar_time corresponds to the transmit
%   event for the first pulse in a waveform/presum group???)
%
% Author: John Paden
% 
% See also: basic_load_mcrds_hdr.m, basic_load_mcrds_time.m,
%   basic_load_mcords.m, basic_load_mcords2.m, basic_load_fmcw.m,
%   basic_load_accum.m

block_length = hdr.IndexRecordStop(end,end);

% Open file
fid = fopen(filename ,'r');
% Seek to first record where comp_time/radar_time stored
status = fseek(fid,hdr.IndexData+4,-1);
% Read in all comp_time/radar_times
raw_input = fread(fid,[4 hdr.NumberRecords],'4*uint32=>double',8+2*block_length);
% Close file
fclose(fid);
% Parse raw inputs
comp_time = raw_input(1,:) + raw_input(2,:)*1e-6;
radar_time = (raw_input(4,:)*2^32 + raw_input(3,:))*1e-7;

return;

