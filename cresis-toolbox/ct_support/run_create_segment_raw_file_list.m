% script run_create_segment_raw_file_list
%
% Script for running create_segment_raw_file_list.m
%
% Author: John Paden

% =========================================================================
%% User Settings
% =========================================================================
% base_path = '/cresis/snfs1/data/Accum_Data/2015_Greenland_Ground/20150430/';
% radar_name = 'accum2';

% base_path = '/process/fmcw/snow/'; % USE create_segment_raw_file_list_v2.m
% radar_name = 'snow3';

% base_path = '/process/fmcw/kuband/'; % USE create_segment_raw_file_list_v2.m
% radar_name = 'kuband3';

% base_path = '/process/mcords/';
% radar_name = 'mcords3';

base_path = '/cresis/snfs1/data/SAR/2008_Greenland/080813/dataINSAR2/';
radar_name = 'mcrds';

min_seg_size = 2;

% Breaking of segments based on file name times uses the median time difference
% between file name time stamps and then adds a guard... sometimes the 
% file writes can get behind schedule/have a lot of variance with a single
% segment and larger values should be used. 0.5 is good when this is not
% the case.
% TIME_SEGMENT_GUARD = 0.5;
TIME_SEGMENT_GUARD = 1.5;

% Choose one method for looking at segment breaks:
% allbig: All big files sizes, defined as 90% of the max file size and
%   occuring at least twice, are considered full files and not a segment
%   break. This is the preferred method for FMCW radars with DDC because
%   their file sizes change often. However, create_segment_raw_file_list_v2.m
%   is now the preferred function to use for these systems.
% file_size_checking = 'allbig';
% anydifference: Any difference in the file size will cause a segment
%    break to happen. This is the preferred method for all other systems.
file_size_checking = 'anydifference';

% =========================================================================
%% Automated Section
% =========================================================================

create_segment_raw_file_list

return;
