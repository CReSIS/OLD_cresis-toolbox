% script array_proc_methods
%
% Script with the integer equivalents of each method
%
% See also: array_proc.m, sim.doa.m, array_proc_methods.m,
% array_proc_method_str.m, array_proc_method_strs.m
%
% Author: John Paden

STANDARD_METHOD = 0;
MVDR_METHOD = 1;
MVDR_ROBUST_METHOD = 2;
MUSIC_METHOD = 3;
EIG_METHOD = 4;
RISR_METHOD = 5;
GEONULL_METHOD = 6;
GSLC_METHOD = 7;
DOA_METHOD_THRESHOLD = 2^16;
MUSIC_DOA_METHOD = 2^16 + 3;
MLE_METHOD = 2^16 + 7;
DCM_METHOD = 2^16 + 8;
PF_METHOD = 2^16 + 9;
WBMLE_METHOD = 2^16 + 10;
