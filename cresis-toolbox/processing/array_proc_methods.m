% script array_proc_methods
%
% Script with the integer equivalents of each method

% Beamformer Methods
STANDARD_METHOD = 0;
MVDR_METHOD = 1;
MVDR_ROBUST_METHOD = 2;
MUSIC_METHOD = 3;
EIG_METHOD = 4;
RISR_METHOD = 5;
GEONULL_METHOD = 6;

% DOA Methods
DOA_METHOD_THRESHOLD = 2^16;
MUSIC_DOA_METHOD = 2^16 + 3;
MLE_METHOD = 2^16 + 7;
DCM_METHOD = 2^16 + 8;
PF_METHOD = 2^16 + 9;
