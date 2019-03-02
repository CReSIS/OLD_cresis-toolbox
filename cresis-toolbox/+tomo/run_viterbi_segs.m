%% Script run_viterbi_segs
%
% Script for running viterbi_segs
%
% Authors: John Paden, Victor Berger
%
% See also: viterbi_segs.m

clear param; clear data;

fprintf('\n\n========================================================\n');
fprintf('run viterbi segs\n');
fprintf('========================================================\n');

%% User Settings

params             = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'),'20161111_05','post');
params.cmd.generic = 1;
params.cmd.frms    = [];
options.name       = 'CSARP_post/standard';
options.compile    = false;

fprintf('\nSTART CLOCK:\n'); clock

tomo.viterbi_segs(params, options);

fprintf('\nEND CLOCK:\n'); clock