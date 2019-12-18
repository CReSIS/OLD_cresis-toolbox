function method_str = array_proc_method_str(method_int)
% method_str = array_proc_method_str(method_int)
%
% Support function for array_proc to convert array proc method integer into
% method string
%
% See also: array_proc.m, sim.doa.m, array_proc_methods.m,
% array_proc_method_str.m, array_proc_method_strs.m
%
% Author: John Paden

array_proc_methods;
switch (method_int)
  case STANDARD_METHOD
    method_str = 'standard';
  case MVDR_METHOD
    method_str = 'mvdr';
  case MVDR_ROBUST_METHOD
    method_str = 'mvdr_robust';
  case MUSIC_METHOD
    method_str = 'music';
  case EIG_METHOD
    method_str = 'eig';
  case RISR_METHOD
    method_str = 'risr';
  case GEONULL_METHOD
    method_str = 'geonull';
  case GSLC_METHOD
    method_str = 'gslc';
  case MUSIC_DOA_METHOD
    method_str = 'music_doa';
  case MLE_METHOD
    method_str = 'mle';
  case DCM_METHOD
    method_str = 'dcm';
  case PF_METHOD
    method_str = 'pf';
  otherwise
    error('Invalid method integer (%d)', method_int);
end
