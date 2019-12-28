function method_strs = array_proc_method_strs(method_ints)
% method_strs = array_proc_method_strs(method_ints)
%
% Support function for generating legends for plotting when multiple
% methods are used. sim.doa_example_* functions use this.
%
% method_ints: list of array processing methods in integer form
%
% method_strs: if method_ints has one element then method_strs returns a
% string, if method_ints has multiple elements, then method_strs returns a
% cell array
%
% See also: array_proc.m, sim.doa.m, array_proc_methods.m,
% array_proc_method_str.m, array_proc_method_strs.m
%
% Author: John Paden

array_proc_methods;

method_strs = cell(size(method_ints));
for method_idx = 1:numel(method_ints)
  method_int = method_ints(method_idx);
  switch (method_int)
    case STANDARD_METHOD
      method_strs{method_idx} = 'standard';
    case MVDR_METHOD
      method_strs{method_idx} = 'mvdr';
    case MVDR_ROBUST_METHOD
      method_strs{method_idx} = 'mvdr_robust';
    case MUSIC_METHOD
      method_strs{method_idx} = 'music';
    case EIG_METHOD
      method_strs{method_idx} = 'eig';
    case RISR_METHOD
      method_strs{method_idx} = 'risr';
    case GEONULL_METHOD
      method_strs{method_idx} = 'geonull';
    case GSLC_METHOD
      method_strs{method_idx} = 'gslc';
    case MUSIC_DOA_METHOD
      method_strs{method_idx} = 'music_doa';
    case MLE_METHOD
      method_strs{method_idx} = 'mle';
    case DCM_METHOD
      method_strs{method_idx} = 'dcm';
    case PF_METHOD
      method_strs{method_idx} = 'pf';
    case DOA_TAGGING_METHOD
      method_strs{method_idx} = 'doa_tag';
    otherwise
      error('Invalid method integer (%d)', method_int);
  end
end
