function param = echo_xcorr_profile(profile_name)
% param = echo_xcorr_profile(profile_name)
%
% Fast-time cross correlation parameter profiles.
%
% INPUTS:
%
% profile_name: string containing profile name
%
% OUTPUTS:
%
% param: parameters for the selected profile
%
% Examples:
%
% See echo_xcorr.m
%
% Author: John Paden
%
% See also: echo_detrend, echo_filt, echo_mult_suppress, echo_noise,
% echo_norm, echo_param, echo_stats, echo_stats_layer, echo_xcorr,
% echo_xcorr_profile

switch (profile_name)
  case 'short_unitstep'
    param.h_filt = fliplr([0.1*ones(1,6) 0.5*ones(1,7)]);
    param.h_filt = param.h_filt-mean(param.h_filt);
    param.h_filt_offset = -4;
  case 'long_unitstep'
    param.h_filt = fliplr([0.1*ones(1,10) 0.5*ones(1,11)]);
    param.h_filt = param.h_filt-mean(param.h_filt);
    param.h_filt_offset = -8;
  case 'xlong_unitstep'
    param.h_filt = fliplr([0.1*ones(1,20) 0.5*ones(1,21)]);
    param.h_filt = param.h_filt-mean(param.h_filt);
    param.h_filt_offset = -16;
  case 'xlong_unitstep_delay'
    param.h_filt = fliplr([0.1*ones(1,20) 0.5*ones(1,21)]);
    param.h_filt = param.h_filt-mean(param.h_filt);
    param.h_filt_offset = -8;
  case 'snow'
    param.h_filt = fliplr([0.1*ones(1,21) 4*ones(1,5) 0.5*ones(1,21)]);
    param.h_filt = param.h_filt-mean(param.h_filt);
    param.h_filt_offset = -20; %-48
  case 'peaky'
    param.h_filt = [0.1 0.1 0.4 0.5 0.4 0.1 0.1];
    param.h_filt = param.h_filt-mean(param.h_filt);
    param.h_filt_offset = -3;
  otherwise
    error('Profile %s does not exist.', profile_name);
end
