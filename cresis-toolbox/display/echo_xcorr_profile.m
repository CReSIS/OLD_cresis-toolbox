function param = echo_xcorr_profile(profile_name)
% param = echo_xcorr_profile(profile_name)
%
% Fast-time cross correlation parameter profiles.
%
% INPUTS:
%
% profile_name: string containing profile name
%
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

switch (profile_name)
  case 'long_unitstep'
    param.h_filt = fliplr([0.1*ones(1,10) 0.5*ones(1,11)]);
    param.h_filt = param.h_filt-mean(param.h_filt);
    param.h_filt_offset = -8;
  case 'peaky'
    param.h_filt = [0.1 0.1 0.4 0.5 0.4 0.1 0.1];
    param.h_filt = param.h_filt-mean(param.h_filt);
    param.h_filt_offset = -3;
  case 'short_unitstep'
    param.h_filt = fliplr([0.1*ones(1,6) 0.5*ones(1,7)]);
    param.h_filt = param.h_filt-mean(param.h_filt);
    param.h_filt_offset = -4;
  otherwise
    error('Profile %s does not exist.', profile_name);
end
