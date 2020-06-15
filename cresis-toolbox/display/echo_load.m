function [mdata,fn_echo] = echo_load(param,echogram_source,frm,img)
% [mdata,fn_echo] = echo_load(param,echogram_source,frm,img)
%
% The trend of the data is estimated using various methods and this trend
% is removed from the data.
%
% INPUTS:
%
% param: radar parameter spreadsheet structure
% echogram_source: 'standard'
%
% param: struct controlling how the normalization is done
%
% OUTPUTS:
%
% mdata: echogram structure
% fn_echo: string containing the filename loaded
%
% Examples:
%
% param = read_param_xls(ct_filename_param('rds_param_2019_Antarctica_Ground.xls'),'20200107_01');
% mdata = echo_load(param,'standard',1);
%
% Author: John Paden
%
% See also: echo_detrend, echo_filt, echo_mult_suppress, echo_noise,
% echo_norm, echo_param, echo_stats, echo_stats_layer, echo_xcorr,
% echo_xcorr_profile

if ischar(param)
  % param is a filename (second and third arguments are ignored)
  fn_echo = param;
  
elseif isstruct(param)
  % param is a parameter spreadsheet structure
  
  % Override default echogram source if second argument defined
  if ~exist('echogram_source','var') || isempty(echogram_source)
    echogram_source = 'standard';
  end
  
  % Override default img if third argument defined
  if ~exist('img','var') || isempty(img)
    % img == 0 is combined
    img = 0;
  end
  
  % Override frames field if third argument defined
  if exist('frm','var') && ~isempty(frm)
    param.cmd.frms = frm;
  elseif isempty(param.cmd.frms)
    param.cmd.frms = 1;
  else
    param.cmd.frms = param.cmd.frms(1);
  end
  
  if img == 0
    % Combined file
    fn_echo = fullfile(ct_filename_out(param,echogram_source),sprintf('Data_%s_%03d.mat',param.day_seg,param.cmd.frms));
  else
    % Image "_img_II" file
    fn_echo = fullfile(ct_filename_out(param,echogram_source),sprintf('Data_img_%02d_%s_%03d.mat',img, param.day_seg,param.cmd.frms));
  end
  
end
mdata = load_L1B(fn_echo);
