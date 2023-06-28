function fname = fname_info_bas(fn)
% fname = fname_info_bas(fn)
%
% Parses the BAS .mat filename
%
% fn = BAS file name
%
% fname = structure containing each of the fields of the MCORDS
%   filename
%
%   For example: f31aBelly_20170122195130_TxP1234_RxB6_C13L4_T01_0044.mat
%  .name: 'f31a'
%  .tx_weights: transmit weights array for elements 1-12
%     TxP1234 corresponds to elements 1-4
%     TxP1234 corresponds to elements 9-12
%  .adc: scalar with ADC
%     B1,B2,B3,...,B9,BA,BB,BC correspond to adc 1 to 12)
%  .zero_pi: 0 or 180
%     0 for 'C', 180 for 'J'
%  .BW = 13e6
%     13 is 13e6
%  .Tpd = 4e-6
%     L1 is 1e-6, L4 is 4e-6
%  .file_idx = 44
%  .datenum = time stamp converted to Matlab date number. Common
%     usage: "[year month day hour min sec] = datevec(fname.datenum)"
%
% Example
%  fn = 'f31aBelly_20170122195130_TxP1234_RxB6_C13L4_T01_0044.mat';
%  fname = fname_info_bas(fn)
%
%  fn = 'f31aBelly_20170122195130_TxP1234_RxB8_C13L1_T01_0044.mat';
%  fname = fname_info_bas(fn)
%
%  fn = 'f31aBelly_20170122195130_TxS9ABC_RxB8_J13L4_T01_0044.mat';
%  fname = fname_info_bas(fn)
%
% Author: John Paden
%
% See also datenum

[path name ext] = fileparts(fn);
fn = [name ext];

[fname.name fn] = strtok(fn,'_');
if fname.name(end-3) == 'P' || fname.name(end-3) == 'S'
  fname.name = fname.name(1:end-4);
else
  fname.name = fname.name(1:end-5);
end

[fname.datenum fn] = strtok(fn,'_');
fname.datenum = datenum(fname.datenum,'yyyymmddHHMMSS');

[fname.wfs.tx_weights fn] = strtok(fn,'_');
if strcmpi(fname.wfs.tx_weights,'TxP1234')
  fname.wfs.tx_weights = [2000 2000 2000 2000 0 0 0 0 0 0 0 0];
elseif strcmpi(fname.wfs.tx_weights,'TxS9ABC')
  fname.wfs.tx_weights = [0 0 0 0 0 0 0 0 2000 2000 2000 2000];
else
  warning('Unknown transmit code.');
  keyboard
end

[fname.wfs.adc fn] = strtok(fn,'_');
fname.wfs.adc = hex2dec(fname.wfs.adc(4));

[tmp fn] = strtok(fn,'_');
if tmp(1) == 'C'
  fname.wfs.wf_adc_sum = 1;
elseif tmp(1) == 'J'
  fname.wfs.wf_adc_sum = -1;
end
fname.wfs.BW = str2double(tmp(2:3)) * 1e6;
fname.wfs.Tpd = str2double(tmp(5)) * 1e-6;

[tmp fn] = strtok(fn,'_');
[fname.file_idx fn] = strtok(fn(2:end),'.');
fname.file_idx = str2double(fname.file_idx);
