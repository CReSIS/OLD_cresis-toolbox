function [hdr,data,fn_name] = noise_meas_2012_Antarctica_DC8_func(param)
%
% This script evaluates 50 ohm term and receive only data.
% It loads all receive channels and then analyzes these data.
%
% For P-3:
% 1. Collect data with ? waveforms in whatever noise configuration
%    you want to measure.
% 2. Characterization should be done for 50 ohm and receive only at
%    least.
% 3. If transmitting, time gate should be large enough to include
%    noise-only data.
%
% Author: John Paden

physical_constants;
tstart = tic;

% =======================================================================
% =======================================================================
% Automated Section
% =======================================================================
% =======================================================================

% =======================================================================
% Load data
% =======================================================================
% fprintf('========================================================\n');
% fprintf('Loading data (%.1f sec)\n', toc(tstart));
% Load the data (disable if you have already loaded)
clear data;
clear num_rec;
if strcmpi(param.radar_name,'mcords')
  for adc_idx = 1:length(param.adcs)
    adc = param.adcs(adc_idx);

    % May need to adjust base_path for non-standard directory structures
    base_path = fullfile(param.base_path, sprintf('chan%d',adc), ...
      param.seg);
    file_midfix = sprintf('r%d-%d.',param.radar_num,adc);
    file_suffix = sprintf('.%04d.dat',param.data_file_num);
%     fprintf('  Path: %s\n', base_path);
%     fprintf('  Match: mcords*%s*%s\n', file_midfix, file_suffix);
    fn = get_filename(base_path,'mcords',file_midfix,file_suffix);
    if isempty(fn)
      fprintf('  Could not find any files which match\n');
      return;
    end
%     fprintf('  Loading file %s\n', fn);
    [hdr,data_tmp] = basic_load_mcords(fn, struct('clk',1e9/9,'first_byte',2^26));
    data(:,:,adc_idx) = data_tmp{param.wf}(1:end-1,1:min(size(data_tmp{param.wf},2),param.rlines(2)));
  end
  data = data - median(data(:,1));
%   basic_remove_mcords_digital_errors;
elseif strcmpi(param.radar_name,'mcords2')
  % test1_1.dat0
  %   testA_N.datB
  %   A = acquisition number
  %   N = file number
  %   B = board number
  % Map ADCs to board numbers
  for board = 0:3
    if any(board == floor((param.adcs-1)/4))
      get_adcs = board*4 + (1:4);
      file_prefix = sprintf('mcords2_%d_',board);
      if isempty(param.acquisition_num)
        file_suffix = sprintf('%04d.bin',param.file_num);
      else
        file_suffix = sprintf('%02d_%04d.bin',param.acquisition_num,param.file_num);
      end
      base_path = fullfile(param.base_path, sprintf('board%d',board), ...
        param.seg);
%       fprintf('  Path: %s\n', base_path);
%       fprintf('  Match: %s*%s\n', file_prefix, file_suffix);
      fn = get_filenames(base_path, file_prefix, '', file_suffix);
      if isempty(fn)
        fprintf('  Could not find any files which match\n');
        return;
      elseif length(fn) == 1
        fn = fn{1};
      else
        fprintf('Select a specific file (fn = fn{1}), then dbcont\n');
        fn
        keyboard
      end
%       fprintf('  Loading file %s\n', fn);
      % Fix get_filenames     'The filename, directory name, or volume label syntax is incorrect.'
      [hdr,data_tmp] = basic_load_mcords2(fn,struct('clk',param.fs));
      for get_adc_idx = 1:length(get_adcs)
        for adc_idx = 1:length(param.adcs)
          if param.adcs(adc_idx) == get_adcs(get_adc_idx)
            if ~exist('num_rec','var')
              % Since each file may have slightly different numbers of
              % records we do this
              num_rec = size(data_tmp{param.wf},2) - 1;
            end
            data(:,:,adc_idx) = data_tmp{param.wf}(:,1:num_rec,get_adc_idx);
          end
        end
      end
    end
  end
end
[fn_dir fn_name] = fileparts(fn);
if ~isfield(param,'seg') || isempty(param.seg)
  param.seg = -1;
end
clear data_tmp;


if ~isfield(param,'rlines') || isempty(param.rlines)
  rlines = 1:size(data,2);
elseif param.rlines(end) > size(data,2)
  rlines = param.rlines(1):size(data,2);
else
  rlines = param.rlines(1):param.rlines(end);
end
if ~isfield(param,'noise_rbins') || isempty(param.noise_rbins)
  noise_rbins = 1:size(data,1);
elseif param.noise_rbins(end) > size(data,1)
  noise_rbins = param.noise_rbins(1):size(data,1);
else
  noise_rbins = param.noise_rbins(1):param.noise_rbins(end);
end

return;
