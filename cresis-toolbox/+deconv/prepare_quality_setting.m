% script read_tom_newman_spreadsheet
%
% Must be run from windows
%
% Sets up "quality" settings field in frames file.
% Optionally reads .xls file with the following format:
%   Worksheet name should match the year of season name.
%   Format of .xls file (NA fields do not matter) should contain a series
%   of entries with the following format:
%
% COLUMN 1    | COLUMN 2                       | NA | NA | NA | COLUMN 6      | NA
% YYYYMMDD_SS | frames in matlab vector format | NA | NA | NA | Artifact code | NA
% NA          | frames in matlab vector format | NA | NA | NA | Artifact code | NA
% NA          | frames in matlab vector format | NA | NA | NA | Artifact code | NA
% NA          | empty or bad string            | NA | NA | NA | Artifact code | NA <-- TERMINATES THE CURRENT SEGMENT
%
% matlab vector format example: "23, 12:14, 74, 2 3"
% Artifact codes: see source code below to see mapping and modify as required

%% User Settings
param.radar_name = 'snow';
% param.season_name = '2009_Greenland_P3';
% param.season_name = '2010_Greenland_DC8';
% param.season_name = '2011_Greenland_P3';
% param.season_name = '2012_Greenland_P3';
% param.season_name = '2009_Antarctica_DC8';
% param.season_name = '2010_Antarctica_DC8';
% param.season_name = '2011_Antarctica_DC8';
% param.season_name = '2012_Antarctica_DC8';
param.season_name = '2017_Antarctica_P3';
% param.season_name = 'SEASON';

if 1
  % xls_fn: quality parameter spreadsheet (leave empty to not use)
%   xls_fn = 'C:\users\dangermo\Desktop\Processing\Snow Radar\snowRadarDeconQc_Newman_final.xlsx';
  xls_fn = 'Y:\2017_Antarctica_P3\2017_Antarctica_P3_snowRadar_DeconQc_snow.xlsx';
else
  xls_fn = '';
%   params = read_param_xls(ct_filename_param('snow_param_2009_Antarctica_DC8.xls'));
%     params = read_param_xls(ct_filename_param('snow_param_2010_Antarctica_DC8.xls'));
%     params = read_param_xls(ct_filename_param('snow_param_2011_Antarctica_DC8.xls'));
    params = read_param_xls(ct_filename_param('snow_param_2012_Antarctica_DC8.xls'));
end

ignore_artifact_type = []; % Do not write these types to the frames file
overwrite_mode = true; % Set to false to have a more interactive experience

%% Automated Section

if ~isempty(xls_fn)
  if ~ispc
    error('XLS reading must be run from windows');
  end
  
  [~,~,alldata] = xlsread(xls_fn,param.season_name(1:4));
else
  for param_idx = 1:length(params)
    param = params(param_idx);
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
      continue;
    end
    frames_fn = ct_filename_support(param,'','frames');
    load(frames_fn);
    
    update_field = 'quality';
    if ~isfield(frames,update_field) || overwrite_mode
      frames.(update_field) = NaN*zeros(size(frames.frame_idxs));
    end
    
    old_frames = frames;
    
    bad_data_mask = frames.proc_mode ~= 0;
    frames.quality(bad_data_mask) = 2^4; % fifth bit: 2^4
    
    unassigned_mask = isnan(frames.quality) & frames.proc_mode == 0;
    frames.quality(unassigned_mask) = 0;
    
    %% Print out frames that changed
    diff_frms = find((frames.(update_field)~=old_frames.(update_field) & ~(isnan(frames.(update_field)) & isnan(old_frames.(update_field)))) ...
      | (frames.nyquist_zone~=old_frames.nyquist_zone & ~(isnan(frames.nyquist_zone) & isnan(old_frames.nyquist_zone))));
    if ~isempty(diff_frms)
      fprintf('%-5s\t%-7s\t%-7s\t%7s\t%7s\n', 'Frm', 'Old', 'New', 'Old NZ', 'New NZ');
      for cur_frm = diff_frms
        fprintf('%03.0f  \t%04.0f   \t%04.0f   \t%7.0f\t%7.0f\n', cur_frm, ...
          old_frames.(update_field)(cur_frm), frames.(update_field)(cur_frm), ...
          old_frames.nyquist_zone(cur_frm), frames.nyquist_zone(cur_frm));
      end
    end
    
    if ~isempty(diff_frms)
      %% Save Output
      if overwrite_mode
        val = 'Y';
      else
        val = input(sprintf('Save frames.%s before quitting (Y/N)? ',update_field),'s');
      end
      if strncmpi(val,'Y',1)
        out_dir = fileparts(frames_fn);
        if ~exist(out_dir,'dir')
          mkdir(out_dir);
        end
        fprintf('  Saving frames file %s\n',frames_fn);
        save(frames_fn,'-v6','frames');
      else
        fprintf('  Not saving (can still manually save by pasting commands)\n');
      end
    end
  end
  return;
end

for row = 1:size(alldata,1)
  day_seg = alldata{row,1};
  if all(~isnan(day_seg)) && length(day_seg) == 11
    param.day_seg = day_seg;
    fprintf('=============================================================\n');
    fprintf('=============================================================\n');
    fprintf('Updating %s\n', param.day_seg);
    
    frames_fn = ct_filename_support(param,'','frames');
    
    load(frames_fn);
    
    update_field = 'quality_snow';
    if ~isfield(frames,update_field) || overwrite_mode
      frames.(update_field) = NaN*zeros(size(frames.frame_idxs));
    end
    
    old_frames = frames;
    
    bad_data_mask = frames.proc_mode ~= 0;
    frames.(update_field)(bad_data_mask) = 2^4; % fifth bit: 2^4
    
    unassigned_mask = isnan(frames.(update_field)) & frames.proc_mode == 0;
    frames.(update_field)(unassigned_mask) = 0;
    
    update_field_mask_len = 8;
    % 1: coherent noise artifact
    % 2: deconvolution artifact
    % 3: elevated noise floor or vertical stripes
    % 4: missing data
    % 5: no good data
    % 6: low SNR
    % 7: unclassified artifact (elevation change, wavy ripples)
    % 8: land or iceberg
    frms = [];
    artifact_type = NaN;
    or_mask = true;
    
    if isempty(xls_fn)
      done = true;
    else
      done = false;
    end
    while ~done
      try
        if isnan(alldata{row,2})
          frms = [];
        elseif isfloat(alldata{row,2})
          frms = alldata{row,2};
        else
          frms = eval(['[' alldata{row,2} ']']);
        end
      catch ME
        done = true;
        continue
      end
      
      if isempty(frms)
        done = true;
        continue
      end
      
      if ~isempty(regexpi(alldata{row,5},'SNR'))
        artifact_type = 6;
      elseif ~isempty(regexpi(alldata{row,5},'land')) || ~isempty(regexpi(alldata{row,5},'ramp')) || ~isempty(regexpi(alldata{row,5},'iceberg'))
        artifact_type = 8;
      elseif ~isempty(regexpi(alldata{row,5},'coher'))
        artifact_type = 1;
      elseif ~isempty(regexpi(alldata{row,5},'vert'))
        artifact_type = 3;
      elseif ~isempty(regexpi(alldata{row,5},'decon')) || ~isempty(regexpi(alldata{row,5},'lobe')) || ~isempty(regexpi(alldata{row,5},'blur'))
        artifact_type = 2;
      elseif ~isempty(regexpi(alldata{row,5},'bad'))
        if ~isempty(regexpi(alldata{row,5},'gap'))
          artifact_type = 4;
        else
          artifact_type = 5;
        end
      else
        warning('Unknown artifact type');
        if ~overwrite_mode
          keyboard
        end
        artifact_type = 7;
      end
      
      alldata{row,2}
      frms
      alldata{row,5}
      artifact_type
      
      if any(artifact_type == ignore_artifact_type) || (length(frms) == 1 && frms == 0)
        row = row + 1;
        continue
      end
      
      for frm = frms
        
        quality_mask = fliplr(dec2bin(frames.(update_field)(frm),update_field_mask_len));
        
        if or_mask
          if quality_mask(artifact_type) == '0'
            frames.(update_field)(frm) = frames.(update_field)(frm) + 2^(artifact_type-1);
          end
        else
          frames.(update_field)(frm) = 2^(artifact_type-1);
        end
        
      end
      
      row = row + 1;
    end
    
    %% Print out frames that changed
    diff_frms = find((frames.(update_field)~=old_frames.(update_field) & ~(isnan(frames.(update_field)) & isnan(old_frames.(update_field)))) ...
      | (frames.nyquist_zone~=old_frames.nyquist_zone & ~(isnan(frames.nyquist_zone) & isnan(old_frames.nyquist_zone))));
    if ~isempty(diff_frms)
      fprintf('%-5s\t%-7s\t%-7s\t%7s\t%7s\n', 'Frm', 'Old', 'New', 'Old NZ', 'New NZ');
      for cur_frm = diff_frms
        fprintf('%03.0f  \t%04.0f   \t%04.0f   \t%7.0f\t%7.0f\n', cur_frm, ...
          old_frames.(update_field)(cur_frm), frames.(update_field)(cur_frm), ...
          old_frames.nyquist_zone(cur_frm), frames.nyquist_zone(cur_frm));
      end
    end
    
    if ~isempty(diff_frms)
      %% Save Output
      if overwrite_mode
        val = 'Y';
      else
        val = input(sprintf('Save frames.%s before quitting (Y/N)? ',update_field),'s');
      end
      if strncmpi(val,'Y',1)
        out_dir = fileparts(frames_fn);
        if ~exist(out_dir,'dir')
          mkdir(out_dir);
        end
        fprintf('  Saving frames file %s\n',frames_fn);
        save(frames_fn,'-v6','frames');
      else
        fprintf('  Not saving (can still manually save by pasting commands)\n');
      end
    end
    
  end
end

return

