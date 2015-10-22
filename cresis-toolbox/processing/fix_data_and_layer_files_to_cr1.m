% Fix data and layers files to the new CR1 format
%
% Converts old files to the new CR1 format
% 1. First create records files
% 2. Run this file which will copy frames, echograms, and layer data

param_fn = '/users/paden/scripts/branch/params-cr1/rds_param_2008_Greenland_TO.xls';

old_records_dir = '/cresis/projects/dev/csarp_support/records/mcrds/2008_Greenland_TO/'
old_frames_dir = '/cresis/projects/dev/csarp_support/frames/mcrds/2008_Greenland_TO/'
old_data_dir = '/cresis/scratch2/mdce/mcrds/2008_Greenland_TO/'

% =========================================================================
%% Automated Section

params = read_param_xls(param_fn,[],'post');

if 0
  for param_idx = 1:length(params)
    param = params(param_idx);
    param.day_seg_old = param.day_seg;
    segment_offset = 0;
    param.day_seg_old(end-1:end) = sprintf('%02d', str2double(param.day_seg_old(end-1:end)) + segment_offset);
    
    if param.cmd.generic
      fprintf('=======================================================\n');
      fprintf('Fixing %s (%s)\n', param.day_seg, param.day_seg_old);
      fprintf('=======================================================\n');
      
      if ~isfield(param.vectors,'tmp') || isempty(param.vectors.tmp)
        param.vectors.tmp = 0;
      end
      
      records_fn = ct_filename_support(param,'','records');
      records = load(records_fn);
      
      fprintf(' NEW: %s to %s\n', datestr(epoch_to_datenum(records.gps_time(1))), ...
        datestr(epoch_to_datenum(records.gps_time(end))));
      
      frames_fn = ct_filename_support(param,'','frames');
      
      old_frames_fn = fullfile(old_frames_dir, ...
        sprintf('frames_%s.mat', param.day_seg_old));
      if exist(old_frames_fn,'file')
        load(old_frames_fn);
      else
        keyboard
      end
      
      old_records_fn = fullfile(old_records_dir, ...
        sprintf('records_%s.mat', param.day_seg_old));
      if ~exist(old_records_fn,'file')
        %% Old records do not exist, use layerData to make frames
        old_base_dir = fullfile(old_data_dir,'CSARP_layerData/',param.day_seg_old);
        old_frames = frames;
        
        for frm = 1:length(frames.frame_idxs)-1
          old_fn = fullfile(old_base_dir,sprintf('Data_%s_%03d.mat', ...
            param.day_seg_old, frm));
          old_gps_time = load(old_fn,'GPS_time');
          old_fn = fullfile(old_base_dir,sprintf('Data_%s_%03d.mat', ...
            param.day_seg_old, frm+1));
          old_gps_time2 = load(old_fn,'GPS_time');
          old_gps_time = (old_gps_time.GPS_time(end) + old_gps_time2.GPS_time(1))/2;
          frames.frame_idxs(frm+1) = find(records.gps_time > old_gps_time,1);
        end
        
        frames.frame_idxs - old_frames.frame_idxs
        
        if all(abs(frames.frame_idxs - old_frames.frame_idxs) < 1500)
          fprintf('  Using old frames\n');
          frames = old_frames;
        else
          keyboard
        end
        if ~isfield(frames,'proc_mode')
          frames.proc_mode = zeros(size(frames.frame_idxs));
        end
        
        save(frames_fn,'frames');
        
      else
        %% Old records file exists (make frames using it)
        old_records = load(old_records_fn);
        if ~isfield(old_records,'records')
          old_records.records.epri = old_records.raw.epri;
        end
        
        if max(frames.frame_idxs) > length(records.lat)
          keyboard;
        end
        along_track = geodetic_to_along_track(records.lat,records.lon,records.elev);
        frame_lengths = diff(along_track(frames.frame_idxs));
        fprintf('  Frames %.1f to %.1f km\n', min(frame_lengths)/1e3, max(frame_lengths)/1e3);
        
        frames_fn_dir = fileparts(frames_fn);
        if ~exist(frames_fn_dir,'dir')
          mkdir(frames_fn_dir);
        end
        
        if frames.frame_idxs(1) ~= 1
          keyboard
        end
        if isfield(old_records.records,'epri')
          % Update frames file
          for idx = 2:length(frames.frame_idxs)
            start_epri = old_records.records.epri(frames.frame_idxs(idx));
            if iscell(records.raw.epri)
              new_record = find(records.raw.epri{1} == start_epri);
            else
              new_record = find(records.raw.epri == start_epri);
            end
            if isempty(new_record)
              keyboard
              if iscell(records.raw.epri)
                new_record = find(records.raw.epri{1} >= start_epri);
              else
                new_record = find(records.raw.epri >= start_epri);
              end
            end
            if new_record ~= frames.frame_idxs(idx)
              fprintf('Updating frame %d\n', idx);
            end
            frames.frame_idxs(idx) == new_record;
          end
          
        else
          % Update frames file
          for idx = 2:length(frames.frame_idxs)
            start_epri = old_records.records.radar_time(frames.frame_idxs(idx));
            new_record = find(records.raw.radar_time == start_epri);
            if isempty(new_record)
              warning('No exact match on radar time for frame %s, finding closest match', idx);
              new_record = find(records.raw.radar_time >= start_epri,1);
            end
            if new_record ~= frames.frame_idxs(idx)
              fprintf('Updating frame %d\n', idx);
            end
            frames.frame_idxs(idx) == new_record;
          end
        end
        
        if ~isfield(frames,'proc_mode')
          frames.proc_mode = zeros(size(frame_idxs));
        end
        save(frames_fn,'frames');
        %system(sprintf('cp %s %s',old_frames_fn,frames_fn));
      end
      
      data_out{1} = ct_filename_out(param,'','CSARP_standard');
      data_out{2} = ct_filename_out(param,'','CSARP_qlook');
      data_out{3} = ct_filename_out(param,'','CSARP_mvdr');
      data_out{4} = ct_filename_out(param,'','CSARP_layerData');
      
      data_dirs{1} = fullfile(old_data_dir,'CSARP_standard/',param.day_seg_old);
      data_dirs{2} = fullfile(old_data_dir,'CSARP_qlook/',param.day_seg_old);
      data_dirs{3} = fullfile(old_data_dir,'CSARP_mvdr/',param.day_seg_old);
      data_dirs{4} = fullfile(old_data_dir,'CSARP_layerData/',param.day_seg_old);
      
      for dir_idx = 1:length(data_dirs)
        if ~exist(data_dirs{dir_idx},'dir')
          continue
        end
        fprintf(' Copying %s\n', data_dirs{dir_idx});
        fns = get_filenames(data_dirs{dir_idx},sprintf('Data_%s',param.day_seg_old),'','.mat');
        using_img01 = false;
        if isempty(fns)
          fns = get_filenames(data_dirs{dir_idx},sprintf('Data_img_01_%s',param.day_seg_old),'','.mat');
          if isempty(fns)
            continue
          else
            using_img01 = true;
          end
        end
        first = load(fns{1},'GPS_time');
        last = load(fns{end},'GPS_time');
        
        fprintf(' OLD: %s to %s\n', datestr(epoch_to_datenum(first.GPS_time(1))), ...
          datestr(epoch_to_datenum(last.GPS_time(end))));
        
        % This code checks to make sure the segment GPS time length has
        % not changed
        if first.GPS_time(1) < records.gps_time(1)
          fprintf('  NONOVERLAPPING_TIME\n');
          if first.GPS_time(1) < records.gps_time(1) - 25
            keyboard
          end
        end
        if last.GPS_time(end) > records.gps_time(end)
          fprintf('  NONOVERLAPPING_TIME\n');
          if last.GPS_time(end) > records.gps_time(end) + 25
            keyboard
          end
        end
        
        if ~exist(data_out{dir_idx},'dir')
          mkdir(data_out{dir_idx});
        end
        % Copy Data Files
        for fn_idx = 1:length(fns)
          fn = fns{fn_idx};
          [fn_dir fn_name fn_ext] = fileparts(fn);
          fn_name(end-14:end-4) = param.day_seg;
          if ~using_img01
            out_fn = fullfile(data_out{dir_idx},[fn_name fn_ext]);
          else
            out_fn = fullfile(data_out{dir_idx},[fn_name([1:5 13:end]) fn_ext]);
          end
          fprintf(' Copying %s to %s\n', fn, out_fn);
          sys_cmd = sprintf('cp %s %s', fn, out_fn);
          system(sys_cmd);
          % copyfile(fn,out_fn);
        end
      end
    end
    
  end
  
end

% =========================================================================
%% Update records with layer information so that the surface and bottom
% in the records file reflect the layerData contents
fprintf('Do you want to update records?');
keyboard

update_records_param.param_fn = param_fn;

% param.skip_phrase: All segments will be skipped with this phrase in their
%   verification field. "do not process" is the standard. Leave this field
%   blank to do all segments.
update_records_param.skip_phrase = 'do not process';

% param.save_changes: Logical, For debugging purposes, you can turn the file save on/off
update_records_param.save_changes = true;

% param.debug_level: default is 1, anything higher is for debugging
update_records_param.debug_level = 1;

update_records_param.layers.update_en = true;
update_records_param.layers.path = '';

update_records_param.gps_time.update_en = false;
update_records_param.gps_time.path = '';
update_records_param.gps_time.time_offset = -14;
update_records_param.gps_time.special_en = false;

update_records_param.gps_source.update_en = false;
update_records_param.gps_source.path = '';

update_records(update_records_param);


