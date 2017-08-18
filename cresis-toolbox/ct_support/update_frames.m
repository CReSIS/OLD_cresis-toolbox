% script update_frames
%
% Creates frames files for FMCW and accumulation radar. Requires
% that the records file be created.
%
% See run_update_frames.m to see how to run.
%
% Prompt
%   F(FFF:PPPP:N[D]):
%     FFF = frame ID
%     PPPP = processing type
%     N = nyquist zone
%     D = default nyquist zone if 'D' present
%
% Commands (case-insensitive except where noted)
% empty or 'n'
%   proceed to next frame
% 'p'
%   proceed to previous frame
% #
%   sets frame processing type to # and proceeds to next frame
% n#
%   sets nyquist zone to # (N# proceeds to next frames after setting)
% c or k
%   runs Matlab eval command in debug prompt
% d
%   prints out differences to the frames file since last save
% i
%   swap image type between mat and param.img_type
% s
%   saves frames file
% q
%   quits (asks for saving)
% j
%   asks which frame to jump to
% Ctrl-C
%   quits without saving
%
% Author: John Paden


%% Automated Section

if strcmpi(update_field_type,'double')
elseif strcmpi(update_field_type,'mask')
  update_field_mask_len = max(3,length(update_field_mask)-1);
else
  error('Invalid field type update_field_type=%s', update_field_type);
end

records_fn = ct_filename_support(param,'','records');

img_type = param.img_type;

num_img = length(param.image_out_dir);

for img = 1:num_img
  figure(img); clf;
  h_axes(img) = axes;
  set(h_axes(img),'Position',[0 0 1 1]);
  set(h_axes(img),'TickDir','out');
  h_image(img) = imagesc(1,'Parent',h_axes(img));
  hold on;
  h_plot(img) = plot(NaN,'k','Parent',h_axes(img));
  hold off;
end

fprintf('Loading records file %s\n',records_fn);
records = load(records_fn,'relative_filename','relative_rec_num');
[records_fn_dir,records_fn_name] = fileparts(records_fn);
param.day_seg = records_fn_name(9:end);

for img=1:num_img
  image_fn_dir{img} = ct_filename_out(param,param.image_out_dir{img}{:});
  mat_fn_dir{img} = ct_filename_out(param,param.mat_out_dir{img}{:});
end

frames_fn = ct_filename_support(param,'','frames');

if exist(frames_fn,'file')
  fprintf('Loading frames file %s\n', frames_fn);
  load(frames_fn);
  if ~isfield(frames,'nyquist_zone')
    frames.nyquist_zone = NaN*zeros(size(frames.frame_idxs));
  end
  if ~isfield(frames,update_field)
    frames.(update_field) = NaN*zeros(size(frames.frame_idxs));
  end
else
  error('No frames files exists yet\n');
end
old_frames = frames;

frm_idx = 1;
if isempty(update_field_match)
  frms = 1:length(frames.frame_idxs);
else
  mask = zeros(1,length(frames.frame_idxs));
  while all(mask==0)
    for idx=1:length(update_field_match)
      if strcmpi(update_field_type,'double')
        mask = mask | frames.(update_field) == update_field_match(idx);
      elseif strcmpi(update_field_type,'mask')
        if update_field_match(idx) == 0
          mask = mask | frames.(update_field) == update_field_match(idx);
        else
          mask = mask | mod(floor(frames.(update_field) / 2^(update_field_match(idx)-1)),2);
        end
      end
    end
    if all(mask==0)
      warning('No frames match. Setting to all frames.');
      mask = ones(1,length(frames.frame_idxs));
    end
  end
  frms = find(mask);
end

frm_idx = 1;
quit_cmd = false;
while ~quit_cmd
  if frm_idx < 1
    frm_idx = 1;
  end
  if frm_idx > length(frms)
    frm_idx = length(frms);
  end
  frm = frms(frm_idx);
  
  % Update frames list if required
  if ~isempty(update_field_match)
    mask = zeros(1,length(frames.frame_idxs));
    while all(mask==0)
      for idx=1:length(update_field_match)
        if strcmpi(update_field_type,'double')
          mask = mask | frames.(update_field) == update_field_match(idx);
        elseif strcmpi(update_field_type,'mask')
          if update_field_match(idx) == 0
            mask = mask | frames.(update_field) == update_field_match(idx);
          else
            mask = mask | mod(floor(frames.(update_field) / 2^(update_field_match(idx)-1)),2);
          end
        end
      end
      if all(mask==0)
        warning('No frames match. Setting to all frames.');
        mask = ones(1,length(frames.frame_idxs));
      end
    end
    frms = find(mask);
    [~,frm_idx] = min(abs(frm-frms));
    frm = frms(frm_idx);
  end
  frm_id = sprintf('%s_%03i',param.day_seg,frm);
  
  for img = 1:num_img
    try
      if strcmpi(img_type,'mat')
        echo_fn{img} = get_filename(mat_fn_dir{img},'Data_',frm_id,'.mat');
      else
        echo_fn{img} = get_filename(image_fn_dir{img},frm_id,'',img_type);
      end
    catch ME
      if strcmpi(ME.identifier,'get_filename:no_file_found')
        if strcmpi(img_type,'mat')
          warning('Missing frame %s in %s', frm_id, mat_fn_dir{img});
        else
          warning('Missing frame %s in %s', frm_id, image_fn_dir{img});
        end
        echo_fn{img} = '';
      else
        rethrow(ME);
      end
    end
    
    try
      set(1,'Name',sprintf('%s_%03d',param.day_seg,frm));
      if ~isempty(echo_fn{img})
        fprintf('  Loading output %s\n', echo_fn{img});
        if strcmpi(img_type,'mat')
          echo = load(echo_fn{img});
          set(h_image(img),'XData',1:size(echo.Data,2));
          set(h_image(img),'YData',echo.Time);
          set(h_image(img),'CData',lp(echo.Data));
          clims = [min(lp(echo.Data(isfinite(lp(echo.Data(:)))))) max(lp(echo.Data(isfinite(lp(echo.Data(:))))))];
          if length(clims) == 2
            caxis(h_axes(img),clims);
          end
          xlim(h_axes(img),[1 size(echo.Data,2)]);
          ylim(h_axes(img),[echo.Time([1 end])]);
          set(h_plot(img),'XData',1:size(echo.Data,2));
          set(h_plot(img),'YData',echo.Surface);
          set(h_axes(img),'YDir','reverse');
          if fmcw_img_debug_mode
            for rline=1:size(echo.Data,2)
              echo.Data(lp(echo.Data(:,rline))<max(lp(echo.Data(:,rline)))+img_sidelobe,rline) = 0;
            end
            set(h_image(img),'CData',lp(echo.Data));
            Nx = size(echo.Data,2);
            max_noise = zeros(1,Nx);
            dt = echo.Time(2)-echo.Time(1);
            for idx = 1:Nx
              surf_t = echo.Surface(1,idx);
              surf_bin = (surf_t - echo.Time(1))/dt + 1;
              noise_surf = max(1, round(surf_bin - noise_time_buffer/dt)+1);
              noise_end = round(noise_surf + noise_time_duration/dt) + 1;
              max_noise(idx) = max(echo.Data(noise_surf:noise_end,idx));
            end
            max_val = lp(max(echo.Data));
            noise_threshold = median(max_noise)+noise_threshold_offset_dB;
            threshold = max(max_val+img_sidelobe,...
              noise_threshold*ones(size(max_noise)));
            if isfinite(noise_threshold)
              img_caxis = noise_threshold + [-6 +12];
              set(1,'Name',sprintf('%s_%03d: caxis threshold %.1f',param.day_seg,frm,noise_threshold));
              caxis(h_axes(img),img_caxis);
              img_cmap = [gray(64); flipud(hsv(128))];
              colormap(img_cmap)
            end
            if 0 && any(medfilt1(double(max_noise > threshold),5))
              figure(2); clf;
              plot(max_noise);
              hold on;
              plot(max_val,'r');
              plot(threshold,'k','LineWidth',2);
              hold off;
              beep
            end
            if 0
              figure(2); clf;
              plot(max_noise);
              hold on;
              plot(max_val,'r');
              plot(threshold,'k','LineWidth',2);
              hold off;
            end
          end
        else
          A = imread(echo_fn{img});
          set(h_image(img),'CData',A);
          set(h_image(img),'XData',1:size(A,2));
          set(h_image(img),'YData',1:size(A,1));
          xlim(h_axes(img),[1 size(A,2)]);
          ylim(h_axes(img),[1 size(A,1)]);
          set(h_plot(img),'XData',NaN);
          set(h_plot(img),'YData',NaN);
          set(h_axes(img),'YDir','reverse');
        end
        %     set(h_image(img),'XData',1:size(A,2));
        %     set(h_image(img),'YData',1:size(A,1));
        axis(h_axes(img),'tight');
        zoom(h_axes(img),'reset');
      end
    catch ME
      warning('Failed to load %s:\n%s', echo_fn{img}, ME.getReport());
      echo_fn{img} = '';
    end
  end
  
  if strcmpi(update_field_type,'double')
    update_field_str = sprintf('%04i',frames.(update_field)(frm));
  elseif strcmpi(update_field_type,'mask')
    if isnan(frames.(update_field)(frm))
      update_field_str = sprintf(sprintf('%%%ds',update_field_mask_len),'NaN');
    else
      update_field_str = dec2bin(frames.(update_field)(frm),update_field_mask_len);
    end
  end
  if audio_tone_for_nonzero_nonisnan && eval(audio_tone_check_code)
    beep;
  end
  
  if strcmpi(img_type,'mat')
    if isnan(frames.nyquist_zone(frm))
      val = input(sprintf('F(%03i:%s:%01iD): ', frm, ...
        update_field_str, echo.param_get_heights.radar.wfs(1).nyquist_zone),'s');
    else
      val = input(sprintf('F(%03i:%s:%01i): ', frm, ...
        update_field_str, frames.nyquist_zone(frm)),'s');
    end
  else
    if isnan(frames.nyquist_zone(frm))
      val = input(sprintf('F(%03i:%s:-D): ', frm, ...
        update_field_str),'s');
    else
      val = input(sprintf('F(%03i:%s:%01i): ', frm, ...
        update_field_str, frames.nyquist_zone(frm)),'s');
    end
  end
  if all(isstrprop(val,'wspace'))
    % Go to next frame
    frm_idx = frm_idx + 1;
  else
    num_val = str2double(val);
    
    if any(strcmpi(val,{'c','k'}))
      fprintf('Variables are frames.%s and frames.nyquist_zone. dbcont when finished.\n', update_field);
      if strcmpi(update_field_type,'mask')
        fprintf('  To update frames.%s, consider using bitor. E.g. bitor(frames.%s,bin2dec(''1000'')).\n',update_field,update_field);
      end
      keyboard
      
    elseif strcmpi(val,'d')
        diff_frms = find((frames.(update_field)~=old_frames.(update_field) & ~(isnan(frames.(update_field)) & isnan(old_frames.(update_field)))) ...
          | (frames.nyquist_zone~=old_frames.nyquist_zone & ~(isnan(frames.nyquist_zone) & isnan(old_frames.nyquist_zone))));
        fprintf('%-5s\t%-7s\t%-7s\t%7s\t%7s\n', 'Frm', 'Old', 'New', 'Old NZ', 'New NZ');
        for cur_frm = diff_frms
          fprintf('%03.0f  \t%04.0f   \t%04.0f   \t%7.0f\t%7.0f\n', cur_frm, ...
            old_frames.(update_field)(cur_frm), frames.(update_field)(cur_frm), ...
          old_frames.nyquist_zone(cur_frm), frames.nyquist_zone(cur_frm));
      end
      
    elseif strcmpi(val,'i')
      if strcmpi(param.img_type,img_type)
        img_type = 'mat';
      else
        img_type = param.img_type;
      end
      
    elseif ~isempty(strfind(upper(val),'J'))
      try
        val = input('    Enter frame to jump to: ');
      catch
        val = [];
      end
      if isempty(val)
        % Ignore command
      else
        [~,new_frm_idx] = min(abs(val-frms));
        if ~isempty(new_frm_idx)
          frm_idx = new_frm_idx;
        end
      end
      
    elseif ~isnan(num_val)
      if strcmpi(update_field_type,'double')
        fprintf('    Changing frames.%s(%d) to %d\n', update_field, frm, num_val);
      elseif strcmpi(update_field_type,'mask')
        if isnan(num_val)
          fprintf('    Changing frames.%s(%d) to %d\n', update_field, frm, num_val);
        else
          num_val = str2double(char(unique(double(sprintf('%d',num_val)))));
          if num_val ~= 0
            num_val = sum(2.^(double(sprintf('%d',num_val))-48-1));
          end
          fprintf('    Changing frames.%s(%d) to b%s (d%d)\n', update_field, ...
            frm, dec2bin(num_val,update_field_mask_len), num_val);
        end
      end
      frames.(update_field)(frm) = num_val;
      frm_idx = frm_idx + 1;
      
    elseif strncmpi(val,'n',1)
      num_val = str2double(val(2:end));
      fprintf('    Changing frames.nyquist_zone to %.0f\n', num_val);
      frames.nyquist_zone(frm) = num_val;
      if val(1) == 'N'
        frm_idx = frm_idx + 1;
      end
      
    elseif strcmpi(val,'p')
      frm_idx = frm_idx - 1;
      
    elseif strcmpi(val,'q')
      quit_cmd = true;
      fprintf('Quitting\n');
      break;
      
    elseif strcmpi(val,'s')
      out_dir = fileparts(frames_fn);
      if ~exist(out_dir,'dir')
        mkdir(out_dir);
      end
      fprintf('  Saving frames file %s\n',frames_fn);
      save(frames_fn,'-v6','frames');
      old_frames = frames;
      
    elseif strcmpi(val,'?')
      fprintf(' <enter>: go to next frame\n');
      fprintf(' p: previous frame\n');
      fprintf(' n#: set nyquist zone to # (e.g. "n3")\n');
      fprintf(' i: swap image type between mat and %s\n', param.img_type);
      fprintf(' d: print out changes to frames so far\n');
      fprintf(' j: jump to frame\n');
      fprintf(' c or k: enter matlab command\n');
      fprintf(' s: save frames file\n');
      fprintf(' q: quit\n');
      if strcmpi(update_field_type,'double')
        fprintf(' #: set frame.%s to #\n', update_field);
      else
        fprintf(' #[#...#]: sets frame.%s masks according to:\n', update_field);
        for idx=2:update_field_mask_len+1
          fprintf('   %d: %s\n',idx-1,update_field_mask{idx});
        end
      end
    end
  end
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
  val = input(sprintf('Save frames.%s before quitting (Y/N)? ',update_field),'s');
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

return;
