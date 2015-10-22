% script update_frames
%
% Creates frames files for FMCW and accumulation radar. Requires
% that the records file be created.
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
% s
%   saves frames file
% q
%   quits (asks for saving)
% j
%   asks which frame to jump to
% c
%   runs Matlab eval command
% Ctrl-C
%   quits without saving
%
% Author: John Paden

clear param;
param.radar_name = 'snow';
param.season_name = '2012_Greenland_P3';
param.day_seg = '20120327_01';
param.post.out_dir = '';
param.post.img_type = 'mat';

records_fn_wild = ct_filename_support(param,'','records');
[records_fn_dir,records_fn_name,records_fn_ext] = fileparts(records_fn_wild);
records_fns = get_filenames(records_fn_dir,[records_fn_name records_fn_ext],'','');

figure(1); clf;
h_axes = axes;
set(h_axes,'Position',[0 0 1 1]);
set(h_axes,'TickDir','out');
h_image = imagesc(1,'Parent',h_axes);
if strcmpi(param.post.img_type,'mat')
  hold on;
  h_plot = plot(1,'k');
  hold off;
end
for rec_idx = 1:length(records_fns)
  records_fn = records_fns{rec_idx};
  fprintf('Loading records file %s\n',records_fn);
  records = load(records_fn,'relative_filename','relative_rec_num');
  [records_fn_dir,records_fn_name] = fileparts(records_fn);
  param.day_seg = records_fn_name(9:end);
  
  if strcmpi(param.post.img_type,'mat')
    image_dir = fullfile(ct_filename_out(param, ...
      param.post.out_dir, 'CSARP_qlook', true),param.day_seg);
  else
    image_dir = fullfile(ct_filename_out(param, ...
      param.post.out_dir, 'CSARP_post', true),'images',param.day_seg);
  end
  
  frames_fn = ct_filename_support(param,'','frames');
  
  if exist(frames_fn,'file')
    fprintf('Loading frames file %s\n', frames_fn);
    load(frames_fn);
    if ~isfield(frames,'nyquist_zone')
      frames.nyquist_zone = NaN*zeros(size(frames.frame_idxs));
    end
  else
    error('No frames files exists yet\n');
  end
  
  frm = 1;
  quit_cmd = false;
  while ~quit_cmd
    if frm > length(frames.frame_idxs)
      frm = length(frames.frame_idxs);
    end
    frm_id = sprintf('%s_%03i',param.day_seg,frm);
    
    if strcmpi(param.post.img_type,'mat')
      echo_fn = get_filename(image_dir,'Data_',frm_id,'.mat');
    else
      echo_fn = get_filename(image_dir,frm_id,'',['echo.' param.post.img_type]);
    end
    
    if isempty(echo_fn)
      warning('Missing frame %s', frm_id);
      % Go to next frame
      frm = frm + 1;
      continue
    end
    
    fprintf('  Loading output %s\n', echo_fn);
    if strcmpi(param.post.img_type,'mat')
      echo = load(echo_fn);
      set(h_image,'XData',1:size(echo.Data,2));
      set(h_image,'YData',echo.Time);
      set(h_image,'CData',lp(echo.Data));
      caxis(h_axes,[min(lp(echo.Data(:))) max(lp(echo.Data(:)))]);
      xlim(h_axes,[1 size(echo.Data,2)]);
      set(h_plot,'XData',1:size(echo.Data,2));
      set(h_plot,'YData',echo.Surface);
    else
      A = imread(echo_fn);
      set(h_image,'CData',A);
    end
    %     set(h_image,'XData',1:size(A,2));
    %     set(h_image,'YData',1:size(A,1));
    axis tight;
    zoom reset;
    if strcmpi(param.post.img_type,'mat')
      if isnan(frames.nyquist_zone(frm))
        val = input(sprintf('F(%03i:%04i:%01iD): ', frm, ...
          frames.proc_mode(frm), echo.param_get_heights.radar.wfs(1).nyquist_zone),'s');
      else
        val = input(sprintf('F(%03i:%04i:%01i): ', frm, ...
          frames.proc_mode(frm), frames.nyquist_zone(frm)),'s');
      end
    else
      if isnan(frames.nyquist_zone(frm))
        val = input(sprintf('F(%03i:%04i:-D): ', frm, ...
          frames.proc_mode(frm)),'s');
      else
        val = input(sprintf('F(%03i:%04i:%01i): ', frm, ...
          frames.proc_mode(frm), frames.nyquist_zone(frm)),'s');
      end
    end
    if isempty(val)
      % Go to next frame
      frm = frm + 1;
    else
      num_val = str2double(val);
      if ~isempty(strfind(upper(val),'J'))
        val = input('    Enter frame to jump to: ');
        if isempty(val)
          % Ignore command
        else
          frm = val;
        end
      elseif ~isnan(num_val)
        fprintf('    Changing frame type to %d\n', num_val);
        frames.proc_mode(frm) = num_val;
        frm = frm + 1;
      elseif strncmpi(val,'n',1)
        num_val = str2double(val(2:end));
        fprintf('    Changing nyquist zone to %.0f\n', num_val);
        frames.nyquist_zone(frm) = num_val;
        if val(1) == 'N'
          frm = frm + 1;
        end
      elseif strcmpi(val,'p')
        frm = frm - 1;
      elseif strcmpi(val,'c')
        sys_cmd = input('Enter system command: ','s');
        try
          eval(sys_cmd);
        end
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
      elseif strcmpi(val,'?')
        fprintf(' <enter>: go to next frame\n');
        fprintf(' #: set frame type to #\n');
        fprintf(' p: previous frame\n');
        fprintf(' n#: set nyquist zone to # (e.g. "n3")\n');
        fprintf(' j: jump to frame\n');
        fprintf(' c: enter matlab command\n');
        fprintf(' s: save frames file\n');
        fprintf(' q: quit\n');
      end
    end
  end
  
  val = input('Save before quitting (Y/N)? ','s');
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
