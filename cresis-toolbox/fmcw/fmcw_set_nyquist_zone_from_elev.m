function nz = fmcw_set_nyquist_zone_from_elev(param,agl_twtt,max_nz)
% nz = fmcw_set_nyquist_zone_from_elev(param,agl_twtt,max_nz)
%
% Prints the recommended default NZ and frames which deviate from this default
%
% param = parameter structure from parameter spreadsheet
%   radar parameters are used to determine nyquist zone boundaries
% agl_twtt = vector of two way travel times to the target
%   one twtt per record in the segment
% max_nz = maximum nyquist zone (cap everything to this)
%
% Called from fmcw_set_nyquist_zone_from_atm, find_gps_offset_using_surface

BW = abs(param.radar.wfs.f1-param.radar.wfs.f0) * param.radar.wfs.fmult;
alpha = BW / param.radar.wfs.Tpd;
nz_twtt = param.radar.fs/2/alpha;

nz = floor((agl_twtt+param.radar.wfs.Tsys) / nz_twtt);
interp_nz = isnan(nz);
nz = round(interp_finite(nz,NaN));

default_nz = mode(nz(nz <= max_nz & ~isnan(nz)));
% fprintf('Segment: %s\n', param.day_seg);
% fprintf('Nz twtt multiple (ns): %.3f\n', nz_twtt*1e9);
fprintf('%s\tRecommended default nz:\t%d\t%d\n', param.day_seg, default_nz, param.radar.wfs.nyquist_zone);
default_nz = param.radar.wfs.nyquist_zone;

nz(isnan(nz)) = default_nz;
nz(nz>max_nz) = max_nz;
nz(nz<0) = 0;

if 0
  figure(5); clf;
  plot(nz);
  hold on
  plot(find(~interp_nz), nz(~interp_nz), 'r.');
  title('Nyquist zone');
end

records_fn = ct_filename_support(param,'','records');
records = load(records_fn);
records.settings.nyquist_zone = nz;
save(records_fn,'-append','-struct','records','settings');
create_records_aux_files(records_fn,false);

spec_nz_frms = cell(max_nz+1,1);
bad_frms = [];

frame_changes = false;
frames_fn = ct_filename_support(param,'','frames');
load(frames_fn);
for frm = 1:length(frames.frame_idxs)
  if frm == length(frames.frame_idxs)
    stop_rec = length(nz);
  else
    stop_rec = frames.frame_idxs(frm+1)-1;
  end
  all_interp = all(interp_nz(frames.frame_idxs(frm):stop_rec));
  frame_nzs = unique(nz(frames.frame_idxs(frm):stop_rec));
  if 0
    % Debug and print out nyquist zones for each frame
    fprintf('Frame %d: ', frm);
    fprintf('%d ', frame_nzs);
    if some_interp
      fprintf('*');
    end
    fprintf('\n');
  end
  
  if numel(frame_nzs) == 1 && ~all_interp && frame_nzs <= max_nz
    spec_nz_frms{frame_nzs+1}(end+1) = frm;
    if frame_nzs ~= frames.nyquist_zone(frm)
      %fprintf('  %d: %d\n', frm, frame_nzs);
      frame_changes = true;
      frames.nyquist_zone(frm) = frame_nzs;
    end
  elseif ~all_interp && all(frame_nzs > max_nz)
    bad_frms(end+1) = frm;
    proc_mode = mod(floor(frames.proc_mode(bad_frms)/10),10);
    if frames.proc_mode(frm) == 0
      %fprintf('  %d: BAD\n', frm);
      frame_changes = true;
      frames.proc_mode(frm) = frames.proc_mode(frm) + 10;
    end
    frames.nyquist_zone(frm) = NaN;
  end
end
% for frame_nzs = 1:min(length(spec_nz_frms),max_nz+1)
%   if frame_nzs ~= default_nz
%     fprintf('  NZ %d: ', frame_nzs-1);
%     fprintf('%d ', spec_nz_frms{frame_nzs});
%     fprintf('\n');
%   end
% end
%
% fprintf('  BAD: ');
% fprintf('%d ', bad_frms);
% fprintf('\n');

% Update frames file
if frame_changes
  save_frames = true;
  if save_frames
    save(frames_fn,'frames');
  end
end

return;
