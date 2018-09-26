function nz = fmcw_set_nyquist_zone_from_elev(param,agl_gps_time,agl_twtt,max_nz)
% nz = fmcw_set_nyquist_zone_from_elev(param,agl_gps_time,agl_twtt,max_nz)
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

%% Setup
wf = 1;
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);
param.load.imgs = {[1 1]};
[wfs,~] = data_load_wfs(param,records);

%% Calculate Nyquist zone based on AGL
BW = abs(wfs(wf).f1-wfs(wf).f0) * wfs(wf).fmult;
alpha = BW / wfs(wf).Tpd;
nz_twtt = param.radar.fs/2/alpha;

nz = floor((agl_twtt+wfs(wf).Tsys-wfs(wf).t_ref) / nz_twtt);
interp_nz = isnan(nz);
nz = round(interp_finite(nz,NaN));

default_nz = mode(nz(nz <= max_nz & ~isnan(nz)));

nz(isnan(nz)) = default_nz;
nz(nz<0) = 0;
nz(nz>max_nz) = max_nz;

if 0
  fprintf('Nz twtt multiple (ns): %.3f\n', nz_twtt*1e9);
  
  figure(5); clf;
  plot(nz);
  hold on
  plot(find(~interp_nz), nz(~interp_nz), 'r.');
  title('Nyquist zone');
end

records.settings.nyquist_zone = interp1(agl_gps_time,nz,records.gps_time,'nearest','extrap');
save(records_fn,'-append','-struct','records','settings');
create_records_aux_files(records_fn,false);
