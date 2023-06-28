% script deconv.write_deconv_wf_to_frames
%
% Useful for update_frames.m so you can view and keep track of suggestions
% for new deconvolution waveforms to use for various frames.

% base_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2009_Greenland_P3/CSARP_deconv/';
% base_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Greenland_DC8/CSARP_deconv/';
% base_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2011_Greenland_P3/CSARP_deconv/';
% base_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Greenland_P3/CSARP_deconv/';
% base_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2009_Antarctica_DC8/CSARP_deconv/';
% base_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2010_Antarctica_DC8/CSARP_deconv/';
% base_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2011_Antarctica_DC8/CSARP_deconv/';
% base_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Antarctica_DC8/CSARP_deconv/';
base_dir = '/cresis/snfs1/dataproducts/ct_data/snow/SEASON/CSARP_deconv/';

% params = read_param_xls(ct_filename_param('snow_param_2009_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_DC8.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2009_Antarctica_DC8.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2010_Antarctica_DC8.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2011_Antarctica_DC8.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2012_Antarctica_DC8.xls'));
params = read_param_xls(ct_filename_param('SPREADSHEET.xls'));

for param_idx = 1:length(params)
  param = params(param_idx);
  
  if ~iscell(param.cmd.generic) && param.cmd.generic
    frames_fn = ct_filename_support(param,'','frames');
    load(frames_fn);
    
    if 0
      if isfield(frames,'deconv_wf')
        frames = rmfield(frames,'deconv_wf');
        save(frames_fn,'-V7','frames');
      end
    else
      if ~isfield(frames,'deconv_wf')
        existed_before = false;
        frames.deconv_wf = NaN*zeros(size(frames.frame_idxs));
      else
        existed_before = true;
      end
      old_frames = frames;
      
      out_dir = fullfile(base_dir,param.day_seg);
      fprintf('%s\n', param.day_seg);
      fns = get_filenames(out_dir,'Data_','','.mat');
      fprintf('%10s\t%10s\n', 'frm', 'deconv wfs');
      for fn_idx = 1:length(fns)
        mdata = load(fns{fn_idx},'custom');
        if isfield(mdata,'custom')
          [~,fn_name] = fileparts(fns{fn_idx});
          frm = str2double(fn_name(end-2:end));
          
          fprintf('%10d',frm);
          unique_deconv = unique(mdata.custom.deconv_filter_idx);
          fprintf('\t%10d',unique_deconv);
          
          fprintf('\n');
          
          frames.deconv_wf(frm) = mode(mdata.custom.deconv_filter_idx);
        end
      end
      
      update_field = 'deconv_wf';
      diff_frms = find((frames.(update_field)~=old_frames.(update_field) & ~(isnan(frames.(update_field)) & isnan(old_frames.(update_field)))) ...
        | (frames.nyquist_zone~=old_frames.nyquist_zone & ~(isnan(frames.nyquist_zone) & isnan(old_frames.nyquist_zone))));
      
      if ~isempty(diff_frms)
        if 0&&existed_before
          fprintf('Frames changed from previous save, are you sure you want to overwrite?\n');
          keyboard
        end
        save(frames_fn,'-V7','frames');
      end
    end
  end
  
end
