function [success] = update_layerdata_format_task(param)
% [success] = update_layerdate_format_task(param)
%
% Cluster task for updateing layer data format. 
%
% param = struct controlling the updates.
%  .load = structure for which records to load
%   .records_fn = filename of records file
%   .recs = current records
%   .imgs = cell vector of images to load, each image is Nx2 array of
%     wf/adc pairs
%     NOTE: wfs/adc pairs are not indices into anything, they are the absolute
%     waveform/adc numbers. The records file will be loaded to decode
%     which index each wf/adc belongs to.
%  .debug_level = debug level (scalar integer)
%
%  .proc = structure containing information about framing
%   .frm = only used to determine the filename of the output
%
% success = boolean which is true when the function executes properly
%   if a task fails before it can return success, then success will be
%   empty
%
% Author: Jilu Li
%
% See also update_layerdata_format.m


%% Load the records
% =========================================================================
records_fn = ct_filename_support(param,'','records');
records = read_records_aux_files(records_fn,param.load.recs);

  

%% Load the old layer data file of the frame
frm = param.load.frm;
layer_fn = fullfile(ct_filename_out(param,param.layerdata_source,''), ...
  sprintf('Data_%s_%03d.mat', param.day_seg, frm));
lay = load(layer_fn);

% Remove data that is not contained within frame boundaries (old version
% layer data file contains overlapped data at boundaries)
if param.frame_overlap_removal
  frm_mask = logical(zeros(size(lay.GPS_time)));
  frm_mask(lay.GPS_time >= records.gps_time(1) & lay.GPS_time <= records.gps_time(end)) = true;
else
  frm_mask = logical(ones(size(lay.GPS_time)));
end

layer.GPS_time = lay.GPS_time(frm_mask);
layer.Pitch = interp1(records.gps_time,records.pitch,layer.GPS_time);
layer.Roll = interp1(records.gps_time,records.roll,layer.GPS_time);
layer.Heading = interp1(records.gps_time,records.heading,layer.GPS_time);
layer.file_version = param.file_version_new;

if strcmp(param.file_version_old,'0')
  for lyr_idx = 1:length(lay.layerData)
    layer.id(lyr_idx,1) = lyr_idx;
    layer.twtt(lyr_idx,:) = lay.layerData{lyr_idx}.value{2}.data(frm_mask);
    layer.quality(lyr_idx,:) = uint8(lay.layerData{lyr_idx}.quality(frm_mask));
    layer.type(lyr_idx,:) = uint32(ones(size(lay.layerData{lyr_idx}.value{2}.data(frm_mask))));
    idxs_manual = find(isfinite(lay.layerData{lyr_idx}.value{1}.data(frm_mask)));
    if ~isempty(idxs_manual)
      layer.type(lyr_idx,idxs_manual) = uint32(0);
    end
  end
else
end

%% Save the new layer data
% =========================================================================
layer_fn = fullfile(ct_filename_out(param,param.layerdata_dest,''), ...
  sprintf('Data_%s_%03d.mat', param.day_seg, frm));
fprintf('  Save %s\n', layer_fn);
ct_save(layer_fn,'-v7.3', '-struct', 'layer');

%% Done
% =========================================================================

fprintf('%s done %s\n', mfilename, datestr(now));

success = true;
