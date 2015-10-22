function imagewin_save_mat_callback(obj,h_obj,event)
% echowin.imagewin_save_mat_callback(obj,h_obj,event)
%
% Support function for imb.echowin class. When the save or open buttons
% are used on the image params (imagewin class), this function is called
% to give the default filename.
%

if any(size(obj.left_panel.imagewin.CData) ~= size(obj.eg.data))
  error('Save MAT only works with TWTT and Range Line axes selections');
end

% Ask for user input for which CSARP_ directory to save to
if 1
  prompt = {'Enter source to save in (E.g. CSARP_imageproc):'};
  dlg_title = 'Input new source';
  num_lines = 1;
  def = {'CSARP_imageproc'};
  answer = inputdlg(prompt,dlg_title,num_lines,def);
  answer = answer{1};
  
  if isempty(answer)
    return;
  end
end

% Get which image source (i.e. the II in "_img_II" part of the filename
% where img == 0 means no "_img_II" in the filename) is being viewed...
% save to that image
sourceMenus = get(obj.left_panel.sourceCM,'Children');
img = length(sourceMenus)-strmatch('on',get(sourceMenus,'Checked'))-3;
if img == 0
  fn_img_str = '';
else
  fn_img_str = sprintf('img_%02d_', img);
end

for cur_frame = obj.eg.frame_idxs
  fn = fullfile(ct_filename_out(obj.eg.cur_sel,'',answer,0), ...
    sprintf('Data_%s%s.mat',fn_img_str,obj.eg.frame_names{cur_frame}));
  fprintf('Saving %s\n', fn);
  
  mdata = [];
  mdata.Time = obj.eg.time;
  
  % Remove any data in echogram that is not part of the frame being loaded
  %   (e.g. some echograms were created with some data from neighboring
  %   frames and this must be removed)
  valid_mask = obj.eg.gps_time >= obj.eg.start_gps_time(cur_frame) ...
    & obj.eg.gps_time < obj.eg.stop_gps_time(cur_frame);
  
  mdata.Data = obj.left_panel.imagewin.CData(:,valid_mask);
  mdata.Latitude = obj.eg.latitude(:,valid_mask);
  mdata.Longitude = obj.eg.longitude(:,valid_mask);
  mdata.Elevation = obj.eg.elevation(:,valid_mask);
  mdata.GPS_time = obj.eg.gps_time(:,valid_mask);
  mdata.Surface = obj.eg.surface(:,valid_mask);
  
  % Save data to file
  fn_dir = fileparts(fn);
  if ~exist(fn_dir,'dir')
    mkdir(fn_dir);
  end
  save(fn,'-v6','-struct','mdata')
end

return;
