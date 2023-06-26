function save_callback(obj,h_obj,event)
% imagewin.save_callback(obj,h_obj,event)
%
% Callback when save button is pushed

if ~isempty(obj.img_fn_save_fh)
  [obj.img_fn,comment] = obj.img_fn_save_fh(obj.img_fn);
else
  comment = '';
end

% Save image data to a file
[fn, fn_dir] = uiputfile( ...
  obj.img_fn_save_formats, ...
  'Save as', obj.img_fn);

if isempty(fn) || fn(1) == 0
  return;
end

obj.img_fn = fullfile(fn_dir,fn);

if ~isempty(obj.save_fh)
  status = obj.save_fh();
else
  status = 1;
end
if status <= 0
  return;
end

[img_fn_dir img_fn_name img_fn_ext] = fileparts(obj.img_fn);

if strcmpi(img_fn_ext,'.png')
  fprintf('Saving %s\n', obj.img_fn);
  
  % Convert image to 16 bit unsigned
  CData = get(obj.img,'CData');
  img_min = min(CData(:));
  img_max = max(CData(:));
  CData = CData - img_min;
  img_data = CData ./ (img_max-img_min) * 65535;
  img_data = uint16(65535 - img_data);
  
  % Store metadata in PNG comments
  start_row = 1;
  start_col = 1;
  comment = sprintf('img_min = %.12f; img_max = %.12f; Row = %d; Col = %d; %s', ...
    img_min, img_max, start_row, start_col, comment);
  
  % Write PNG file
  imwrite(img_data, obj.img_fn,'Comment',comment);
  
else
  error('Invalid file format %s detected in filename %s', img_fn_ext, obj.img_fn);
  
end

end
