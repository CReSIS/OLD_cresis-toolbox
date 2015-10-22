function open_callback(obj,h_obj,event)
% imagewin.open_callback(obj,h_obj,event)
%
% Callback when open button is pushed

if ~isempty(obj.img_fn_open_fh)
  [obj.img_fn,comment] = obj.img_fn_open_fh(obj.img_fn);
else
  comment = '';
end

[fn, fn_dir, filterindex] = uigetfile( ...
  obj.img_fn_open_formats, ...
  'Pick a file', ...
  obj.img_fn, 'MultiSelect', 'on');

if isempty(fn) || fn(1) == 0
  return;
end

obj.img_fn = fullfile(fn_dir,fn);

if ~isempty(obj.open_fh)
  status = obj.open_fh();
else
  status = 1;
end
if status <= 0
  return;
end

[img_fn_dir img_fn_name img_fn_ext] = fileparts(obj.img_fn);

if any(strcmpi(img_fn_ext,{'.png','.jpg','.jpeg','.jp2','.jpf','.jpx','.j2c','.j2k','.tif','.tiff'}))
  fprintf('Importing %s', obj.img_fn);
  img_data = imread(obj.img_fn);
  img_info = imfinfo(obj.img_fn);
  fprintf('  depth (%d) class (%s)\n', size(img_data,3), class(img_data));

  try
    eval(img_info.Comment);
  catch ME
    warning('Evaluation of image comments (%s) failed with error:\n %s', ...
      img_info.Comment, ME.getReport);
  end

  % Write variables in case img_info.Comment failed
  if ~exist('img_min','var')
    img_min = 0;
  end
  if ~exist('img_max','var')
    img_max = 100;
  end
  if ~exist('start_row','var')
    start_row = 1;
  end
  if ~exist('start_col','var')
    start_col = 1;
  end
  
  % Get current image data
  CData = get(obj.img,'CData');
  
  % Check size of img_data and truncate if necessary
  img_data = img_data(1 : min(size(CData,1)-start_row+1, size(img_data,1)), ...
    1 : min(size(CData,2)-start_col+1, size(img_data,2)), :);
  
  % Convert image data from raw image format
  if size(img_data,3) ~= 1 && isa(img_data,'uint8')
    CData(start_row+(0:size(img_data,1)-1), start_col+(0:size(img_data,2)-1)) ...
      = 10.^((mean(double(255-img_data),3)/255 * (img_max-img_min) + img_min)/10);
  elseif size(img_data,3) == 1 && isa(img_data,'uint8')
    CData(start_row+(0:size(img_data,1)-1), start_col+(0:size(img_data,2)-1)) ...
      = 10.^((double(255-img_data)/255 * (img_max-img_min) + img_min)/10);
  elseif size(img_data,3) == 1 && isa(img_data,'uint16')
    CData(start_row+(0:size(img_data,1)-1), start_col+(0:size(img_data,2)-1)) ...
      = 10.^((double(65535-img_data)/65535 * (img_max-img_min) + img_min)/10);
  else
    error('Unsupported image type depth(%d) class(%s) for file %s\', ...
      size(img_data,3), class(img_data), obj.img_fn);
  end
  
  % Set current image data
  obj.set_cdata(CData);
  
else
  error('Invalid file format %s detected in filename %s', img_fn_ext, obj.img_fn);
  
end

end
