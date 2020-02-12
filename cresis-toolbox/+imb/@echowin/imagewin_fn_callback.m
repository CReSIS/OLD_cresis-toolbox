function [fn,comment] = imagewin_fn_callback(obj,fn)
% [fn,comment] = echowin.imagewin_fn_callback(obj,fn)
%
% Support function for imb.echowin class. When the save or open buttons
% are used on the image params (imagewin class), this function is called
% to give the default filename.
%

fn_dir = fileparts(fn);

fn = fullfile(fn_dir, sprintf('Data_%s_%03d_%03d.png', ...
  obj.eg.cur_sel.day_seg, obj.eg.frms(1), obj.eg.frms(end)));
comment = '';

return;
