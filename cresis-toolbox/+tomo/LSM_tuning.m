function [matrix_x,matrix_y,flag]= LSM_tuning (param, options, data)
% function LSM_tuning (params, options)
% gets matrix_x and matrix_y and performs interpolation if possible
% See also: run_LSM_tuning.m
%
global gRadar;
clear('param_override');
data_fn_dir = ct_filename_out(param, options.name, '');
labels.success = 1;
fprintf('\nLSM: Running frame %s_%03d (%s)\n',param.day_seg, data.frm, datestr(now,'HH:MM:SS'));
lsm_tic = tic;

data_fn_name = sprintf('Data_%s_%03d.mat',param.day_seg,data.frm);
data_fn      = fullfile(data_fn_dir, data_fn_name);
imds         = imageDatastore(data_fn, 'FileExtensions', '.mat');
obj          = tomo.LSMObject_tuning(imds.Files);

obj.setLSMOptions('y',options.lsm.y,'dy',options.lsm.dy,'outerIter',400);

[flag, Labels.top, Labels.bot, matrix_x, matrix_y] = obj.runLSM();

temp = size(matrix_y,3);

for i =1:2
  for idx = 1:temp
    try
      matrix_y(i,:,idx) = interp1(matrix_x(i,:,idx),matrix_y(i,:,idx),1:length(data.Bottom),'linear','extrap');
    catch ME
      flag(idx) = 0;
      continue
    end
  end
end

end