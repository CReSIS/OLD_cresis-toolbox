% script basic_remove_mcords_digital_errors
%
% Attempts to detect and remove (by setting to zero) digital errors
% The specific digital error that it looks for is a sequence of 3 or 4
% very large values.
%
% Requires that data be stored in a variable called "data" which will
% be overwritten. Data may be a cell array of 2-D/3-D arrays or just
% a single 2-D/3-D array.
%
% Author: John Paden

if ~exist('tstart','var')
  fprintf('Removing mcords digital errors\n');
else
  fprintf('Removing mcords digital errors (%.1f sec)\n', toc(tstart));
end
  
if iscell(data)
  for wf_idx = 1:length(data)
    data{wf_idx} = permute(data{wf_idx},[3 1 2]);
    for adc_idx = 1:size(data{wf_idx},1)
      data_pow = abs(data{wf_idx}(adc_idx,:,:)).^2;
      
      data_pow = permute(data_pow,[2 3 1]);
      cfar_threshold = medfilt2(data_pow,[5 3]);
      cfar_threshold(:,1:3) = repmat(cfar_threshold(:,5),[1 3]);
      cfar_threshold(:,end-2:end) = repmat(cfar_threshold(:,end-4),[1 3]);
      cfar_threshold(1:5,:) = repmat(cfar_threshold(7,:),[5 1]);
      cfar_threshold(end-4:end,:) = repmat(cfar_threshold(end-6,:),[5 1]);
      data_pow = permute(data_pow,[3 1 2]);
      cfar_threshold = permute(cfar_threshold,[3 1 2]);
      
      data{wf_idx}(adc_idx,data_pow > cfar_threshold*1000) = 0;
    end
    data{wf_idx} = permute(data{wf_idx},[2 3 1]);
  end
else
  data = permute(data,[3 1 2]);
  for adc_idx = 1:size(data,1)
    data_pow = abs(data(adc_idx,:,:)).^2;

    data_pow = permute(data_pow,[2 3 1]);
    cfar_threshold = medfilt2(data_pow,[5 3]);
    cfar_threshold(:,1:3) = repmat(cfar_threshold(:,5),[1 3]);
    cfar_threshold(:,end-2:end) = repmat(cfar_threshold(:,end-4),[1 3]);
    cfar_threshold(1:5,:) = repmat(cfar_threshold(7,:),[5 1]);
    cfar_threshold(end-4:end,:) = repmat(cfar_threshold(end-6,:),[5 1]);
    data_pow = permute(data_pow,[3 1 2]);
    cfar_threshold = permute(cfar_threshold,[3 1 2]);

    data(adc_idx,data_pow > cfar_threshold*1000) = 0;
  end
  data = permute(data,[2 3 1]);
end

return;
