function calArray = takeNAcal(obj,Sparam_flag)
% calArray = takeNAcal(obj,Sparam_flag)
%
% obj is gpib object already opened
%
% Typical usage:
% getNAdata(1,1);
% putNAsettings(NA)
% putNAcal(NA.calArray)

% 8722C

fprintf(obj,'FORM5;');

if Sparam_flag == 0
  total_cal_params = 12;
else
  total_cal_params = 3;
end

for ind = 1:total_cal_params
    sprintf('OUTPCALC%02.0f;',ind);
    fprintf(obj,sprintf('OUTPCALC%02.0f;',ind));
    char(fread(obj,2,'char'));
    len = fread(obj,1,'uint16')/8;
    if len == 0
      fprintf('Calibration not available\n');
      calArray = [];
      return;
    end
    tmp = fread(obj,2*len,'float32');
    tmp = reshape(tmp,[2 len]).';
    calArray(:,ind) = tmp(:,1) + j*tmp(:,2);
end

return;
