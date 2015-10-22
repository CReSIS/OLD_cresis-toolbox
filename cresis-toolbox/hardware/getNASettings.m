function NA = getNASettings(obj)
% obj is gpib object already opened

format long; format compact;

% 8722C

fprintf(obj,'SWET?;');
NA.sweepTime = fscanf(obj,'%f');
fprintf(obj,'POIN?;');
NA.numOfPoints = fscanf(obj,'%f');
fprintf(obj,'STAR?;');
NA.startFreq = fscanf(obj,'%f');
fprintf(obj,'STOP?;');
NA.stopFreq = fscanf(obj,'%f');
fprintf(obj,'IFBW?;');
NA.ifBW= fscanf(obj,'%f');
fprintf(obj,'POWE?;');
NA.power = fscanf(obj,'%f');

return;
