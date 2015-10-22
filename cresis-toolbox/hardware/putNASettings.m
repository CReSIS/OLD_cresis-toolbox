function putNASettings(NA)
% obj is gpib object already opened
%
% Typical usage:
% getNAdata(1,1);
% putNAsettings(NA)
% putNAcal(NA.calArray)

format long; format compact;

% 8722C

global obj;
obj = gpib('ni',0,16);
set(obj,'InputBufferSize',2*1601*8);
set(obj,'Timeout',100);

fopen(obj);
fprintf(obj,sprintf('POIN %.0d;',NA.numOfPoints));
fprintf(obj,sprintf('STAR %.12d;',NA.startFreq));
fprintf(obj,sprintf('STOP %.12d;',NA.stopFreq));
fprintf(obj,sprintf('IFBW %.0d;',NA.ifBW));
fprintf(obj,sprintf('POWE %.2d;',NA.power));
fprintf(obj,sprintf('SWET %.12d;'));

fclose(obj);

return;
