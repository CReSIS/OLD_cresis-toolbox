function putNAcal(calArray)
% putNAcal(calArray)
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
set(obj,'OutputBufferSize',2*1601*8);
set(obj,'Timeout',100);

fopen(obj);

fprintf(obj,'*IDN?');
idn = fscanf(obj)

fprintf(obj,'FORM5;');
fprintf(obj,'CALIFUL2;');
for ind = 1:12
    set(obj,'EOIMode','off');
    sprintf('INPUCALC%02.0f;',ind)
    fprintf(obj,sprintf('INPUCALC%02.0f;',ind));
    fwrite(obj,'#A','schar')
    fwrite(obj,size(calArray,1)*8,'uint16')
    outData = reshape([real(calArray(:,ind))'; imag(calArray(:,ind))'],[size(calArray,1)*2 1]);
    set(obj,'EOIMode','on');
    fwrite(obj,outData,'float32')
end
fprintf(obj,'SAVC;');

fclose(obj);

return;
