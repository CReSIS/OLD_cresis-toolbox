function getNAdata(numAve,numMeas,Sparam_flag)
% getNAdata(numAve,numMeas,Sparam_flag)
%
% numAve = hardware averages in NA
% numMeas = separate measurements stored to file
% Sparam_flag = 0 for full 2-port measurement (default)
%   1 for S11 measurement
%   2 for S22 measurement
%
% Typical usage:
% putNAsettings(NA)
% putNAcal(NA.calArray)
% getNAdata(1,1);
%
% NOTE: If the program crashes, it will not close the gpib
% object. The GPIB object must be closed before running the
% program again.  The GPIB object is "obj" and is set to be
% a global variable in this program. So from command line:
%   global obj; fclose(obj);

if ~exist('Sparam_flag','var')
  Sparam_flag = 0;
end

format compact; format long;
% Hewlett Packard (Agilent) 8722C, 8722D, and 8753D all appear to work

global obj;
% Open 'ni' National Instrument card, board 0, GPIB interface 16
%   Verify that your NA is set to address 16 (under "Local" button on NA)
obj = gpib('ni',0,16);
% Need to increase input buffer size:
%   Maximum number of points = 1601 points
%   Real + Imag = 2 doubles per point
%   8 bytes per double
set(obj,'InputBufferSize',2*1601*8);
% Increase time out to allow the instrument a long time to take data
%   (for really slow datasets, you may have to increase this more)
set(obj,'Timeout',100);

try
  fopen(obj);

  fprintf(obj,'CORR ON');
  fprintf(obj,'AVEROON');
  fprintf(obj,sprintf('AVERFACT %.0d',numAve));
  % The following settings should be set by the user already
  % (left here only for future reference)
  % fprintf(obj,'POIN 1601');
  % fprintf(obj,'STAR 300e6');
  % fprintf(obj,'STOP 2e9');
  % fprintf(obj,'IFBW 1000');
  % fprintf(obj,'POWE -10');

  for ind = 1:numMeas

    fprintf('Measurement %.0d of %.0d\n',ind,numMeas);
    fprintf(obj,sprintf('NUMG %.0d',numAve));

    if Sparam_flag == 0 || Sparam_flag == 1
      fprintf(obj,'S11;');
      fprintf(obj,'FORM5;');
      fprintf(obj,'OUTPDATA;');
      hdr = char(fread(obj,2,'char')).';
      if ~strcmp(hdr,'#A')
        fprintf('Error in header\n');
      end
      len = fread(obj,1,'ushort')/8;
      tmp = fread(obj,2*len,'float');
      tmp = reshape(tmp,[2 len]).';
      S11(:,ind) = tmp(:,1) + j*tmp(:,2);
    end

    if Sparam_flag == 0
      fprintf(obj,'S21;');
      fprintf(obj,'FORM5;');
      fprintf(obj,'OUTPDATA;');
      char(fread(obj,2,'char'));
      len = fread(obj,1,'ushort')/8;
      tmp = fread(obj,2*len,'float');
      tmp = reshape(tmp,[2 len]).';
      S21(:,ind) = tmp(:,1) + j*tmp(:,2);
      fprintf(obj,'S12;');
      fprintf(obj,'FORM5;');
      fprintf(obj,'OUTPDATA;');
      char(fread(obj,2,'char'));
      len = fread(obj,1,'ushort')/8;
      tmp = fread(obj,2*len,'float');
      tmp = reshape(tmp,[2 len]).';
      S12(:,ind) = tmp(:,1) + j*tmp(:,2);
    end

    if Sparam_flag == 0 || Sparam_flag == 2
      fprintf(obj,'S22;');
      fprintf(obj,'FORM5;');
      fprintf(obj,'OUTPDATA;');
      char(fread(obj,2,'char'));
      len = fread(obj,1,'ushort')/8;
      tmp = fread(obj,2*len,'float');
      tmp = reshape(tmp,[2 len]).';
      S22(:,ind) = tmp(:,1) + j*tmp(:,2);
    end
  end

  NA = getNASettings(obj);
  NA.ave = numAve;
  NA.calArray = getNAcal(obj,Sparam_flag);
  fprintf(obj,'*IDN?');
  NA.idn = fscanf(obj)
  fclose(obj);
catch ME
  fclose(obj);
  rethrow(ME);
end

filename = input('Input filename: ','s');
switch Sparam_flag
  case 0
    save(filename,'S11','S21','S12','S22','NA');
  case 1
    save(filename,'S11','NA');
  case 2
    save(filename,'S22','NA');
end

return;
