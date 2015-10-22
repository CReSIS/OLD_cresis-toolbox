function [TimeComputer TimeRadar] = load_time_acords(filename,Header)

TimeComputer = [];
TimeRadar = [];

% OPEN FILE     ==========================================================
fid = fopen( filename ,'r');
status = fseek(fid,Header.IndexData,'bof');

% START READING TIMES     ================================================
block_length = Header.IndexRecordStop(end,end);
recordindex = 1;
TimeComputer = [];
TimeRadar = [];
while recordindex <= Header.NumberRecords  

    datatype = fread(fid,1,'int32');      %read datatype
    
    seconds = fread(fid,1,'int32');
    useconds = fread(fid,1,'int32');
    TimeComputer = [TimeComputer seconds+(1e-6*useconds)];

    TimeRadar_temp = fread(fid,1,'int64');
    TimeRadar = [TimeRadar TimeRadar_temp/10e6];
   
    % Skip past DAQ error bytes and data
    status = fseek(fid,4+2*block_length,'cof');
      
    recordindex = recordindex + 1;
    
end

fclose(fid);


