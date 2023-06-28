function [Data] = load_data_mcrds(filename,Header,start_record,stop_record,channel,waveform)

Data = [];

% OPEN FILE     ==========================================================
fid = fopen( filename ,'r');
status = fseek(fid,Header.IndexData,'bof');

% SKIP PAST UNWANTED RECORDS    =========================================
status = fseek(fid,(start_record - 1)*[4+4+4+8+4+2*Header.IndexRecordStop(end,end)],'cof');

% START READING RECORDS     ==============================================
temp = zeros(Header.IndexRecordStop(end,end),(stop_record - start_record + 1));
recordindex = 1;
while recordindex <= (stop_record - start_record + 1)  
  
    % Skip past datatype,time and DAQ error bytes
    status = fseek(fid,4+4+4+8+4,'cof');
    
    temp(:,recordindex) = fread(fid,Header.IndexRecordStop(end,end),'uint16');
    
    recordindex = recordindex + 1;
    
end

fclose(fid);

Data = temp(Header.IndexRecordStart(channel,waveform):Header.IndexRecordStop(channel,waveform),:);

% data = data-mean(mean(data));
% 
% if waveform == 1
%     Data = [zeros(1,stop_record-start_record+1);data];
% else
%     Data = data;
% end


