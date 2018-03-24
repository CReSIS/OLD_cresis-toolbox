function [] = vectors_line(dateorstruc, seg_id, adc, start_idx, ...
    stop_idx, base_dir, adc_folder_name, file_prefix, gps_time_offset)

%Takes the inputs and prints a line of the vectors worksheet
    %Date needs to be a string following format YYYYMMDD
    %seg_id should be an integer
    %radar_num is only necessary for mcords radars. Input 0 to leave blank
        %Don't include as field in structure if unused.
    %adc is typically just 1 or another number
    %start_idx should be an integer
    %stop_idx should be an integer larger than start_idx
    %base_dir should be a string ending with a /
    %adc_folder_name should be a string with no / at the beginning or end
    %file_prefix should be a string 
    %gps_time_offset is an floating point number
        %Don't include as field in structure if unused.
    %gps_utc_time_halved is a logical (i.e. 1 or 0)
        %Only kaband radars should have a 1
        %Don't include as field in structure if unused.
    
if nargin <2
    S = dateorstruc; %Pull from structure
    date = S.date;
    seg_id = S.seg_id;
    if ~isfield(S,'adc')
        adc = 1;%Typical value for this column
    else
        adc = S.adc;
    end
    start_idx = S.start_idx;
    stop_idx = S.stop_idx;
    base_dir = S.base_dir;
    adc_folder_name = S.adc_folder_name;
    file_prefix = S.file_prefix;
    if ~isfield(S,'gps_time_offset')
        gps_time_offset = 0; %assumed 0
    else
        gps_time_offset = S.gps_time_offset;
    end
else
    date = dateorstruc;
    
end


fprintf('%s\t%02d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\n', date, seg_id, ...
        adc, start_idx, stop_idx, base_dir, adc_folder_name,file_prefix,...
        gps_time_offset);
end