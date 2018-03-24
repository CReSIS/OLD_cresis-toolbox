function [] = records_line(dateorstruct,seg_id,file_adcs,geotiff_fn,...
    frame_mode,file_version)

%Takes the inputs and prints a line of the records worksheet
    %date needs to be a string following format YYYYMMDD
    %seg_id should be an integer
    %file.adcs A row vector which specified which ADCs will be used in this
        %segment and the order that they will be referenced. 
        %Usually this is all of the ADCs
    %gps_en Include GPS data in the records file 
        %generally always true unless you are collecting lab data
    %geotiff_fn This is the background geotiff that is used 
        %by the GUI during frames creation
    %frame_mode: Decimal scalar interpreted as a bit mask. 
        %For example to set bit 0 and bit 2, the field would be set to 1+4=5
        %bit 0: [0] Use create_frames.m GUI or [1] use autogenerate_frames.m
        %bit 1: [1] Enable special mode for snow2/kuband2 Antarctica 
            %2012 DC8 only which uses the settings field to create frames 
            %if loopback mode or nyquist zone settings change.
            %No affect for other datasets.
        %bit 2: [1] Overwrite frame without asking
        %bit 3-7: Reserved/ignored
    %debug_level: Leave on 1
if nargin <2
    S = dateorstruct;
    date = S.date;
    seg_id = S.seg_id;
    if ~isfield(S,'file_adcs')
        file_adcs = '1';
    else
        file_adcs = S.file_adcs;
    end
    geotiff_fn = S.geotiff_fn;
    if ~isfield(S,'frame_mode')
        frame_mode = 1;
    else
        frame_mode = S.frame_mode;
    end
    file_version = S.file_version;
else
    date = dateorstruct;
end
fprintf('%s\t%.0f\t%s\t%s\t%.0f\t%.0f\n',...
    date,seg_id,file_adcs,geotiff_fn,frame_mode,file_version)
end

    