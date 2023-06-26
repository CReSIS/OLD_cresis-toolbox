function [flags,info,eof_error] = TDMS_processLeadIn(fid,lastLetter)
%TDMS_processLeadIn  Reads the "Lead In" portion of the data segment
%
%   [flags,info] = TDMS_processLeadIn(fid)
%
%   OUTPUTS
%   ===============================
%   flags (struct)  :
%           .hasMetaData    - true if segment has meta data
%           .kTocNewObjList - true if there is a new object order for raw data
%           .hasRawData     - true if segment has raw data
%           .isInterleaved  - true if raw data is interleaved
%           .isBigEndian    - true if bigEndian
%           .hasRawDaqMX    - true if DaqMX raw data is present
%   info (struct)   :
%           .verNumber   - ill defined, current versions should have no
%                          impact on the ability to read the file according
%                          to current specifications (Jan 2010), i.e. no
%                          processing change is needed if the version is one 
%                          or the other
%           .segLength   - length to the next segment in bytes (after lead in)
%           .metaLength  - length of meta data in bytes

    %LEAD IN PROCESSING
    %================================================
    %1) TDSm - (indicats lead in)
    Ttag = fread(fid,1,'uint8');
    Dtag = fread(fid,1,'uint8');
    Stag = fread(fid,1,'uint8');
    mtag = fread(fid,1,'uint8');

    if ~(Ttag == 84 && Dtag == 68 && Stag == 83 && mtag == lastLetter)
        error('Unexpected lead in header')
    else
        %2) Kind of Data
        tocMask = fread(fid,1,'uint32');    
        flags = struct(...
            'hasMetaData',      bitget(tocMask,2),...
            'kTocNewObjList',   bitget(tocMask,3),...
            'hasRawData',       bitget(tocMask,4),...
            'isInterleaved',    bitget(tocMask,6),... %false, contiguous
            'isBigEndian',      bitget(tocMask,7),... %false, little-endian
            'hasRawDaqMX',      bitget(tocMask,8));
        
        if flags.isBigEndian
            error('Currently code is unable to handle Big-Endian format')
        end
        
        if flags.hasRawDaqMX
            error('Currently code is unable to ignore/handle Raw Daq MX data')
        end

        %3) Version Number
        info = struct(...
            'verNumber',  fread(fid,1,'uint32'),...  %This is not very well defined
            'segLength',  fread(fid,1,'uint64'),...
            'metaLength', fread(fid,1,'uint64'));  
    
    %NOT THROWING THIS ERROR RIGHT NOW
    %======================================================================
    %This should really quit instead of throwing an error
    eof_error = info.segLength == 2^64-1;
%         if info.segLength == 2^64-1
%             error('File got corrupted when saving')
%         end
    
        %This most likely suggests an error in the reading
        if ~eof_error && flags.hasMetaData ~= (info.metaLength ~= 0)
            error('Flags suggest presence of meta data but no meta is present according to length value')
        end
    end