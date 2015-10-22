function [dbfData, dbfFieldNames] = dbfread(filename, records2read, requestedFieldNames)
%DBFREAD Read the specified records and fields from a DBF file.
%
%   [DATA, NAMES] = DBFREAD(FILE) reads numeric, float, character and date
%   data and field names from a DBF file, FILE.
%
%   [DATA, NAMES] = DBFREAD(FILE, RECORDS2READ) reads only the record
%   numbers specified in RECORDS2READ, a scalar or vector.
%
%   [DATA, NAMES] = DBFREAD(FILE, RECORDS2READ, REQUESTEDFIELDNAMES) reads
%   the data from the fields, REQUESTEDFIELDNAMES, for the specified
%   records. REQUESTEDFIELDNAMES must be a cell array. The fields in the
%   output will follow the order given in REQUESTEDFIELDNAMES.
%
%   Examples:
%
%       % Get all records and a list of the field names from a DBF file.
%       [DATA,NAMES] = dbfread('c:\matlab\work\mydbf')
%
%       % Get data from records 3:5 and 10 from a DBF file.
%       DATA = dbfread('c:\matlab\work\mydbf',[3:5,10])
%
%       % Get data from records 1:10 for three of the fields in a DBF file.
%       DATA = dbfread('c:\matlab\work\mydbf',1:10,{'FIELD1' 'FIELD3' 'FIELD5'})
%
%   See also XLSREAD, DLMREAD, DLMWRITE, LOAD, FILEFORMATS, TEXTSCAN.

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.0 $   $Date: 2008/04/18 05:58:17 $

[pathstr,name,ext] = fileparts(filename);

dbfFileId = fopen(filename,'r','ieee-le');
if (dbfFileId == -1)
    dbfFileId = fopen(fullfile(pathstr, [name '.dbf']),'r','ieee-le');
end
if (dbfFileId == -1)
    dbfFileId = fopen(fullfile(pathstr, [name '.DBF']),'r','ieee-le');
end

if (dbfFileId == -1)
    eid = sprintf('MATLAB:%s:missingDBF', mfilename);
    msg = sprintf('Failed to open file %s.dbf or file %s.DBF.',...
            name, name);
    error(eid,'%s',msg)
end

info = dbfinfo(dbfFileId);
if ~exist('requestedFieldNames','var')
    dbfFieldNames = {info.FieldInfo.Name};
    requestedFieldNames = dbfFieldNames;
else
    dbfFieldNames = (info.FieldInfo(matchFieldNames(info,requestedFieldNames)).Name);
end
fields2read = matchFieldNames(info,requestedFieldNames);

% The first byte in each record is a deletion indicator
lengthOfDeletionIndicator = 1;

if ~exist('records2read','var')
    records2read = (1:info.NumRecords);
elseif max(records2read) > info.NumRecords
    eid = sprintf('MATLAB:%s:invalidRecordNumber', mfilename);
    msg = sprintf('Record number %d does not exist, please select from the range 1:%d.',...
        max(records2read), info.NumRecords);
    error(eid,'%s',msg)
end

% Loop over the requested fields, reading in the data
dbfData = cell(numel(records2read),numel(fields2read));
for k = 1:numel(fields2read),
    n = fields2read(k);
    fieldOffset = info.HeaderLength ...
                  + sum([info.FieldInfo(1:(n-1)).Length]) ...
                  + lengthOfDeletionIndicator;
    fseek(dbfFileId,fieldOffset,'bof');
    formatString = sprintf('%d*uint8=>char',info.FieldInfo(n).Length);
    skip = info.RecordLength - info.FieldInfo(n).Length;
    data = fread(dbfFileId,[info.FieldInfo(n).Length info.NumRecords],formatString,skip);
    dbfData(:,k) = feval(info.FieldInfo(n).ConvFunc,(data(:,records2read)'));
%     dbfData(:,k) = info.FieldInfo(n).ConvFunc(data(:,records2read)');
end

fclose(dbfFileId);

%--------------------------------------------------------------------------
function fields2read = matchFieldNames(info, requestedFieldNames)
% Determine which fields to read.

allFieldNames = {info.FieldInfo.Name};
if isempty(requestedFieldNames)
    if ~iscell(requestedFieldNames)
        % Default case: User omitted the parameter, return all fields.
        fields2read = 1:info.NumFields;
    else
        % User supplied '{}', skip all fields.
        fields2read = [];
    end
else
    % Match up field names to see which to return.
    fields2read = [];
    for k = 1:numel(requestedFieldNames)
        index = strmatch(requestedFieldNames{k},allFieldNames,'exact');
        if isempty(index)
            wid = sprintf('MATLAB:%s:nonexistentDBFName',mfilename);
            wrn = sprintf('DBF name ''%s'' %s\n%s',requestedFieldNames{k},...
                     'doesn''t match an existing DBF name.',...
                     '         It will be ignored.');
            warning(wid,wrn)
        end
        for l = 1:numel(index)
            % Take them all in case of duplicate names.
            fields2read(end+1) = index(l);
        end
    end
end

%--------------------------------------------------------------------------

function info = dbfinfo(fid)
%DBFINFO Read header information from DBF file.
%   FID File identifier for an open DBF file.
%   INFO is a structure with the following fields:
%      Filename       Char array containing the name of the file that was read
%      DBFVersion     Number specifying the file format version
%      FileModDate    A string containing the modification date of the file
%      NumRecords     A number specifying the number of records in the table
%      NumFields      A number specifying the number of fields in the table
%      FieldInfo      A 1-by-numFields structure array with fields:
%         Name        A string containing the field name 
%         Type        A string containing the field type 
%         ConvFunc    A function handle to convert from DBF to MATLAB type
%         Length      A number of bytes in the field
%      HeaderLength   A number specifying length of the file header in bytes
%      RecordLength   A number specifying length of each record in bytes

%   Copyright 1996-2005 The MathWorks, Inc.
%   $Revision: 1.1.10.4 $  $Date: 2005/11/15 01:07:13 $

[version, date, numRecords, headerLength, recordLength] = readFileInfo(fid);
fieldInfo = getFieldInfo(fid);

info.Filename     = fopen(fid);
info.DBFVersion   = version;
info.FileModDate  = date;
info.NumRecords   = numRecords;
info.NumFields    = length(fieldInfo);
info.FieldInfo    = fieldInfo;
info.HeaderLength = headerLength;
info.RecordLength = recordLength;

%----------------------------------------------------------------------------
function [version, date, numRecords, headerLength, recordLength] = readFileInfo(fid)
% Read from File Header.

fseek(fid,0,'bof');

version = fread(fid,1,'uint8');

year  = fread(fid,1,'uint8') + 1900;
month = fread(fid,1,'uint8');
day   = fread(fid,1,'uint8');

dateVector = datevec(sprintf('%d/%d/%d',month,day,year));
dateForm = 1;% dd-mmm-yyyy
date = datestr(dateVector,dateForm);

numRecords   = fread(fid,1,'uint32');
headerLength = fread(fid,1,'uint16');
recordLength = fread(fid,1,'uint16');

%----------------------------------------------------------------------------
function fieldInfo = getFieldInfo(fid)
% Form FieldInfo by reading Field Descriptor Array.
%
% FieldInfo is a 1-by-numFields structure array with the following fields:
%       Name      A string containing the field name 
%       Type      A string containing the field type 
%       ConvFunc  A function handle to convert from DBF to MATLAB type
%       Length    A number equal to the length of the field in bytes

lengthOfLeadingBlock    = 32;
lengthOfDescriptorBlock = 32;
lengthOfTerminator      =  1;
fieldNameOffset         = 16;  % Within table field descriptor
fieldNameLength         = 11;

% Get number of fields.
fseek(fid,8,'bof');
headerLength = fread(fid,1,'uint16');
numFields = (headerLength - lengthOfLeadingBlock - lengthOfTerminator)...
               / lengthOfDescriptorBlock;

% Read field lengths.
fseek(fid,lengthOfLeadingBlock + fieldNameOffset,'bof');
lengths = fread(fid,[1 numFields],'uint8',lengthOfDescriptorBlock - 1);

% Read the field names.
fseek(fid,lengthOfLeadingBlock,'bof');
data = fread(fid,[fieldNameLength numFields],...
             sprintf('%d*uint8=>char',fieldNameLength),...
             lengthOfDescriptorBlock - fieldNameLength);
data(data == 0) = ' '; % Replace nulls with blanks
names = cellstr(data')';

% Read field types.
fseek(fid,lengthOfLeadingBlock + fieldNameLength,'bof');
dbftypes = fread(fid,[numFields 1],'uint8=>char',lengthOfDescriptorBlock - 1);

% Convert DBF field types to MATLAB types.
typeConv = dbftype2matlab(upper(dbftypes));

% Return a struct array.
fieldInfo = cell2struct(...
    [names;  {typeConv.MATLABType}; {typeConv.ConvFunc}; num2cell(lengths)],...
    {'Name', 'Type',                'ConvFunc',          'Length'},1)';

%----------------------------------------------------------------------------
function typeConv = dbftype2matlab(dbftypes)
% Construct struct array with MATLAB types & conversion function handles.

typeLUT = ...
    {'N', 'double', @str2double2cell;...   % DBF numeric
     'F', 'double', @str2double2cell;...   % DBF float
     'C', 'char',   @cellstr;...           % DBF character
     'D', 'char',   @cellstr};             % DBF date

unsupported = struct('MATLABType', 'unsupported', ...
                     'ConvFunc',   @cellstr);

% Unsupported types: Logical,Memo,N/ANameVariable,Binary,General,Picture

numFields = length(dbftypes);
if numFields ~= 0
  typeConv(numFields) = struct('MATLABType',[],'ConvFunc',[]);
end
for k = 1:numFields
    idx = strmatch(dbftypes(k),typeLUT(:,1));
    if ~isempty(idx)
        typeConv(k).MATLABType = typeLUT{idx,2};
        typeConv(k).ConvFunc   = typeLUT{idx,3};
    else
        typeConv(k) = unsupported;
    end
end

%----------------------------------------------------------------------------
function out = str2double2cell(in)
% Translate IN, an M-by-N array of class char, to an M-by-1 column vector
% OUT, of class double.  IN may be blank- or null-padded. If IN(k,:) does
% not represent a valid scalar value, then OUT(k) has value NaN.
if isempty(in)
    out = {[NaN]};
    return
end

% Use sprintf when possible, but fall back to str2double for unusual cases.
fmt = sprintf('%%%df',size(in,2));
[data count] = sscanf(reshape(in',[1 numel(in)]),fmt);
if count == size(in,1)
    out = cell(count,1);
    for k = 1:count
        out{k} = data(k);
    end
else
    out = num2cell(str2double(cellstr(in)));
end