%The following is some example code of calls that can be made

MY_TEST_FILE_PATH = '';
testDirectory = '';
saveDirectory = '';


%% 1) Batch functionality - get all TDMS files in a directory and save the
%processed versions to another directory
%==========================================================================
fileStruct    = dir(fullfile(testDirectory,'*.tdms'));
fileNames     = {fileStruct.name};

nFiles       = length(fileNames);
successFlags = true(1,nFiles);
meCA         = cell(1,nFiles);
for iFile = 1:nFiles
   try
       temp          = TDMS_readTDMSFile(fullfile(testDirectory,fileNames{iFile}));
       [output,nav]  = TDMS_dataToGroupChanStruct_v2(tempOutput); %#ok<*NASGU,*ASGLU>
       [~,nameClean] = fileparts(fileNames{iFile});
       save(fullfile(saveDirectory,nameClean),'output','nav')
   catch ME
       meCA{iFile} = ME;
       successFlags(iFile) = false;
   end
end

allPassed = all(successFlags);

%% 2) Requesting of limited data - retrieve no data (does return properties)
%==========================================================================
filePath = MY_TEST_FILE_PATH;

params = {'GET_DATA_OPTION','getnone'};
tempOutput = TDMS_readTDMSFile(filePath,params{:});
[output,nav] = TDMS_dataToGroupChanStruct_v2(tempOutput);
%NOTE: In this case nav might be the most interesting

%% 3) Requesting of limited data - get a subset of data
%==========================================================================
filePath = MY_TEST_FILE_PATH;

params = {'GET_DATA_OPTION','getSubset','OBJECTS_GET',...
    struct('groupsKeep',{{'Marker 0' 'Marker 1'}})};
%NOTE: Make sure to double the {} so as to not make the structure be a
%structure array
tempOutput = TDMS_readTDMSFile(filePath,params{:});
[output,nav] = TDMS_dataToGroupChanStruct_v2(tempOutput);
%NOTE: In this case nav might be the most interesting

%% 4) Requesting of limited data - ignore a subset of data
%==========================================================================
filePath = MY_TEST_FILE_PATH;

params = {'GET_DATA_OPTION','ignoreSubset','OBJECTS_IGNORE',...
    struct('groupsIgnore',{{'Marker 0' 'Marker 1'}})};
%NOTE: Make sure to double the {} so as to not make the structure be a
%structure array
tempOutput = TDMS_readTDMSFile(filePath,params{:});
output = TDMS_dataToGroupChanStruct_v1(tempOutput);

%% 5) Use of other input paramters
%==========================================================================
filePath = MY_TEST_FILE_PATH;

params = {...
    'UTC_DIFF'  -6  ...%Refer timestamps to U.S. Central Time (UTC - 6)
    'USE_INDEX' false ... %index may be corrupted, use tdms file, not tdms_index
    'MAX_NUM_OBJECTS' 1100 ... %We have about 1000 objects, +100 for padding
    'DATE_STR_FORMAT' 'dd-mmm-yyyy' ... %Only return data info for properties
    %NOTE: timestamp data is returned in such a format that datestr() works
    %on it
    };
tempOutput = TDMS_readTDMSFile(filePath,params{:});
output = TDMS_dataToGroupChanStruct_v1(tempOutput);