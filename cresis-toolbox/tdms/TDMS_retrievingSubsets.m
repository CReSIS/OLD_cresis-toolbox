function TDMS_retrievingSubsets
%TDMS_retrievingSubsets  Documentation for retrieving subsets of raw data from TDMS_readTDMSFile
%
%   I created this file as a place to document retrieval of subsets of
%   data, instead of making the documentation in the main function rather
%   long.  All property/value pairs would still get passed into the
%   TDMS_readTDMSfile.
%
%   NOTE: FOR MULTIPLE CALLS TO THE SAME FILE
%   It is recommended that you pass the meta_struct from an initial read
%   back into the function.  See examples below.
%
%   There are currently 3 methods for filtering data:
%
%   PROPERTY NOTES:
%   *  indicates input dependence, may be required, depends on input
%   ** indicates an applicable but always optional property
%   All properties are defined 
%
%   1) OBJECT INCLUSION, EXCLUSION, NO RAW DATA
%   =======================================================================
%   Specificy which channels to retrieve or not to retrieve, or whether
%   or not to just get properties and no raw data.  This is implemented via
%   the GET_DATA_OPTION property, in conjunction with OBJECTS_GET and
%   OBJECTS_IGNORE.  In this case, if raw data is retrieved, then all of
%   the raw data available for each non-ignored object is retrieved.
%
%   PROPERTIES: (defined in separate section)
%       - GET_DATA_OPTION
%       - * OBJECTS_GET
%       - * OBJECTS_IGNORE
%       - ** META_STRUCT
%
%   EXAMPLE:
%   a) [~,metaStruct] = TDMS_readTDMSfile(fileName,'GET_DATA_OPTION','getnone')
%
%      %some code ...
%      objStruct = struct;
%      objStruct.groupsKeep = {'myGroup1' 'myGroup2'};
%   
%      data = TDMS_readTDMSfile(fileName,'META_STRUCT',metaStruct,...
%               'GET_DATA_OPTION','getSubset','OBJECTS_GET',objStruct);
%
%
%   2) SINGLE SUBSET OF RAW DATA, FOR ALL NON-IGNORED OBJECTS
%   =======================================================================
%   This method was designed as an extension to the specification above.
%   Use of the SUBSET_GET property, allows additional filtering for those
%   channels that are not excluded.  The same amount of data is read from
%   each channel, and only one subset can be specified. 
%   
%   PROPERTIES: 
%       - SUBSET_GET
%       - ** GET_DATA_OPTION
%       - * OBJECTS_GET
%       - * OBJECTS_IGNORE
%       - ** META_STRUCT
%       - ** SUBSET_IS_LENGTH
%       
%   EXAMPLES:
%   a) [~,metaStruct] = TDMS_readTDMSfile(fileName,'GET_DATA_OPTION','getnone')
%
%   %some code ...
%   
%   data = TDMS_readTDMSfile(metaStruct,'SUBSET_GET',[100 2000],...
%                       'SUBSET_IS_LENGTH',true)
%
%   This reads the same amount of data from all objects with data.  The
%   GET_DATA_OPTION has not been used to limit which objects to retrieve.
%   The SUBSET_IS_LENGTH = true (defined below), indicates that the 2000
%   value indicates the # of samples to read, not the stopping index.
%
%   3) COMPLETE SPECIFICATION OF DATA TO RETRIEVE
%   =======================================================================
%   This was designed to offer even more flexibility.  In this case the
%   GET_DATA_OPTION and SUBSET_GET are IGNORED.  A new variable,
%   GET_INDICES needs to be specified.
%
%   PROPERTIES:
%       - GET_INDICES
%       - ** META_STRUCT
%       - ** SUBSET_IS_LENGTH
%
%   EXAMPLES:
%   a) [~,metaStruct] = TDMS_readTDMSfile(fileName,'GET_DATA_OPTION','getnone')
%
%   %some code
%   
%   getStruct = struct('group',myGroup,'channel',myChannel,...
%       indices,[1 5; 10 5; 20 5; 30 5]);
%
%   data = TDMS_readTDMSfile(metaStruct,'GET_INDICES',getStruct);
%
%   SIDENOTE ON IMAGES:
%   Since tdms does not support 2d writing, image data itself can be 
%   represented in 1d with multiple lines.  Reading multiple contiguous 
%   frames of image data would only require a single read.  If however, for
%   each of those images, a subset of that frame is needed, then data is
%   needed in a non-contiguous fashion. This may be acheivable via 
%   indexing in Matlab.  When memory is a concern, this might be better to
%   do when originally reading in the data. This was the original motivation
%   for developing this read specification.
%
%   =======================================================================
%                         PROPERTIES - DEFINITIONS
%   =======================================================================
%   NOTE: These values get passed in as property/value pairs into the 
%   TDMS_readTDMSfile
%
%   GET_DATA_OPTION = (default 'getAll')
%                       NOTE: This applies specifically to the channels and
%                       groups, not to the data that is in those objects
%               'getall'       - retrieve raw data from all objects
%               'getnone'      - retrieve raw data from no objects
%               'ignoreSubset' - requires population of OBJECTS_IGNORE
%                                retrieves raw data for all objects except
%                                those specified
%               'getSubset'    - requires population of OBJECTS_GET
%                                retrieves raw data only from objects
%                                specified
%
%   OBJECTS_GET     = must be specified for GET_DATA_OPTION = 'getSubset'
%               (structure) with fields:
%               .fullPathsKeep - (cell array of strings) the full object
%                                path name like /'myGroup'/'myChan' (NOTE:
%                                the apostrophes are needed)
%               .groupsKeep    - (cell array of strings) the name of any
%                                 groups to keep -> 'myGroup', all channels
%                                 in that group are retrieved
%   OBJECTS_IGNORE  = must be specified for GET_DATA_OPTION = 'ignoreSubset'
%               (structure with fields:
%               .fullPathsIgnore - (see above), specified paths don't
%                                  retrieve raw data
%               .groupsIgnore    - (see above), raw data not retrieved for 
%                                  any of the channels in the groups mentioned
%
%   SUBSET_GET     = (2 column array), cols =>
%                    [start sample, # of samples] 
%                               or 
%                    [start sample, stop sample] 
%                               for SUBSET_IS_LENGTH = false
%                    The default, [], indicates to retrieve all raw data 
%                    for an object.  WARNING: Only implemented for
%                    decimated data, not interleaved
%
%   SUBSET_IS_LENGTH = (default true), if true, then the 2nd input for
%                       SUBSET_GET, and the 2nd column for the 'indices'
%                       field in GET_INDICES, is treated as the # of
%                       samples to read, if false, it is treated as the
%                       stopping index
%
%   GET_INDICES     = (structure array), length is equal to the number
%                      of channels to retrieve
%               .group   - group name of object to retrieve
%               .channel - channel name of object to retrieve
%               .indices - (rows by 2 columns), first column indicates the 
%                          the first sample to retreive (1 based), 2nd
%                          columns indicates either the # of samples to
%                          retrieve (SUBSET_IS_LENGTH = true) or the last
%                          sample to grab (SUBSET_IS_LENGTH = false)
%                          WARNINGS: 
%                            - Only implemented for decimated data
%                            - indices must be sorted and non-overlapping
%                            i.e. for SUBSET_IS_LENGTH = false,
%                               .indices = [5 10; 15 20] is fine
%                               .indices = [5 10; 10 20] fine ...
%                               NOTE: 10 will only be returned once
%                               .indices = [5 10; 8 15] not ok, (overlapping)
%                               .indices = [15 20; 5 10] not ok, (not sorted)
%   
%   META_STRUCT     = The meta struct from a previous parsing. Passing this
%                     in when reading the same file that generated the meta
%                     struct previously can significantly speed up
%                     processing time.
%
%
%   See Also:
%       TDMS_readTDMSFile