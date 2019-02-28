%                   TDMS_VERSION_INFO
%
%   Current Version: 2.5
%   Date    : 7/28/2012
%   Authors : James Hokanson
%
%   The most notable reason for this change was to make a major bug fix.
%   Major only in the sense that it could fail without the user knowing it
%   (fail silently), giving the illusion of working. 
%
%   In addition I now allow files that were not closed properly to be read,
%   although I only read up to the final segment. This is different from
%   the NI drivers which seem to be able to read parts of the final segment
%   as well. A warning is shown if this happens. At some point I would like
%   to set this as an option (to throw an error or warning) with the
%   default being an error. I would also like to mimic the drivers behavior
%   of parsing out the information in the last segment.
%
%   BUG FIX: POTENTIALLY SILENT BUG *************** 
%   Misinterpreted specification regarding channel data available in a
%   segment, specifically the meaning of "Raw data index" having a value of
%   0. In some situations it is possible that bug could have failed
%   silently (i.e. letting someone thing that the read worked when it
%   didn't) but the liklihood of that having happened is VERY small. This
%   bug was discovered because of a non-silent failure which caused the
%   person's computer to run in a nearly infinite loop.
%   Many thanks to Dedric Xu for reporting the problem.
%   See his post: 22 May 2012
%
%   NOTE 1: This has resulted in NI stating it will update its
%   documentation to be more precise (specific). As of the time of this
%   version being published it has yet to be updated.
%   NOTE 2: Since the specifics on how NI writes the TDMS file are not
%   publically available, I can't specify exactly how the bug would have
%   failed silently. Failure would have resulted in data from channels
%   in certain write cycles being swapped. 
%
%
%   Current Version: 2.4
%   Date    : 10/12/2011
%   Authors : James Hokanson
%   
%       Fixed some bugs and tried to improve documentation. I had been
%       holding off on this as I was hoping to release a few more features
%       but I decided I need to get these fixes out.
%
%   NON-SILENT ERRORS RELATED TO INTERLEAVED DATA
%   ----------------------------------------------
%   BUG FIX: Fix variable naming for reading interleaved data. This one was
%   an error in copy/paste on my part but caused an error if used.
%
%   BUG FIX: Fixed reading of interleaved booleans and timestamps. I wasn't
%   reading the correct number of samples, this also caused an error when
%   run.
%
%   BUG FIX: I didn't update a variable that indicated the numeber of
%   samples read so an error was being needlessly thrown for interleaved
%   data
%
%   SILENT ERRORS: 
%   -------------------------
%   BUG FIX: Fixed reading of timestamps for dates prior to 1904 when written
%   as a channel (not a problem for properties). My guess is almost no one
%   would ever use this but it came up in testing.
%   
%   Current Version: 2.3
%   Date    : 7/14/2011
%   Authors : James Hokanson
%
%       Fixed a few small bugs and tried to further improve documentation.
%       Also allows now for only passing in .tdms_index, for debugging purposes
%
%       BUG FIX: Fixed problem with a check on whether or not all data
%       requested was actually returned to the user. Thanks goes to 
%       Juha Suomalainen for pointing out the problem.
%
%       BUG FIX: Parsing of objects with / characters in their names was 
%       incorrect.  I have generalized the parsing of the objects to allow
%       for ' and / characters in the name.  Thanks goes to 
%       Craig R. Smith for pointing out the problem.
%
%   Current Version: 2.2
%   Date    : 5/17/2011
%   Authors : James Hokanson
%
%       - Added GET_INDICES input with enhanced flexibility for only 
%        retrieving subsets of the data
%
%       - Changed code to allow faster retrieval of subsets. This is
%       accomplished by documenting that start position of each bit of data
%       that belongs to a channel. I had originally been hesitant to do
%       this as this requires more memory, but using this information
%       allows for a significant reduction in the # of required freads and
%       fseeks when multiple channels are available
%
%       - Moved all functions that should not be directly called into a 
%       subfolder, "tdmsSubfunctions" which needs to be added to the path 
%       as well
%
%   Current Version: 2.1
%   Date    : 4/29/2011
%   Authors : James Hokanson
%             assistance by Preston K. Manwaring
%
%       - BUG FIX: The subset parsing was not being handled correctly in
%         Version 2.0, should be fixed in this version
%
%   Version: 2.0
%   Date    : 4/28/2011
%   Authors : James Hokanson
%             assistance by Preston K. Manwaring
%   Updates :
%       - MAJOR UPDATE: Added on the ability to only get a subset of the
%       data from channels specified as a start and a length, i.e. if a
%       channel has 1 million values you could choose to read the data 10x,
%       each time grabbing a new set of 100,000 values for doing whatever
%       processing before moving onto the next chunk
%       
%       This feature is implemented as a option called GET_SUBSET, and
%       applies to all objects that are being retrieved.
%
%       This only works for decimated data currently.
%
%       NOTE: To facilitate quick multiple reads, either make a habit of
%       defragmenting your TDMS files using Labview, or capture the
%       metaStruct output when first processing a file, and pass it back in
%       on subsequent runs -> pass in via optional parameter, META_STRUCT
%
%       - NEW FILE: TDMS_readChannelOrGroup
%       - NEW FILE: TDMS_dataToGroupChanStruct_v4
%       - NEW FILE: TDMS_getStruct
%       
%
%   Version: 1.2
%   Date    : 3/19/2011
%   Author  : James Hokanson
%   Updates :
%       - BUG FIX: conversion of a timestamp property to value was incorrect 
%       I needed a '+' sign instead of a '-' sign
%       Thanks goes to Ed Zechmann for pointing out this bug
%       
%       - Added on function TDMS_dataToGroupChanStruct_v3 which now casts
%       objects AND properties to field names in a struct
%
%   -----------------------------------------------------------------------
%   Version 1.1
%   Author  : James Hokanson
%   Updates : 
%       - MAJOR BUG FIX: fixed the way unicode strings are handled
%
%       Matlab specifies the # of characters to read for unicode strings
%       Labview specifies the # of bytes to read for unicode strings
%
%       Thus one needs to read in the desired # of bytes, and then convert
%       those bytes to characters
%
%       http://www.mathworks.com/matlabcentral/newsreader/view_thread/302145
%
%       --------------------
%       - BUG FIX: fixed case sensitivity on {'GET_DATA_OPTION','getnone'}
%       - BUG FIX: fixed bug with ignoring data retrieval on certain
%           objects
%       - BUG FIX: fixed skipping bug on timestamp data
%       - all structure fields no longer necessary for 
%          ignoreSubset or getSubset GET_DATA_OPTIONs
%       - added function dataToGroupChanStruct_v2
%       - added function TDMS_exampleFunctionCalls
%
%   -----------------------------------------------------------------------
%   Version 1.0
%   Author: 
%       James Hokanson
%       University of Pittsburgh
%       Graduate Student Researcher
%   Date : January 13, 2011